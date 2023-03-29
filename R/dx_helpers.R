#' @export dxh_prepare_plot_data
dxh_prepare_plot_data <- function(data) {
  list(
    backend_diagnostics =
      data$backend_diagnostics %>%
      dplyr::group_by(sim_desc) %>%
      dplyr::summarise(failed_chain_pct = ceiling(100 - 25*(sum(successful_chains)/dplyr::n())),
                       max_chain_time = max(max_chain_time, na.rm = TRUE),
                       max_divergent = max(n_divergent, na.rm = TRUE),
                       max_treedepth = max(n_max_treedepth, na.rm = TRUE),
                       min_bfmi = min(min_bfmi, na.rm = TRUE),
                       max_rejects = max(n_rejects, na.rm = TRUE)) %>%
      dplyr::mutate(dplyr::across(-c(sim_desc), \(x)ifelse(is.infinite(x), NA, x)),
                    failed_chain_pct = ifelse(failed_chain_pct == 0, NA, failed_chain_pct),
                    variable = sort(data$fit_diagnostics$variable)[1]),
    fit_diagnostics =
      data$fit_diagnostics %>%
      dplyr::group_by(sim_desc, variable) %>%
      dplyr::summarise(rmse = sqrt(mean((mean - simulated_value)**2, na.rm = TRUE)),
                       rhat = mean(rhat, na.rm = TRUE),
                       ess_bulk = mean(ess_bulk, na.rm = TRUE),
                       ess_tail = mean(ess_tail, na.rm = TRUE))
  )
}

dxh_get_num_chains <- function(x, ...) {
  UseMethod("dxh_get_num_chains", x)
}

#' @method dxh_get_num_chains brmsfit
#' @export
dxh_get_num_chains.brmsfit <- function(x, ...)x$fit@sim$chains

#' @method dxh_get_num_chains blavaan
#' @export
dxh_get_num_chains.blavaan <- function(x, ...)x@external$mcmcout@sim$chains

#' @method dxh_get_num_chains lavaan
#' @export
dxh_get_num_chains.lavaan <- function(x, ...) {
  conv <- lavaan::inspect(x, what = "converged")
  if(conv)return(NA)
  0
}

#' @export dxh_fits_summary
dxh_fits_summary <- function(sim_desc, x) {
  stats <- dplyr::bind_cols(
    sim_desc = sim_desc,
    x$stats
  )
  backend <- tibble::tibble(
    sim_desc = sim_desc,
    sim_id = seq_along(x$fits),
    successful_chains = purrr::map_int(x$fits, \(x) {
      if(is.null(x))return(0)
      dxh_get_num_chains(x)
    })
  ) %>%
    dplyr::left_join(x$backend_diagnostics, by = "sim_id")

  x <- list(backend_diagnostics = backend, fit_diagnostics = stats)
  x
}

#' @export dxh_prep_stats
dxh_prep_stats <- function(x) {
  if(!"sim_desc" %in% names(x)) {
    warning("sim_desc column not found, assuming all stats belongs to same sim")
    x$sim_desc <- ""
  }
  x %>%
    dplyr::transmute(
      sim_desc, variable,
      bias = mean - simulated_value,
      rhat, ess_bulk, ess_tail) %>%
    tidyr::pivot_longer(!c(sim_desc, variable)) %>%
    dplyr::group_by(sim_desc, variable, name) %>%
    dplyr::filter(!variable %in% c("lp__", "lprior"))
}

dxh_rescale <- function(x, breaks) {
  u <- max(breaks)
  l <- min(breaks)
  x <- pmin(u, x, na.rm = FALSE)
  x <- pmax(l, x, na.rm = FALSE)
  (x-l)/(u-l)
}

dxh_scale_rhat <- scales::trans_new(
  name = "pseudolog10",
  transform = function (x)asinh((x-1.05)*50)/log(10),
  inverse = function (y)  sinh(y * log(10))/50+1.05,
  breaks = function(...) c(1.01, 1.05, 1.2)
)
