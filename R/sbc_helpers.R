#' @export sbch_dx_plot
sbch_dx_plot <- function(x, scale_x = TRUE) {
  x %>%
    dplyr::mutate(p50_err = (median-simulated_value)) %>%
    tidyr::pivot_longer(cols = c(p50_err, ess_bulk, ess_tail, rhat)) %>%
    ggplot2::ggplot(aes(x = value, y = variable, fill = setting_id)) +
    ggplot2::geom_boxplot(aes(color = setting_id), alpha = 0.33,
                          outlier.color = NULL) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::facet_grid(cols = vars(name), scales = "free") -> p

  if(scale_x)p <- p + ggplot2::scale_x_continuous(trans = ggallin::ssqrt_trans)

  p
}

#' Construct `cmdstanr` backend for SBC
#'
#' @param model_path Stan model file location
#' @param exe_path Compiled model executable location
#'
#' @return Object of class `SBC_backend_cmdstan_sample` with an added `model_name` attribute
#' @export
sbch_cmdstan_backend <- function(model_path, exe_path) {
  SBC::SBC_backend_cmdstan_sample(cmdstanr::cmdstan_model(stan_file = model_path, exe_file = exe_path)) -> backend_obj
  attr(backend_obj, "model_name") <- basename(exe_path)
  backend_obj
}

#' Construct `brms` backend for SBC with more options
#'
#' The most important difference with the vanilla SBC constructor is that this will let you store the compiled model to avoid recompilation when re-running a `targets` pipeline.
#'
#' @param ... Arguments passed to `brm`
#' @param template_data Sample of data `brms` will work with
#'
#' If a function is passed as `formula` or `prior` arguments via `...`, this
#' helper will send the result of applying the function on `template_data` to
#' the backend instead.
#'
#' @export
sbch_brms_backend <- function(..., template_data) {
  dots <- list(...)
  save_stanvars <- dots$stanvars
  if(is.function(dots$formula))dots$formula <- dots$formula(template_data)
  if(is.function(dots$prior))dots$prior <- dots$prior(template_data)
  if(is.function(dots$stanvars)) {
    dots$stanvars <- dots$stanvars(template_data)
  }

  backend_obj <- do.call(SBC::SBC_backend_brms,
                         c(dots, list(template_data = template_data)))

  backend_obj$args$stanvars <- save_stanvars
  backend_obj
}

#' Construct SBC generator from a function and arguments list
#'
#' @param gen_fun The generator function
#' @param arg_list List of named arguments that will be passed to the generator
#' @param ... Other parameters passed to the generator
#'
#' @return An object of class `SBC_generator_function`
#' @export
sbch_generator <- function(gen_fun, arg_list, ...) {
  stopifnot({is.function(gen_fun);is.list(arg_list)})
  do.call(SBC::SBC_generator_function, c(gen_fun, arg_list, ...))
}

#' Adds more flexibility to the SBC generator
#'
#' @param generator SBC generator function
#' @param n_sims Number of datasets to simulate (usual behaviour)
#' @param n_reps Number of times to repeat a single simulated dataset (for other MCMC assessment purposes)
#' @param mvbrms `FALSE` will use the default `SBC::generate_datasets` function. By default, we use `sbch_generate_semdata`.
#' @param makena A character vector with the names of variables that will be replaced with `NA_real_`.
#'
#' @return Datasets in the format SBC expects
#' @export
sbch_generate <- function(generator, ..., n_sims = NULL, n_reps = NULL,
                          mvbrms = TRUE) {
  if(is.null(n_sims)&is.null(n_reps))stop("Must specify n_sims OR n_reps")
  if(!is.null(n_sims)&!is.null(n_reps))stop("Set only ONE of n_sims or n_reps")
  generate_fun <- sbch_generate_semdata
  if(!mvbrms){
    warning("mvbrms = FALSE untested in sbch_generate")
    generate_fun <- SBC::generate_datasets
  }

  datasets <- generate_fun(generator, c(n_sims, n_reps))

  if(!is.null(n_reps)) {
    draws_attr <- attributes(datasets$variables)
    datasets$variables <- datasets$variables[rep(1, n_reps),]
    attr(datasets$variables, "nchains") <- draws_attr$nchains
    attr(datasets$variables, "class") <- draws_attr$class
    attr(datasets$variables, "dimnames")$draw <- as.character(1:n_reps)
    datasets$generated <- datasets$generated[rep(1, n_reps)]
  }

  datasets
}

#' A concise execution of the SBC workflow
#'
#' @param arg_row a tibble row with colum names that correspond to arguments in the generator function
#' @param n_reps,n_sims the amount of
#' @param generator some function that produces datasets based on the values in `arg_row`
#' @param backend model sampler made by some `SBC_backend_*`function
#' @param keep_fits passed to `compute_SBC`, defaults to `TRUE`
#' @param keep_data append the datasets used to the `compute_SBC` result, defaults to `FALSE`
#'
#' @export
sbch_run <- function(arg_row, ..., n_reps = NULL, n_sims = NULL,
                     generator, backend, keep_fits = TRUE,
                     keep_data = FALSE) {
  unpacked_args <- purrr::map(arg_list, unlist)
  sbc_gen <- sbch_generator(generator, unpacked_args)
  datasets <- sbch_generate(sbc_gen, n_reps, n_sims)

  sbc_obj <- SBC::compute_SBC(datasets, backend, keep_fits, cache_mode = "none")
  if(keep_data)sbc_obj$datasets <- datasets

  sbc_obj
}

#' Modified `SBC::generate_datasets` to make SEM data from `brms` models
#'
#' @param generator hmm
#' @param n_sims Number of simulated datasets to produce
#'
#' @return Object of class `SBC_datasets` (contains generated data and true parameter values)
#' @export
#'
#' @examples print("nope")
sbch_generate_semdata <- (\(){
  gd_tmp <- SBC:::generate_datasets.SBC_generator_brms
  body(gd_tmp)[[10]][[4]][[3]][[4]][[3]][[2]] <- substitute({
    generated <- brmsh_pp_iter(prior_fit_brms)
    break
  })
  body(gd_tmp)[[11]] <- substitute({
    nonparameter_names <- "(\\.\\.FILLERforSAMPLING)|(^.+LVi[0-9]+$)"
    drop_generated <- grepl(nonparameter_names, names(generated[[1]]))
    lv_names <- colnames(generated[[1]][,!drop_generated,drop = FALSE])
    gen_n <- nrow(generated[[1]])
    draws <- posterior::as_draws_matrix(prior_fit_brms$fit)
    draws <- draws[,!grepl("^(sigma_FILLERforSAMPLING)|(Intercept_.+)$",
                           colnames(draws))]
    for(i in 1:n_sims) {
      if (generator$generate_lp) {
        ll <- log_lik(prior_fit_brms, newdata = generated[[i]],
                      draw_ids = i, cores = 1)
        log_likelihoods[i] <- sum(ll)
      }
      sim_values <-unlist(generated[[i]][,!drop_generated])
      names(sim_values) <- paste0(
        rep(paste0("Ymi_", lv_names), each = gen_n),
        "[", 1:gen_n, "]")
      draws[i, which_order(names(sim_values), colnames(draws))] <- sim_values
    }
  })
  gd_tmp
})()

# dont export, called from bdlvm
compute_SBC_stanvars <- (\(){
  func <- SBC::compute_SBC
  body(func)[[10]] <- substitute(
    results_raw <- future.apply::future_lapply(
      vars_and_generated_list,
      lunafuns:::stanvars_SBC_single, backend = backend, cores = cores_per_fit,
      keep_fit = keep_fits, thin_ranks = thin_ranks, ensure_num_ranks_divisor = ensure_num_ranks_divisor,
      dquants = dquants, future.seed = TRUE, future.globals = future.globals,
      future.chunk.size = chunk_size))
  func
})()

stanvars_SBC_single <- (\(){
  func <- SBC:::compute_SBC_single
  body(func)[[4]][[3]][[2]][[2]][[3]][[2]][[2]] <- substitute({
    backend$args$stanvars <- backend$args$stanvars(generated)
    fit <- SBC_fit(backend, generated, cores = cores)
  })
  func
})()
