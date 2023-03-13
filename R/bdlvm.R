#' @name bdlvm
#' @rdname bdlvm
#'
#' @title `bdlvm_*` functions for streamlined model generation/fitting
#'
#' @description
#' These functions have been designed with the following scenario in mind:
#' - Someone has a specific SEM in mind they wish to test
#' - They may want to generate data from a variety of parameter values
#' - They may want to fit to that data using a variety of model specifications
#'
#' Therefore, it makes sense to break up each of those steps into different functions, so that e.g. a model needs to be specified only once, while generation and fitting are each specified separately. In the general case, then one can have a total of `S*D*F` sets of results, where `S` is the number of SEMs, `D` is the number of datasets and `F` is the number of fitting specifications.
#'
#' See details for the workflow in terms of specific functions. Note these package is currently designed to work solely on top of `brms` and `SBC`.
#'
#' @details
#' The specific workflow is as follows:
#' 1. Provide a description of the SEM "core" (i.e. exogenous and latent variables only) via a `brms` formula. Any variables tagged with the `mi()` function are assumed to be latent variables.
#' 2. Feed the "base formula" into `bdlvm_model_setup` to produce a model object that will have a "full formula" and template dataset that has added variables corresponding to the items for each latent variable.
#' 3. Feed the model object to `bdlvm_make_data` together with the number of desired datasets to simulate the data (this is a SLOW STEP as it uses `brms::posterior_predict` iteratively to populate the entire graph).
#' 4. Feed the dataset and model objects to `bdlvm_compute` in order to fit the model to each dataset via a modified version of `SBC::compute_SBC`. The same generative priors will be used for fitting by default, but see the arguments for ways to specify different priors. EXPERIMENTAL: the backend object SBC requires will be generated at this point in order to guarantee compilation of the latest model specification. YES I LEARNED THIS IN A VERY PAINFUL MANNER.
#'
NULL

#' @describeIn bdlvm Make a MODEL OBJECT out of a BASE FORMULA
#'
#' @param formula `brmsformula` specifying the "core" model.
#' @param data If there is no exogenous variables, an integer suffice to indicate you want a dataframe with the amount of rows as part of the output. Otherwise, a dataframe with the values of exogenous variables is expected and its number of rows will be used for the template data.
#' @param M The number of items to be added on each latent variable. Currently a single number i.e. all the latent variables will get the same amount of items.
#'
#' @export
#'
#' @examples
#' suppressMessages(library(brms))
#' bdlvm_model_setup(bf(y1|mi()~1), data = 50, M = 3)
bdlvm_model_setup <- function(formula, data, M) {
  # This is a separate step from data generation in case one wishes to
  # add actual mi() variables to the model,

  # Identify LVs (formula is not checked for consistent use of mi() on RHS)
  f_list <- brmsh_get_formulas(formula)
  lv_names <- get_lhs_vars(f_list)[get_lhs_mi(f_list)]
  n_lv <- length(lv_names)

  # Dataframe size
  if(is.data.frame(data)){N<-nrow(data)}else{N<-data}

  out <- list()
  # Construct template dataset
  out$data <- as.data.frame(matrix(0, nrow = N, ncol = n_lv*(1+M)))
  colnames(out$data) <- rep(lv_names,
                            each = M+1) %p0% c("", rep("LVi", M) %p0% 1:M)
  i_names <- colnames(out$data)[-(1+(0:(n_lv-1))*(M+1))]

  if(is.data.frame(data)) {
    out$data <- cbind(data[,setdiff(colnames(data), lv_names)], out$data)
  }

  # Extend formula with LV items
  out$formula <- formula
  for(lv in lv_names) {
    for(m in 1:M) {
      out$formula <- out$formula +
        brms::bf(paste0(lv, "LVi", m, "| mi() ~ mi(", lv, ")"))
    }
  }

  # Make prior dataframe
  out$prior <- brms::get_prior(out$formula + set_rescor(F), out$data)
  out$prior <- prepare_prior_df(out$prior, lv_names, i_names)

  out$lv_names <- lv_names
  out$i_names <- i_names
  out
}

prepare_prior_df <- function(df, idL, idI) {
  df <- with(df, df[resp!="" & !(class == "b" & coef == ""), ])
  df[, "order"] <- 1:nrow(df)

  with(df, {
    item_inter <- intersect(which_order(idI, resp),
                            which(class == "Intercept"))
    item_sigma <- intersect(which_order(idI, resp),
                            which(class == "sigma"))
    item_loadf <- intersect(which_order(idI, resp),
                            which(class == "b"))
    lv_inter <- intersect(which_order(idL, resp),
                          which(class == "Intercept" & dpar == "" & nlpar == ""))
    lv_sigma <- intersect(which_order(idL, resp),
                          which(class == "sigma" & nlpar == ""))
    lv_betas <- {
      tmp <- df[which_order(idL, resp),]
      tmp[intersect(which_order("mi" %p0% idL, tmp$coef),
                    which(tmp$class == "b" & tmp$dpar == "" & tmp$nlpar == "")),
          "order"]
    }
    lv_dist_inter <- intersect(which_order(idL, resp),
                               which(class == "Intercept" & dpar != "" & nlpar == ""))
    lv_dist_betas <- {
      tmp <- df[which_order(idL, resp),]
      tmp[intersect(which_order("mi" %p0% idL, tmp$coef),
                    which(tmp$class == "b" & tmp$dpar != "" & tmp$nlpar == "")),
          "order"]
    }

    df[item_inter, "sem_role"] <- sem_prior_type[1]
    df[item_sigma, "sem_role"] <- sem_prior_type[2]
    df[item_loadf, "sem_role"] <- sem_prior_type[3]
    df[lv_inter, "sem_role"] <- sem_prior_type[4]
    df[lv_sigma, "sem_role"] <- sem_prior_type[5]
    df[lv_betas, "sem_role"] <- sem_prior_type[6]
    df[lv_dist_inter, "sem_role"] <- sem_prior_type[7]
    df[lv_dist_betas, "sem_role"] <- sem_prior_type[8]

    lv_other <- intersect(which_order(idL, resp),
                          which(is.na(df$sem_role) & dpar == ""))
    lv_d_other <- intersect(which_order(idL, resp),
                            which(is.na(df$sem_role) & dpar != ""))
    df[lv_other, "sem_role"] <- sem_prior_type[9]
    df[lv_d_other, "sem_role"] <- sem_prior_type[10]

    df <- rbind(df[item_inter,], df[item_sigma,], df[item_loadf,],
                df[lv_inter,], df[lv_sigma,], df[lv_betas,],
                df[lv_dist_inter,], df[lv_dist_betas,],
                df[lv_other,], df[lv_d_other,],
                df[is.na(df$sem_role),])

    attr(df, "sem_role") <- df$sem_role
    rownames(df) <- NULL
    df$sem_role <- NULL
    df$order <- NULL

    df
  })
}

#' @export bdlvm_make_data
bdlvm_make_data <- function(model_obj, ...) {
  bdlvm_generate_datasets(bdlvm_generator_setup(model_obj, ...), ...)
}

bdlvm_generator_setup <- function(model_obj, priors = NULL, ...) {

  sem_pars <- attr(model_obj$prior, "sem_role")
  if(is.null(sem_pars)|!any(sem_pars%in%sem_prior_type))stop(
    'No SEM parameters found in attr(., "sem_role")')

  for(type in sem_prior_type) {
    model_obj$prior <- default_prior(model_obj$prior, priors, type)
  }

  filler_row <- brms::set_prior("normal(0, 1.67)", class = "sigma",
                                lb = "0", ub = "", resp = "FILLERforSAMPLING")
  filler_row$source <- "bdlvm_generator_required"
  model_obj$prior <- rbind(filler_row, model_obj$prior)
  attr(model_obj$prior, "sem_role") <- c(NA, sem_pars)
  class(model_obj$prior) <- c("brmsprior", "data.frame")

  model_obj$formula <- model_obj$formula +
    brms::bf(..FILLERforSAMPLING ~ 0) +
    brms::set_rescor(FALSE)
  model_obj$data <- cbind(model_obj$data,
                          ..FILLERforSAMPLING = rnorm(nrow(model_obj$data)))
  model_obj
}

default_prior <- function(df, priors, type) {

  sem_role <- attr(df, "sem_role")

  if(!any(sem_role == type)) {
    if(!is.null(priors[[type]])) {
      warning("Passed priors for "%p0%type%p0%
                " parameters, but not found in model.")
    }
    return(df)
  }

  df[sem_role == type, "prior"] <- {
    if(is.null(priors[[type]])) {
      df[sem_role == type, "source"] <- "bdlvm_generator_default"
      if(type %in% sem_prior_type[9:10])stop("No defaults for "%p0%type%p0%
                                               " parameters.")
      if(type %in% sem_prior_type[c(1, 4, 7)]){"constant(0)"}else{"constant(1)"}
    } else {
      n_pars <- nrow(df[sem_role == type,])
      if(!length(priors[[type]]) %in% c(1, n_pars)) {
        stop("Priors for "%p0%type%p0%
               " parameters must be length 1 or match "%p0%n_pars)
      }
      df[sem_role == type, "source"] <- "bdlvm_generator_user"
      priors[[type]]
    }
  }

  df
}

bdlvm_generate_datasets <- function(model_obj, ...) {
  generator_args <- list(formula = model_obj$formula,
                         data = model_obj$data,
                         prior = model_obj$prior,
                         thin = 50, warmup = 10000, refresh = 5000)

  dots <- list(...)

  generator <- list(do.call(SBC::SBC_generator_brms, generator_args))
  datasets <- do.call(sbch_generate, c(generator, dots))

  if(is.null(dots$lv_names))dots$lv_names <- model_obj$lv_names
  lv_names <- dots$lv_names

  generated <- purrr::list_transpose(lapply(
    datasets$generated, \(x) {
      lv_values <- x[,which(colnames(x) %in% lv_names), drop = FALSE]
      x[,which(colnames(x) %in% lv_names)] <- NA_real_
      list(data = x, lv_values = lv_values)}
    ))

  datasets$generated <- generated$data

  list(datasets = datasets, lv_values = generated$lv_values, prior = model_obj$prior)
}

#' Title
#'
#' Bottom text
#'
#' @details
#'
#' Names in the list passed via the `priors` argument can refer to SEM-specific parameters by being one of
#' `item_intercept`, `item_sigma`, `item_loading`, `lv_intercept`, `lv_sigma`, `lv_to_lv_coef`,
#' `lv_d_intercept`, `lv_to_lv_d_coef`, `lv_other`, `lv_d_other`
#'
#' One can pass a `sem_global` argument to assign a prior on all non-intercept
#' parameters that were not specifically named. Additionally, the name `non_sem`
#' can be used to refer to any other parameters that are not part of the basic SEM
#' structure; this include exogenous variables and `I()` transformations of latent variables.
#'
#' @param model_obj An object made with `bdlvm_model_setup`
#' @param data_obj An object made with `bdlvm_generate_datasets`
#' @param ... Arguments passed to `brms::brm`
#' @param priors Named list passed to `make_fit_prior`, see details for accepted names
#' @param cache_mode Disabling SBC cache by default
#' @param cache_location Disabling SBC cache by default
#'
#' @return Just a regular `SBC_results`-class object.
#' @export bdlvm_compute
#'
#' @examples
#' suppressMessages({
#'   library(brms)
#'   library(lunafuns)
#' })
#'
#' # general model setup
#' base_formula <- bf(y1|mi()~1) +
#'   bf(y21|mi()~mi(y1)) + bf(y22|mi()~mi(y1)) +
#'   bf(y3|mi()~mi(y21)+mi(y22))
#'
#' model_object <- bdlvm_model_setup(base_formula, data = 50, M = 3)
#'
#' # data generation setup
#' dataset_object <- bdlvm_make_data(model_object, n_sims = 2)
#'
#' # fit (with a potentially different prior)
#' fits_object <- bdlvm_compute(model_object, dataset_object,
#'                              priors = list(
#'                                global = "normal(0,5)"
#'                              ))
bdlvm_compute <- function(model_obj, data_obj, ..., debug_mode = FALSE,
                          priors = NULL, var_constraint = FALSE,
                          cache_mode = "none", cache_location = NULL) {

  backend_obj <- make_backend(model_obj, data_obj$prior, ...,
                              priors = priors, var_constraint = var_constraint)

  if(debug_mode)return(list(data = data_obj, backend = backend_obj))

  compute_SBC_stanvars(
    datasets = data_obj$datasets,
    backend = backend_obj,
    cache_mode = cache_mode,
    cache_location = cache_location)
}
make_backend <- function(model_obj, prior_obj, ...,
                         backend_args = NULL, stanvar_func = NULL) {

  fit_prior <- make_fit_prior(model_obj, prior_obj, ...)
  if(is.null(stanvar_func))stanvar_func <- stanvar_func_datavar(model_obj)

  backend_args <- as.list(backend_args)
  if(is.null(backend_args$chains))backend_args$chains <- 4
  if(is.null(backend_args$backend))backend_args$backend <- "cmdstanr"

  backend_args <- c(list(formula = model_obj$formula + brms::set_rescor(FALSE),
                         template_data = model_obj$data, prior = fit_prior,
                         stanvars = stanvar_func),
                    backend_args)

  do.call(sbch_brms_backend, backend_args)
}
make_fit_prior <- function(model_obj, prior_obj, ...,
                           priors = NULL, var_constraint = TRUE,
                           lf1_fixto1 = TRUE, as_is = FALSE) {

  if(as_is)return(prior_obj)

  out_prior <- prior_obj
  out_prior$sem_role <- attributes(out_prior)$sem_role
  out_prior <- out_prior[out_prior$resp != "FILLERforSAMPLING",]

  if(lf1_fixto1) {
    lf1_idx <- grepl("LVi1$", out_prior$resp) & out_prior$sem_role == "item_loading"
    out_prior[lf1_idx, "prior"] <- "constant(1)"
    out_prior[lf1_idx, "source"] <- "bdlvm_compute_fix_item1"
  }

  # other values can be added here when e.g. we add the option to fix variances
  is_fixed <- out_prior$source == "bdlvm_compute_fix_item1"

  if("item_sigma" %in% names(priors) & var_constraint) {
    stop("Must set var_constraint = FALSE to place priors on item_sigma")
  }

  if("sem_global" %in% names(priors)) {

    intercept_or_assigned <- out_prior$sem_role %in% names(priors) |
      out_prior$sem_role %in% sem_prior_type[c(1,4,7)]
    free_pars <- !is_fixed & !intercept_or_assigned

    if(all(intercept_or_assigned)) {
      stop("No non-intercept SEM parameters left for global assignment")
    }
    if(!length(priors$sem_global) %in% c(1, sum(free_pars))) {
      stop("sem_global must be length 1 or "%p0%sum(free_pars))
    }

    out_prior[free_pars, "prior"] <- priors$sem_global
    out_prior[free_pars, "source"] <- "bdlvm_compute_global"
  }

  for(type in sem_prior_type) {
    if(is.null(priors[[type]]))next

    type_rows <- out_prior$sem_role == type & !is_fixed
    if(!length(priors[[type]]) %in% c(1, sum(type_rows))) {
      stop(type%p0%" prior must have length 1 or equal to parameters of its type")
    }

    out_prior[type_rows, "prior"] <- dplyr::coalesce(priors[[type]],
                                                     out_prior[type_rows, "prior"])
    out_prior[type_rows, "source"] <- ifelse(!is.na(priors[[type]]),
                                             "bdlvm_compute_user",
                                             out_prior[type_rows, "source"])
  }

  if(var_constraint) {
    lv_names <- model_obj$lv_names
    i_names <- model_obj$i_names
    M <- length(i_names)/length(lv_names)

    isigma_idx <- intersect(which_order(i_names, out_prior$resp),
                            which(out_prior$sem_role == "item_sigma"))

    out_prior[isigma_idx, "prior"] <- "constant(sqrt(datavar_" %p0%
      i_names %p0% " - (bsp_" %p0% out_prior[isigma_idx, "resp"] %p0%
      "[1]^2)*variance(Ymi_" %p0% rep(lv_names, each=M) %p0% ")))"
    out_prior[isigma_idx, "source"] <- "bdlvm_compute_var_constraint"
  }

  if("non_sem" %in% names(priors)) {
    non_sem <- is.na(out_prior$sem_role)
    if(!length(priors$non_sem) %in% c(1, sum(non_sem))) {
      stop("non_sem must be length 1 or "%p0%sum(non_sem))
    }

    out_prior[non_sem, "prior"] <- priors$non_sem
    out_prior[non_sem, "source"] <- "bdlvm_compute_user"
  }

  out_prior$sem_role <- NULL
  out_prior[!(out_prior$class == "Intercept" &
                out_prior$source == "bdlvm_generator_default"), ]
}
stanvar_func_datavar <- function(model_obj) {
  function(data) {
    var_data <- brms::stanvar(block = "functions", scode = "// Empty stanvars")
    for(i in model_obj$i_names) {
      var_data <- var_data +
        brms::stanvar(x = var(data[,i]),
                      name = paste0("datavar_", i),
                      block = "data")
    }
    var_data
  }
}

#' @export bdlvm_lavaan
bdlvm_lavaan <- function(formula, model_obj, data_obj, ...) {
  fit_obj <- \(x)lavaan::lavaan(formula, data = x, ...)
  data_obj$datasets$generated <- purrr::map(data_obj$datasets$generated,
                                            \(x)x[,!names(x) %in% model_obj$lv_names])
  bdlvm_fit(fit_obj, data_obj$datasets, model_obj)
}

#' @export bdlvm_blavaan
bdlvm_blavaan <- function(formula, model_obj, data_obj, ...,
                          n.chains = 4, bcontrol = list(refresh = 0)) {
  fit_obj <- \(x)blavaan::blavaan(formula, data = x,
                                  n.chains = n.chains, bcontrol = bcontrol, ...)
  data_obj$datasets$generated <- purrr::map(data_obj$datasets$generated,
                                            \(x)x[,!names(x) %in% model_obj$lv_names])
  bdlvm_fit(fit_obj, data_obj$datasets, model_obj)
}

#' @export bdlvm_fit
bdlvm_fit <- function(fit_obj, data_obj, ...) {
  fits <- purrr::map(data_obj$generated, fit_obj, .progress = TRUE)

  stats <- dplyr::full_join(
    tibble::tibble(
      sim_id = rep(seq_along(fits), each = ncol(data_obj$variables)),
      variable = rep(colnames(data_obj$variables), times = length(fits)),
      simulated_value = as.vector(t(data_obj$variables))
    ),
    bdlvm_get_stats(fits, ...),
    by = c("sim_id", "variable"))

  backend_diagnostics <- bdlvm_get_backdiag(fits)

  list(stats = stats, fits = fits, backend_diagnostics = backend_diagnostics)
}

bdlvm_get_backdiag <- function(x, ...) {
  UseMethod("bdlvm_get_backdiag", x)
}

#' @method bdlvm_get_backdiag list
#' @export
bdlvm_get_backdiag.list <- function(x, ...) {
  x <- purrr::map(x, bdlvm_get_backdiag, ...)
  dplyr::bind_rows(purrr::map2(x, seq_along(x), \(x, y)dplyr::bind_cols(sim_id = y, x)))
}

#' @method bdlvm_get_backdiag blavaan
#' @export
bdlvm_get_backdiag.blavaan <- function(x, ...) {
  fit <- x@external$mcmcout
  data.frame(max_chain_time = max(rowSums(rstan::get_elapsed_time(fit))),
             n_divergent = rstan::get_num_divergent(fit), n_max_treedepth = rstan::get_num_max_treedepth(fit),
             min_bfmi = min(rstan::get_bfmi(fit)))
}

#' @method bdlvm_get_backdiag lavaan
#' @export
bdlvm_get_backdiag.lavaan <- function(x, ...)lavaan::inspect(x, what = "timing")$total

bdlvm_get_stats <- function(x, ...) {
  UseMethod("bdlvm_get_stats", x)
}

#' @method bdlvm_get_stats list
#' @export
bdlvm_get_stats.list <- function(x, ...) {
  x <- purrr::map(x, bdlvm_get_stats, ...)
  dplyr::bind_rows(purrr::map2(x, seq_along(x), \(x, y)dplyr::bind_cols(sim_id = y, x)))
}

#' @method bdlvm_get_stats blavaan
#' @export
bdlvm_get_stats.blavaan <- function(fit_obj, model_obj) {
  lv_names <- model_obj$lv_names
  i_names <- model_obj$i_names

  lavaan::partable(fit_obj)[,c("lhs", "op", "rhs", "pxnames", "mat", "free", "est")] %>%
    dplyr::group_by(mat) %>%
    dplyr::transmute(lhs, op, rhs,
                     est = ifelse(free == 0, est, NA_real_),
                     variable = pxnames, pardim = dplyr::n()) %>%
    dplyr::mutate(variable = ifelse(pardim!=1, variable, gsub(".{3}$", "", variable))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mat, -pardim) %>%
    dplyr::left_join(summary(posterior::as_draws_rvars(fit_obj@external$mcmcout)),
                     by = "variable") %>%
    dplyr::transmute(variable = dplyr::case_when(
      op == "=~" & lhs %in% lv_names ~ "bsp_" %p0% rhs %p0% "_mi" %p0% lhs,
      op == "~1" & lhs %in% c(lv_names, i_names) ~ "b_" %p0% lhs %p0% "_Intercept",
      op == "~~" & lhs %in% c(lv_names, i_names) & lhs == rhs ~ "sigma_" %p0% lhs,
      op == "~" & lhs %in% lv_names ~ "bsp_" %p0% lhs %p0% "_mi" %p0% rhs,
      TRUE ~ paste(lhs, op, rhs)),
      mean = dplyr::coalesce(est, mean), median = dplyr::coalesce(est, median),
      sd = ifelse(is.na(est), 0, sd), mad = ifelse(is.na(est), 0, mad),
      q5 = dplyr::coalesce(est, q5), q95 = dplyr::coalesce(est, q95),
      rhat, ess_bulk, ess_tail)
}

#' @method bdlvm_get_stats lavaan
#' @export
bdlvm_get_stats.lavaan <- function(fit_obj, model_obj) {

  lv_names <- model_obj$lv_names
  i_names <- model_obj$i_names

  lavaan::partable(fit_obj)[,c("lhs", "op", "rhs", "free", "est", "se")] %>%
    dplyr::transmute(variable = dplyr::case_when(
      op == "=~" & lhs %in% lv_names ~ "bsp_" %p0% rhs %p0% "_mi" %p0% lhs,
      op == "~1" & lhs %in% c(lv_names, i_names) ~ "b_" %p0% lhs %p0% "_Intercept",
      op == "~~" & lhs %in% c(lv_names, i_names) & lhs == rhs ~ "sigma_" %p0% lhs,
      op == "~" & lhs %in% lv_names ~ "bsp_" %p0% lhs %p0% "_mi" %p0% rhs,
      TRUE ~ paste(lhs, op, rhs)), mean = est, sd = se)
}

#' @export bldvm_parametric_constraint
bldvm_parametric_constraint <- function(lv_names, M) {
  "constant(sqrt(datavar_" %p0%
    lv_names %p0% "LVi" %p0% 1:M %p0%
    " - (bsp_" %p0% lv_names %p0% "LVi" %p0% 1:M %p0%
    "[1]^2)*sigma_" %p0% lv_names %p0% "^2))"
}
