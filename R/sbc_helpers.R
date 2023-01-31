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

#' Construct `brms` backend for SBC
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
  if(is.function(dots$formula))dots$formula <- dots$formula(template_data)
  if(is.function(dots$prior))dots$prior <- dots$prior(template_data)

  do.call(SBC::SBC_backend_brms,
          c(dots, list(template_data = template_data)))
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
#'
#' @return Datasets in the format SBC expects
#' @export
sbch_generate <- function(generator, ..., n_sims = NULL, n_reps = NULL,
                          mvbrms = TRUE) {
  if(is.null(n_sims)&is.null(n_reps))stop("Must specify n_sims OR n_reps")
  if(!is.null(n_sims)&!is.null(n_reps))stop("Set only ONE of n_sims or n_reps")

  generate_fun <- sbch_generate_semdata
  if(!mvbrms)generate_fun <- SBC::generate_datasets

  if(!is.null(n_sims))return(generate_fun(generator, n_sims))

  single_dataset <- generate_fun(generator, 1)

  draws_attr <- attributes(single_dataset$variables)
  single_dataset$variables <- single_dataset$variables[rep(1, n_reps),]
  attr(single_dataset$variables, "nchains") <- draws_attr$nchains
  attr(single_dataset$variables, "class") <- draws_attr$class
  attr(single_dataset$variables, "dimnames")$draw <- as.character(1:n_reps)

  single_dataset$generated <- single_dataset$generated[rep(1, n_reps)]

  single_dataset
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
sbch_run <- function(arg_row, n_reps = NULL, n_sims = NULL,
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
#' @param generator A generator built with `bdlvm_brms_pregen` and `SBC_generator_brms`
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
    discard_filler <- which(grepl("(..FILLERforSAMPLING)|(^.+LVi[0-9]+$)",
                                  names(generated[[1]])))
    mi_names <- colnames(generated[[1]][,-discard_filler])
    gen_n <- nrow(generated[[1]])
    draws <- posterior::as_draws_matrix(prior_fit_brms$fit)
    for(i in 1:n_sims) {
      if (generator$generate_lp) {
        ll <- log_lik(prior_fit_brms, newdata = generated[[i]],
                      draw_ids = i, cores = 1)
        log_likelihoods[i] <- sum(ll)
      }
      sim_values <-unlist(generated[[i]][,-discard_filler])
      names(sim_values) <- paste0(
        rep(paste0("Ymi_", mi_names), each = gen_n),
        "[", 1:gen_n, "]")
      draws[i, which_order(names(sim_values), colnames(draws))] <- sim_values
    }
  })
  gd_tmp
})()
