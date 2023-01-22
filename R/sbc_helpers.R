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
  stopifnot({is.function(gen_fun);is.list(arg_list[[1]])})
  do.call(SBC::SBC_generator_function, c(gen_fun, arg_list[[1]], ...))
}

#' Adds more flexibility to the SBC generator
#'
#' @param generator SBC generator function
#' @param n_sims Number of datasets to simulate (usual behaviour)
#' @param n_reps Number of times to repeat a single simulated dataset (for other MCMC assessment purposes)
#'
#' @return Datasets in the format SBC expects
#' @export
sbch_generate <- function(generator, n_sims = NULL, n_reps = NULL) {
  # probably add way to execute both behaviors via a vector of values
  # that determine what to do for each generator that comes in
  # NO! actually, just detect the name of the argument coming in
  if(is.null(n_sims)&is.null(n_reps))stop("Must specify n_sims OR n_reps")
  if(!is.null(n_sims)&!is.null(n_reps))stop("Set only ONE of n_sims or n_reps")

  if(!is.null(n_sims))return(SBC::generate_datasets(generator, n_sims))

  single_dataset <- SBC::generate_datasets(generator, 1)

  draws_attr <- attributes(single_dataset$variables)
  single_dataset$variables <- single_dataset$variables[rep(1, n_reps),]
  attr(single_dataset$variables, "nchains") <- draws_attr$nchains
  attr(single_dataset$variables, "class") <- draws_attr$class
  attr(single_dataset$variables, "dimnames")$draw <- as.character(1:n_reps)

  single_dataset$generated <- single_dataset$generated[rep(1, n_reps)]

  single_dataset
}

#' `compute_SBC` wrapper
#'
#' @param datasets an object of class `SBC_datasets`
#' @param backend model sampler made by some `SBC_backend_*`function
#'
#' @export
sbch_run <- function(datasets, backend) {
  SBC::compute_SBC(datasets, backend, keep_fits = TRUE, cache_mode = "none")
}
