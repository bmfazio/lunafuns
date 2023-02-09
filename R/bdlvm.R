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
#' 2. Feed the "base formula" into `bdlvm_add_items` to produce a model object that will have a "full formula" and template dataset that has added variables corresponding to the items for each latent variable.
#' 3. Feed the model object to `bdlvm_make_data` together with the number of desired datasets to simulate the data (this is a SLOW STEP as it uses `brms::posterior_predict` iteratively to populate the entire graph).
#' 4. Feed the dataset and model objects to `bdlvm_compute_SBC` in order to fit the model to each dataset via a modified version of `SBC::compute_SBC`. The same generative priors will be used for fitting by default, but see the arguments for ways to specify different priors. EXPERIMENTAL: the backend object SBC requires will be generated at this point in order to guarantee compilation of the latest model specification. YES I LEARNED THIS IN A VERY PAINFUL MANNER.
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
#' bdlvm_add_items(bf(y1|mi()~1), data = 50, M = 3)
bdlvm_add_items <- function(formula, data, M, ...) {
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

  out$lv_names <- lv_names
  out$i_names <- i_names
  out
}

#' @describeIn bdlvm Make GENERATIVE PRIORS out of a base model object (formula+data)
#'
#' @param ... Look at defaults for available parameters and defaults.
#' @param explain Use this to print an explanation of the way order in which values with length > 1 are assigned to specific parameters.
#'
#' @return `bdlvm_add_prior`: Object with edited formula and data + priors unless prior_only = TRUE
#' @export
#'
#' @examples
#' suppressMessages(library(brms))
#' bdlvm_add_items(bf(y1|mi()~1), data = 50, M = 3) -> sem_formula
#' bdlvm_add_prior(sem_formula, prior_only = TRUE)
bdlvm_add_prior <- function(model_obj,
                            mu = 0, sigma = 1, betas = 1,
                            lf = 1, noise = 1, item_mu = 0,
                            explain = FALSE, prior_only = FALSE, ...) {

  model_obj$formula <- model_obj$formula + brms::bf(..FILLERforSAMPLING ~ 0) + set_rescor(FALSE)
  model_obj$data <- cbind(model_obj$data,
                          ..FILLERforSAMPLING = rnorm(nrow(model_obj$data)))
  model_obj$prior <- brms::get_prior(model_obj$formula, model_obj$data)
  lv_names <- model_obj$lv_names
  n_lv <- length(lv_names)
  i_names <- model_obj$i_names
  M <- length(i_names)/n_lv

  # Translate hyperparameter args into full vectors
  mu <- brmsh_constprior(mu)
  sigma <- brmsh_constprior(sigma)
  lf <- brmsh_constprior(lf)
  noise <- brmsh_constprior(noise)
  item_mu <- brmsh_constprior(item_mu)
  betas <- brmsh_constprior(betas)

  mu <- rep(mu, n_lv/length(mu))
  sigma <- rep(sigma, n_lv/length(sigma))
  lf <- rep(lf, times = n_lv*M/length(lf))
  noise <- rep(noise, times = n_lv*M/length(noise))
  item_mu <- rep(item_mu, times = n_lv*M/length(item_mu))
  # Beta length is not checked

  # Prepare indexes for value assignment
  mu_idx <- intersect(which_order(lv_names, model_obj$prior$resp),
                      which(model_obj$prior$class == "Intercept"))
  sigma_idx <- intersect(which_order(lv_names, model_obj$prior$resp),
                         which(model_obj$prior$class == "sigma"))
  lf_idx <- intersect(which_order(i_names, model_obj$prior$resp),
                      which(model_obj$prior$class == "b" &
                              model_obj$prior$coef != ""))
  noise_idx <- intersect(which_order(i_names, model_obj$prior$resp),
                         which(model_obj$prior$class == "sigma"))
  itmu_idx <- intersect(which_order(i_names, model_obj$prior$resp),
                        which(model_obj$prior$class == "Intercept"))
  betas_idx <- {
    tmp <- cbind(1:nrow(model_obj$prior), model_obj$prior)
    tmp <- tmp[which_order(lv_names, tmp$resp),]
    tmp[which_order(paste0("mi", lv_names), tmp$coef),1]
  }

  model_obj$prior[mu_idx, "prior"] <- mu
  model_obj$prior[sigma_idx, "prior"] <- sigma
  model_obj$prior[lf_idx, "prior"] <- lf
  model_obj$prior[noise_idx, "prior"] <- noise
  model_obj$prior[itmu_idx, "prior"] <- item_mu
  model_obj$prior[betas_idx, "prior"] <- betas
  model_obj$prior[model_obj$prior$resp == "FILLERforSAMPLING",
                  "prior"] <- "normal(0,1.67)"
  model_obj$prior$source <- "bdlvm_add_prior"

  # Add explainer function if requested via explain = TRUE
  if(explain) {
    cat("\nPriors are being assigned in the following order:\n",
        "< LV intercepts & variances>\n",
        paste(model_obj$prior[mu_idx, "resp"], collapse = ", "), "\n",
        "< Item loading factors, intercepts & variances>\n",
        paste(model_obj$prior[lf_idx, "resp"], collapse = ", "), "\n",
        "< LV-to-LV regression coefficients >\n",
        paste(substring(model_obj$prior[betas_idx, "coef"], 3),
              "->", model_obj$prior[betas_idx, "resp"], collapse = "\n"),
        sep = "")
  }

  if(prior_only)return(model_obj$prior)
  model_obj
}

bdlvm_generate_datasets <- function(model_obj, ...) {
  generator_args <- list(formula = model_obj$formula,
                         data = model_obj$data,
                         prior = model_obj$prior)

  # PENDING: set some brms values to minimise parameter correlations
  datasets <- sbch_generate(do.call(SBC::SBC_generator_brms, generator_args),
                            makena = model_obj$lv_names, ...)

  list(datasets = datasets, prior = model_obj$prior)
}

#' @export bdlvm_make_data
bdlvm_make_data <- function(model_obj, ...) {
  bdlvm_generate_datasets(bdlvm_add_prior(model_obj, ...), ...)
}

#' Title
#'
#' Bottom text
#'
#' @param model_obj An object made with `bdlvm_add_items`
#' @param data_obj An object made with `bdlvm_generate_datasets`
#' @param ... Arguments passed to `brms::brm`
#' @param fit_prior_args List of named arguments passed to `make_fit_prior`
#' @param cache_mode Disabling SBC cache by default
#' @param cache_location Disabling SBC cache by default
#'
#' @return Just a regular `SBC_results`-class object.
#' @export bdlvm_compute_SBC
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
#' model_object <- bdlvm_add_items(base_formula, data = 5, M = 5)
#'
#' # data generation setup
#' generator_arguments <- bdlvm_add_prior(model_object)
#' dataset_object <- bdlvm_generate_datasets(generator_arguments, n_sims = 2)
#'
#' # fit (with a potentially different prior)
#' fits_object <- bdlvm_compute_SBC(model_object, dataset_object,
#'                                  fit_prior_args = list(
#'                                    varmethod = "varmi",
#'                                    global_prior = "normal(0,3)"
#'                                  ))
bdlvm_compute_SBC <- function(model_obj, data_obj, ...,
                              fit_prior_args = list(varmethod = "varmi"),
                              cache_mode = "none", cache_location = NULL) {

  backend_obj <- make_backend(model_obj, data_obj$prior, ...,
                              fit_prior_args = fit_prior_args)

  compute_SBC_stanvars(
    datasets = data_obj$datasets,
    backend = backend_obj,
    cache_mode = cache_mode,
    cache_location = cache_location)
}
make_backend <- function(model_obj, prior_obj, ...,
                         backend = "cmdstanr", chains = 4, stanvar_func = NULL) {

  fit_prior <- do.call(make_fit_prior,
                       c(list(model_obj, prior_obj), list(...)$fit_prior_args))
  if(is.null(stanvar_func))stanvar_func <- stanvar_func_datavar(model_obj)

  arg_list <- list(formula = model_obj$formula + brms::set_rescor(FALSE),
                   template_data = model_obj$data, prior = fit_prior,
                   backend = backend, chains = chains,
                   stanvars = stanvar_func, ...)

  do.call(sbch_brms_backend, arg_list)
}
make_fit_prior <- function(model_obj, prior_obj, varmethod, ...,
                               lf1_fixto1 = TRUE, global_prior = NULL,
                               insert_prior = NULL) {

  if(varmethod == "copy")return(prior_obj)

  lv_names <- model_obj$lv_names
  n_lv <- length(lv_names)
  i_names <- model_obj$i_names
  M <- length(i_names)/n_lv
  mu_idx <- intersect(which_order(lv_names, prior_obj$resp),
                      which(prior_obj$class == "Intercept"))
  itmu_idx <- intersect(which_order(i_names, prior_obj$resp),
                        which(prior_obj$class == "Intercept"))

  out_prior <- prior_obj
  out_prior <- out_prior[out_prior$resp != "FILLERforSAMPLING" &
                           out_prior$prior != "",]

  lf1_idx <- out_prior$resp %in%
    i_names[1+0:(n_lv-1)*M] &
    out_prior$coef %in% paste0("mi", lv_names)

  if(lf1_fixto1)out_prior[lf1_idx, "prior"] <- "constant(1)"
  out_prior[lf1_idx, "source"] <- "item1_lf"

  isigma_idx <- intersect(which_order(i_names, out_prior$resp),
                          which(out_prior$class %in% "sigma"))

  if(varmethod == "varmi") {
    out_prior[isigma_idx, "prior"] <- "constant(sqrt(datavar_" %p0%
      i_names %p0% " - (bsp_" %p0% out_prior[isigma_idx, "resp"] %p0%
      "[1]^2)*variance(Ymi_" %p0% rep(lv_names, each=M) %p0% ")))"
    out_prior[isigma_idx, "source"] <- "varvarmethod_constrain_variance"
  } else if(varmethod == "sigma") {
    if(length(lv_names)>1)warning("Item variance constraint is not properly implemented for dependent latent variables!")
    out_prior[isigma_idx, "prior"] <- "constant(sqrt(datavar_" %p0%
      i_names %p0% " - (bsp_" %p0% out_prior[isigma_idx, "resp"] %p0%
      "[1]^2)*sigma_" %p0% rep(lv_names, each=M) %p0% "^2))"
    out_prior[isigma_idx, "source"] <- "varvarmethod_constrain_sigma"
  } else if(varmethod == "default") {
    warning("Assuming t scale as 2.5 but brms adjusts based on data!")
    out_prior[isigma_idx, "prior"] <- "student_t(3, 0, 2.5)"
    out_prior[isigma_idx, "source"] <- "varvarmethod_imitate_default"
  } else {
    out_prior[isigma_idx, "prior"] <- varmethod
    out_prior[isigma_idx, "source"] <- "varvarmethod_custom_string"
  }

  if(!is.null(global_prior)) {
    global_idx <- out_prior$class != "Intercept" & out_prior$source == "bdlvm_add_prior"
    out_prior[global_idx, "prior"] <- global_prior
    out_prior[global_idx, "source"] <- "global_prior"
  }

  if(!is.null(insert_prior)) {
    stopifnot(is.data.frame(insert_prior))
    match_cols <- out_prior[,setdiff(colnames(insert_prior), "prior")]
    not_priors <- colnames(insert_prior) != "prior"
    for(i in 1:nrow(insert_prior)) {
      match_rows <- apply(
        match_cols, 1, identical, unlist(insert_prior[i, not_priors]))

      out_prior[match_rows, "prior"] <- insert_prior[i, "prior"]
      out_prior[match_rows, "source"] <- "insert_prior_row_" %p0% i
    }
  }

  out_prior <- rbind(out_prior[lf1_idx,],
                     out_prior[isigma_idx,],
                     out_prior[-union(which(lf1_idx), isigma_idx),])

  out_prior[out_prior$class != "Intercept" | out_prior$source != "bdlvm_add_prior", ]
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
#' Prepare arguments required for a SEM `SBC` generator via `brms`
#'
#' @param formula SEM formula specifying exogenous and latent variables. Any variable using the `mi()` syntax is assumed to be a latent variable and will be assigned receive priors and item variables according to the arguments that follow.
#' @param data A template dataframe that will be passed to `brm`. The number of rows determines the number of observations that will be generated. Only values for exogenous variables are used, all latent variables will be overwritten by the generated values. If no exogenous variables are present, a single integer can be passed to indicate the number of observations to generate.
#' @param M The amount of items per latent variable. It is currently assumed to be identical for all.
#' @param mu,sigma Intercepts and variances for the latent variables. Can be fixed values or strings describing valid prior distributions to be used by `brms` for parameter draws. Must be of length 1 or equal to the number of latent variables in the formula.
#' @param betas Coefficients for regressions between latent variables. Currently only supporting linear  predictors on the location parameters. Fixed values or strings, as above. Length not checked but should be equal to amount of relationships in the model, see `explain` argument for more.
#' @param lf,noise,item_mu Loading factors, variances and intercepts for the items. Fixed values or strings, as above. Must be of length 1, M or equal to the total number of items in the model.
#' @param add_items `TRUE` by default. Will cause item terms and priors to be added to the output model formula, dataframe and priors.
#' @param mode Default value of `"predict"`. Generates output with the expectation that `brms` will be used to draw parameter values and then generate data via `posterior_predict`. Can be changed to `"fit"`, which will instead fit the full model to generate values via `mi()` imputation.
#' @param explain Will print text explaining the order in which priors will be assigned to model parameters.
#' @param ... Currently unused.
#'
#' @return A list with `formula`, `data` and `prior` components, ready to be passed into `SBC_generator_brms`. Be aware that priors for exogenous variables are left unspecified, so the user must finish setting them up in order to produce sensible results.
#' @export
#'
#' @examples
#' library(brms)
#' bdlvm_brms_pregen(
#'  bf(y1|mi()~1) +
#'    bf(y21|mi()~mi(y1)) + bf(y22|mi()~mi(y1)) +
#'    bf(y3|mi()~mi(y21)+mi(y22)),
#'  M = 3, data = mtcars[1:5,1:2],
#'  sigma = c(10, 2, 2, 2), noise = 0.01,
#'  explain = TRUE) -> .bdlvm_pregen_example
#'
#'  .bdlvm_pregen_example$genargs$formula
#'  head(.bdlvm_pregen_example$prior, 10)
bdlvm_brms_pregen <- function(formula, data = 1, M = 1,
                              mu = 0, sigma = 1, betas = 1,
                              lf = 1, noise = 1, item_mu = 0,
                              add_items = TRUE, mode = "predict",
                              explain = FALSE, ...) {

  # Identify LVs (formula is not checked for appropriate use of mi() on RHS)
  f_list <- brmsh_get_formulas(formula)
  lv_names <- get_lhs_vars(f_list)[get_lhs_mi(f_list)]

  # Data checks
  if(is.data.frame(data)){N<-nrow(data)}else{N<-data}
  n_lv <- pregen_check_args(lv_names, N, M, add_items, mu, sigma,
                            lf, noise, betas, item_mu, mode)

  # Prepare output
  out <- list()
  # 1. Construct template dataset
  out$data <- as.data.frame(matrix(0, nrow = N, ncol = n_lv*(1+M)))
  colnames(out$data) <- rep(lv_names,
                            each = M+1) %p0% c("", rep("LVi", M) %p0% 1:M)
  i_names <- colnames(out$data)[-(1+(0:(n_lv-1))*(M+1))]

  if(is.data.frame(data)) {
    out$data <- cbind(data[,setdiff(colnames(data), lv_names)], out$data)
  }
  if(mode == "predict"){out$data <- cbind(out$data, ..FILLERforSAMPLING = rnorm(N))}

  # 2. Extend formula with LV items (for now assume always add_items = TRUE)
  out$formula <- formula
  for(lv in lv_names) {
    for(m in 1:M) {
      out$formula <- out$formula +
        brms::bf(paste0(lv, "LVi", m, "| mi() ~ mi(", lv, ")"))
    }
  }
  fit_formula <- out$formula + set_rescor(FALSE)
  if(mode == "predict"){out$formula <- out$formula + brms::bf(..FILLERforSAMPLING ~ 0)}
  out$formula <- out$formula + set_rescor(FALSE)

  # 3. Construct priors
  out$prior <- brms::get_prior(out$formula, out$data)
  ## Translate hyperparameter args into full vectors
  if(is.numeric(mu))mu <- paste0("constant(", mu, ")")
  if(is.numeric(sigma))sigma <- paste0("constant(", sigma, ")")
  if(is.numeric(lf))lf <- paste0("constant(", lf, ")")
  if(is.numeric(noise))noise <- paste0("constant(", noise, ")")
  if(is.numeric(item_mu))item_mu <- paste0("constant(", item_mu, ")")
  if(is.numeric(betas))betas <- paste0("constant(", betas, ")")

  mu <- rep(mu, n_lv/length(mu))
  sigma <- rep(sigma, n_lv/length(sigma))
  lf <- rep(lf, times = n_lv*M/length(lf))
  noise <- rep(noise, times = n_lv*M/length(noise))
  item_mu <- rep(item_mu, times = n_lv*M/length(item_mu))
  ## Checking for correct beta length is non-trivial but we hope
  ## errors are caught when non-multiple lengths are entered.

  ## Prepare indexes for value assignment
  mu_idx <- intersect(which_order(lv_names, out$prior$resp),
                      which(out$prior$class == "Intercept"))
  sigma_idx <- intersect(which_order(lv_names, out$prior$resp),
                         which(out$prior$class == "sigma"))
  lf_idx <- intersect(which_order(i_names, out$prior$resp),
                      which(out$prior$class == "b" & out$prior$coef != ""))
  noise_idx <- intersect(which_order(i_names, out$prior$resp),
                         which(out$prior$class == "sigma"))
  itmu_idx <- intersect(which_order(i_names, out$prior$resp),
                        which(out$prior$class == "Intercept"))
  betas_idx <- {
    tmp <- cbind(1:nrow(out$prior), out$prior)
    tmp <- tmp[which_order(lv_names, tmp$resp),]
    tmp[which_order(paste0("mi", lv_names), tmp$coef),1]
  }

  out$prior[mu_idx, "prior"] <- mu
  out$prior[sigma_idx, "prior"] <- sigma
  out$prior[lf_idx, "prior"] <- lf
  out$prior[noise_idx, "prior"] <- noise
  out$prior[itmu_idx, "prior"] <- item_mu
  out$prior[betas_idx, "prior"] <- betas
  if(mode == "predict"){out$prior[out$prior$resp == "FILLERforSAMPLING",
                                  "prior"] <- "normal(0,1.67)"}
  out$prior$source <- "pregen"

  # Add explainer function if requested via explain = TRUE
  if(explain) {
    cat("\nPriors are being assigned in the following order:\n",
        "< LV intercepts & variances>\n",
        paste(out$prior[mu_idx, "resp"], collapse = ", "), "\n",
        "< Item loading factors, intercepts & variances>\n",
        paste(out$prior[lf_idx, "resp"], collapse = ", "), "\n",
        "< LV-to-LV regression coefficients >\n",
        paste(substring(out$prior[betas_idx, "coef"], 3),
              "->", out$prior[betas_idx, "resp"], collapse = "\n"),
        sep = "")
  }

  out$refresh <- 0

  out <- list(genargs = out,
              fit_formula = fit_formula,
              lv_names = lv_names,
              prior_df = (\(){
                out_df <- out$prior[-c(mu_idx, itmu_idx),]
                out_df <- out_df[out_df$resp != "FILLERforSAMPLING" &
                                   out_df$prior != "",]

                lf1_idx <- out_df$resp %in% i_names[1+0:(n_lv-1)*M] &
                  out_df$coef %in% paste0("mi", lv_names)
                out_df[lf1_idx, "source"] <- "item1_lf"

                isigma_idx <- intersect(which_order(i_names, out_df$resp),
                                        which(out_df$class %in% "sigma"))
                out_df[isigma_idx, "prior"] <- "constant(sqrt(datavar_" %p0%
                  i_names %p0% " - (bsp_" %p0%
                  out_df[isigma_idx, "resp"] %p0% "[1]^2)*variance(Ymi_" %p0%
                  rep(lv_names, each=M) %p0% ")))"
                out_df[isigma_idx, "source"] <- "var_constraint"
                # prior(constant(var_item1-bsp_item1[1]*sigma_lv^2), class = "sigma", resp = "item1")
                # sigma_<LVname>: pattern for each LVs variance
                # bsp_<itemname>[1]: pattern for each LV -> item coefficient name
                # sigma <itemname>: pattern for each item's variance
                # datavar_<itemname>: pattern for empirically calculated variance for stanvars

                out_df <- rbind(out_df[lf1_idx,],
                                out_df[isigma_idx,],
                                out_df[-union(which(lf1_idx), isigma_idx),])
                out_df
              })(),
              make_datavars = (function(data) {
                var_data <- brms::stanvar(block = "functions", scode = "// Empty stanvars")
                for(i in i_names) {
                  var_data <- var_data +
                    brms::stanvar(x = var(data[,i]),
                                  name = paste0("datavar_", i),
                                  block = "data")
                }
                var_data
              }
              ))
  out
}

bdlvm_check_args <- function(f, N, M, add_i, mu, sigma,
                             lf, noise, betas, item_mu, mode) {
  check_length(mode, 1);
  if(!mode%in%c("predict",
                "estimate"))stop("mode must be 'predict' or 'estimate'")
  L <- length(f)
  check_length(N, 1); check_length(M, 1)
  check_length(mu,1,L);check_length(sigma,1,L)
  check_length(lf,1,M,L*M);check_length(noise,1,M,L*M)
  check_length(item_mu,1,M,L*M)
  L
}
