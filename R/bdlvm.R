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

pregen_check_args <- function(f, N, M, add_i, mu, sigma,
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

#' @export
bdlvm_generate_sem <- function(pregen, ...) {
  sbch_generate(do.call(SBC_generator_brms, pregen$genargs),
                makena = pregen$lv_names, ...)
}
