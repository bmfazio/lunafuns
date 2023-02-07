brmsh_constprior <- function(x){
  if(is.numeric(x))return(paste0("constant(", x, ")"))
  x
}

#' @export
brmsh_get_formulas <- function(x) {
  if(!"bform" %in% class(x))stop("Expected brms::bf() formula")
  if(is.null(x$forms))return(list(x$formula))
  x <- unlist(x$forms)
  x[which(sapply(x, rlang::is_formula))]
}

#' Get fully propagated predicitons from a `brmsfit` object
#'
#' @param brmsfit_obj An object of class `brmsfit`
#'
#' @return A list with one data.frame per posterior draw, each of which holds exogenous variables and propagated predictions of endogenous variables.
#' @export
#'
#' @examples
#' library(brms)
#' library(dplyr)
#'
#' {
#'   test_data <- data.frame(x = 1:5, y = 999, z1 = 0, z2 = 0)
#'   formulas <- bf(y ~ 0+x) + bf(z1 ~ 0+y) + bf(z2 ~ 0+y) + set_rescor(FALSE)
#'   get_prior(formulas, test_data) %>%
#'   arrange(resp) %>%
#'     mutate(prior = case_when(
#'       coef == "x" ~ "constant(1)", class == "sigma" ~ "normal(0, 0.1)",
#'       resp == "z1" ~ "normal(100, 10)", resp == "z2" ~ "normal(-100, 10)",
#'      TRUE ~ prior)) %>% slice(-c(1:2,5,8)) -> priors
#'    pars <- brm(formulas, test_data, prior = priors, sample_prior = "only",
#'               chains = 1, iter = 5150, warmup = 5000, thin = 50,
#'               silent = 0, refresh = 0)
#'   brmsh_pp_iter(pars)
#' }
brmsh_pp_iter <- function(brmsfit_obj) {
  forms_brms <- brmsh_get_formulas(brmsfit_obj$formula)
  forms_brms <- depmat_ordlist(formlst_depmat(forms_brms))
  forms_brms[[1]] <- forms_brms[[1]][forms_brms[[1]] != "..FILLERforSAMPLING"]

  lhs_vars <- unlist(forms_brms)
  init_frame <- brmsfit_obj$data
  init_frame <- cbind(init_frame[,setdiff(colnames(init_frame),lhs_vars),
                                 drop=FALSE],
                      init_frame[,lhs_vars])

  n_row <- nrow(init_frame)
  n_iter <- nrow(brms:::as.data.frame.brmsfit(brmsfit_obj))
  pred_frames <- list()

  for(i in 1:n_iter) {
    pred_frames[[i]] <- init_frame
    for(v in forms_brms) {
      array(brms::posterior_predict(
        brmsfit_obj,
        pred_frames[[i]],
        resp = v,
        draw_ids = i
      ), dim = c(1, n_row, length(v))
      )[1,,] -> pred_frames[[i]][,v]
    }
  }
  pred_frames
}
