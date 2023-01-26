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
#' test_data <- data.frame(x = 1:5, y = 999, z1 = 0, z2 = 0)
#' formulas <- bf(y ~ 0+x) + bf(z1 ~ 0+y) + bf(z2 ~ 0+y) + set_rescor(FALSE)
#' get_prior(formulas, test_data) %>%
#' arrange(resp) %>%
#'   mutate(prior = case_when(
#'     coef == "x" ~ "constant(1)", class == "sigma" ~ "normal(0, 0.1)",
#'     resp == "z1" ~ "normal(100, 10)", resp == "z2" ~ "normal(-100, 10)",
#'     TRUE ~ prior)) %>% slice(-c(1:2,5,8)) -> priors
#'   pars <- brm(formulas, test_data, prior = priors, sample_prior = "only",
#'             chains = 1, iter = 10500, warmup = 10000, thin = 50,
#'             silent = 0, refresh = 0)
#'  brmsh_pp_iter(pars)
brmsh_pp_iter <- function(brmsfit_obj) {
  forms_brms <- unlist(brmsfit_obj$formula$forms)
  forms_brms <- forms_brms[which(sapply(forms_brms, rlang::is_formula))]
  forms_brms <- depmat_ordlist(formlst_depmat(forms_brms))

  lhs_vars <- unlist(forms_brms)
  init_frame <- brmsfit_obj$data
  init_frame <- cbind(init_frame[,setdiff(colnames(init_frame),lhs_vars),
                                 drop=FALSE],
                      init_frame[,lhs_vars])

  n_iter <- nrow(brms:::as.data.frame.brmsfit(brmsfit_obj))
  pred_frames <- list()

  for(i in 1:n_iter) {
    pred_frames[[i]] <- init_frame
    for(v in forms_brms) {
      brms::posterior_predict(
        brmsfit_obj,
        pred_frames[[i]],
        draw_ids = i
      )[1,,][,v,drop=FALSE] -> pred_frames[[i]][,v]
    }
  }
  pred_frames
}
