#' @export dxp_rhat
dxp_rhat <- function(x, group_size, exclude_vars = c("lprior", "lp__")) {
  fit <- x$fit_diagnostics
  x <- dxh_prepare_plot_data(x)
  summ_fit <- x$fit_diagnostics
  summ_back <- x$backend_diagnostics

  rhat_bounds <- c(1, 1.5); scale_exponent <- 0.25
  fit$rhat <- dxh_rescale(fit$rhat, rhat_bounds)**scale_exponent
  summ_fit$rhat <- dxh_rescale(summ_fit$rhat, rhat_bounds)**scale_exponent

  values <- c(0, 0.02, 0.1, 1)**scale_exponent
  labels <- c("1", "1.01", "1.05", "1.5+")
  colors <- c("blue", "white", "red", "black")

  n_groups <- nrow(summ_back)/group_size
  if(n_groups%%1 != 0)stop("group_size must be multiple of summary table rows")

  list(
    heat =
      summ_fit %>%
      dplyr::filter(!variable %in% exclude_vars) %>%
      ggplot2::ggplot(ggplot2::aes(x = variable, y = sim_desc, fill = rhat)) +

      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(colors = colors, values = values, breaks = values,
        labels = labels, limits = c(0,1)) +

      ggplot2::geom_label(ggplot2::aes(y = sim_desc, x = 1, label = failed_chain_pct),
                          data = summ_back, color = "red", fill = "black") +

      ggplot2::geom_hline(yintercept = group_size*seq_len(n_groups-1) + 0.5,
                          linetype = "dashed", color = "gray80", size = 1) +
      ggplot2::scale_x_discrete(expand = c(0,0)) +
      ggplot2::scale_y_discrete(expand = c(0,0)) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth=2)) +
      ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::guides(linetype = "none")
    ,
    violin =
      fit %>%
      dplyr::filter(!variable %in% exclude_vars, !is.na(rhat)) %>%
      dplyr::mutate(y = "required by geom_violin") %>%
      ggplot2::ggplot(ggplot2::aes(x = rhat, y = y, color = rhat)) +
      ggplot2::geom_violin(alpha = 0.5, fill = "gray40", color = "gray50") +
      ggbeeswarm::geom_quasirandom(size = 1, method = "tukey", width = 0.25) +
      ggplot2::scale_x_continuous(position = "top", breaks = values[2:3], labels = NULL) +
      ggplot2::scale_color_gradientn(colors = colors, values = values, breaks = values,
                                     labels = labels, limits = c(0,1)) +

      ggplot2::geom_label(ggplot2::aes(x = 0.5, y = 1, label = failed_chain_pct),
                          data = summ_back %>% dplyr::filter(variable == sort(variable)[1]),
                          color = "red", fill = "black") +
      ggplot2::facet_grid(rows = dplyr::vars(forcats::fct_rev(sim_desc)),
                          cols = dplyr::vars(variable),
                          switch = "both") +

      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust = 1),
                     strip.text.y.left = ggplot2::element_text(angle = 0),
                     strip.text.x = ggplot2::element_text(angle = 90),
                     axis.title = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_line(color = "gray80"),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.spacing.x = unit(0.1, "lines"), panel.spacing.y = unit(0.1, "lines"))
  )
}
