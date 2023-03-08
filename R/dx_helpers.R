#' @export dxh_statsp
dxh_statsp <- function(x, dx = c("rhat", "ess_bulk", "ess_tail")) {
  names(dx) <- dx
  x <- dxh_prep_stats(x)

  list( # make function specification flexible to add/remove plots at will
    heat = lapply(dx, (\(y)dxh_heat(x, y))),
    boxp = lapply(dx, (\(y)dxh_boxp(x, y)))
  )
}

#' @export dxh_prep_stats
dxh_prep_stats <- function(x) {
  if(is.null(x$sim_desc))x$sim_desc <- ""
  x %>%
    dplyr::transmute(
      sim_desc, variable,
      rhat, ess_bulk, ess_tail) %>%
    tidyr::pivot_longer(!c(sim_desc, variable)) %>%
    dplyr::group_by(sim_desc, variable, name) %>%
    dplyr::filter(!variable %in% c("lp__", "lprior"))
}

#' @export dxh_heat
dxh_heat <- function(x, dx) {
  if(!dx%in%c("ess_bulk","ess_tail", "rhat"))stop("Wrong dx")

  colors <- c("red", "white", "blue")
  values <- c(0, 0.125, 1)
  breaks <- c(0, 250, 1000, 2000)
  tfun <- function(x)psych::winsor.mean(x, trim = 0.1)
  trans <- "identity"

  if(dx == "rhat") {
    colors <- c(rev(colors), "black")
    values <- c(0, sqrt(0.02), sqrt(0.1), 1)
    breaks <- c(1, 1.01, 1.05, 1.5)
    tfun <- function(x)exp(mean(log(x)))
    trans <- "sqrt"
  }

  x %>%
    dplyr::summarise(value = tfun(value), .groups = "drop") %>%
    dplyr::filter(name == dx) %>%
    dplyr::mutate(value = dxh_rescale(value, breaks)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sim_desc, y = variable, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = colors, values = values,
      breaks = (breaks-min(breaks))/(max(breaks)-min(breaks)),
      labels = breaks %p0% c(rep("", length(breaks)-1), "+"),
      trans = trans, limits = c(0,1)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.border = ggplot2::element_rect(colour = "black",
                                                        fill=NA, linewidth=2)) +
    ggplot2::ggtitle(dx) +
    ggplot2::xlab("")
}

#' @export dxh_boxp
dxh_boxp <- function(x, dx) {
  if(!dx%in%c("ess_bulk","ess_tail", "rhat"))stop("Wrong dx")

  x %>%
    dplyr::filter(name == dx) %>%
    dplyr::mutate(valuex = value - 1.05) %>%
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(y = 0),
                         height = 0.2, alpha = 0.5, color = "gray42") +
    ggplot2::scale_x_continuous(
      position = "top",
      trans = if(dx=="rhat"){dxh_scale_rhat}else{"identity"}) +
    ggplot2::ggtitle(dx) +
    ggplot2::facet_grid(cols = dplyr::vars(sim_desc),
                        rows = dplyr::vars(forcats::fct_rev(variable)),
                        switch = "both") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust = 1),
                   strip.text.y.left = ggplot2::element_text(angle = 0),
                   strip.text.x = ggplot2::element_text(angle = 45),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

#' @export dxh_viop
dxh_viop <- function(x, dx) {
  if(!dx%in%c("ess_bulk","ess_tail", "rhat"))stop("Wrong dx")

  x %>%
    dplyr::filter(name == dx) %>%
    dplyr::mutate(valuey = "whatever") %>%
    ggplot2::ggplot(ggplot2::aes(x = value, y = valuey)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(height = 0.2, alpha = 0.5, color = "gray42") +
    ggplot2::scale_x_continuous(
      position = "top",
      trans = if(dx=="rhat"){dxh_scale_rhat}else{"identity"}) +
    ggplot2::ggtitle(dx) +
    ggplot2::facet_grid(cols = dplyr::vars(sim_desc),
                        rows = dplyr::vars(forcats::fct_rev(variable)),
                        switch = "both") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust = 1),
                   strip.text.y.left = ggplot2::element_text(angle = 0),
                   strip.text.x = ggplot2::element_text(angle = 45),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
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

#' @export dxh_heat2
dxh_heat2 <- function(x, dx) {
  if(!dx%in%c("ess_bulk","ess_tail", "rhat"))stop("Wrong dx")

  colors <- c("red", "white", "blue")
  values <- c(0, 0.125, 1)
  breaks <- c(0, 250, 1000, 2000)
  tfun <- function(x)psych::winsor.mean(x, trim = 0.1)
  trans <- "identity"

  if(dx == "rhat") {
    colors <- c(rev(colors), "black")
    values <- c(0, sqrt(0.02), sqrt(0.1), 1)
    breaks <- c(1, 1.01, 1.05, 1.5)
    tfun <- function(x)exp(mean(log(x)))
    trans <- "sqrt"
  }

  x %>%
    dplyr::summarise(value = tfun(value), .groups = "drop") %>%
    dplyr::filter(name == dx) %>%
    dplyr::mutate(value = dxh_rescale(value, breaks)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sim_desc, y = variable, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(data = x %>%
                          dplyr::filter(variable == "b_y1_Intercept",
                                        sim_mis != 0),
                        mapping = ggplot2::aes(label = sim_mis), color = "red") +
    ggplot2::scale_fill_gradientn(
      colors = colors, values = values,
      breaks = (breaks-min(breaks))/(max(breaks)-min(breaks)),
      labels = breaks %p0% c(rep("", length(breaks)-1), "+"),
      trans = trans, limits = c(0,1)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.border = ggplot2::element_rect(colour = "black",
                                                        fill=NA, linewidth=2)) +
    ggplot2::ggtitle(dx) +
    ggplot2::xlab("")
}
