#' @export
#' @noRd
print.pc_single = \(x,...) {
  print(x$summary)
}

#' @export
#' @noRd
print.pc_boot = \(x,...) {
  print(x$causality)
}

#' @export
#' @noRd
print.pc_ops = \(x,...) {
  res = x$param
  cat(paste0("The suggested E, k, tau is ", res[1], ", ", res[2], " and ", res[3]), "\n")
}

#' @export
#' @noRd
plot.pc_boot = \(x, family = "serif",
                 xbreaks = NULL, xlimits = NULL,
                 ybreaks = seq(0, 1, by = 0.1), 
                 ylimits = NULL, ...){
  causdf = x$causality
  if (ncol(causdf) == 3) names(causdf) = c("libsizes", "type", "q50")
  if(is.null(xbreaks)) xbreaks = causdf$libsizes
  if(is.null(xlimits)) xlimits = c(min(xbreaks) - 1,max(xbreaks) + 1)
  if(is.null(ylimits)) ylimits = range(causdf$q50) + c(-0.05, 0.05)

  ggplot2::ggplot(causdf, ggplot2::aes(x = libsizes, y = q50, color = type)) +
    ggplot2::geom_line(linewidth = 1.25) +
    ggplot2::scale_color_manual(name = NULL,
                                values = c("dark" = "#6A0DAD", 
                                           "negative" = "#FF4136", 
                                           "positive" = "#0074D9")) +
    ggplot2::scale_x_continuous(breaks = xbreaks, limits = xlimits,
                                expand = c(0, 0), name = "Library Size") +
    ggplot2::scale_y_continuous(breaks = ybreaks, limits = ylimits,
                                expand = c(0, 0), name = "Causal Strength") +
    ggplot2::theme_bw(base_family = family) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
