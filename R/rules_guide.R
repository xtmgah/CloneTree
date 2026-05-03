#' Plot a simple guide to CloneTree QC rules
#'
#' Draws a small explanatory figure for the two pairwise constraints used by
#' [clone_qc()]: the pigeonhole principle and crossing rule. The red outcome
#' highlights the case CloneTree treats as incompatible: two clusters that must
#' branch by crossing but cannot both fit as separate branches in at least one
#' tumor region by pigeonhole.
#'
#' @param filename Optional output path. When supplied, the plot is saved and
#'   returned invisibly.
#' @param width,height Optional output dimensions for `filename`.
#'
#' @return A `ggplot` object, or invisibly the filename if saved.
#' @export
plot_qc_rules_guide <- function(filename = NULL, width = NULL, height = NULL) {
  font_family <- ct_font_family()
  clone_a <- "#A72B30"
  clone_b <- "#68AFC3"
  red <- "#C33A3A"
  grey <- "#F2F2F2"

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0.2, xmax = 4.8, ymin = 3.8, ymax = 7.0, fill = grey, color = NA) +
    ggplot2::annotate("rect", xmin = 5.2, xmax = 9.8, ymin = 3.8, ymax = 7.0, fill = grey, color = NA) +
    ggplot2::annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.2, ymax = 3.3, fill = "#FBFBFB", color = NA) +
    ggplot2::annotate("text", x = 2.5, y = 6.65, label = "Pigeonhole principle", family = font_family, fontface = "bold", size = 4.2) +
    ggplot2::annotate("text", x = 7.5, y = 6.65, label = "Crossing rule", family = font_family, fontface = "bold", size = 4.2) +
    ggplot2::annotate("text", x = 5.0, y = 3.0, label = "Incompatible pair", family = font_family, fontface = "bold", size = 4.2, color = red) +
    ggplot2::annotate("text", x = 5.0, y = 2.55, label = "If both rules are triggered, these two clusters cannot co-occur in one solution.", family = font_family, size = 3.5) +
    ggplot2::annotate("text", x = 5.0, y = 2.05, label = "CloneTree excludes at least one and keeps all best compatible solutions.", family = font_family, size = 3.5) +
    ggplot2::annotate("text", x = 5.0, y = 1.15, label = "PH + X", family = font_family, fontface = "bold", size = 5.0, color = red) +
    ggplot2::annotate("text", x = 5.8, y = 1.15, label = "=", family = font_family, size = 5.0) +
    ggplot2::annotate("text", x = 6.7, y = 1.15, label = "exclude one", family = font_family, fontface = "bold", size = 4.6, color = red) +
    ggplot2::annotate("rect", xmin = 1.1, xmax = 2.0, ymin = 4.45, ymax = 6.05, fill = clone_a, alpha = 0.85) +
    ggplot2::annotate("rect", xmin = 2.25, xmax = 3.15, ymin = 4.45, ymax = 5.80, fill = clone_b, alpha = 0.85) +
    ggplot2::annotate("segment", x = 0.95, xend = 3.3, y = 4.45, yend = 4.45, linewidth = 0.35) +
    ggplot2::annotate("text", x = 1.55, y = 4.25, label = "A 65%", family = font_family, size = 3.4, color = clone_a) +
    ggplot2::annotate("text", x = 2.70, y = 4.25, label = "B 55%", family = font_family, size = 3.4, color = clone_b) +
    ggplot2::annotate("text", x = 2.5, y = 6.25, label = "65% + 55% = 120%", family = font_family, size = 3.7) +
    ggplot2::annotate("text", x = 2.5, y = 3.95, label = "separate branches cannot both fit", family = font_family, size = 3.2, color = red) +
    ggplot2::annotate("segment", x = 6.0, xend = 9.0, y = 4.55, yend = 4.55, linewidth = 0.35) +
    ggplot2::annotate("segment", x = 6.0, xend = 6.0, y = 4.55, yend = 6.05, linewidth = 0.35) +
    ggplot2::annotate("segment", x = 9.0, xend = 9.0, y = 4.55, yend = 6.05, linewidth = 0.35) +
    ggplot2::annotate("segment", x = 6.0, xend = 9.0, y = 5.90, yend = 4.80, linewidth = 1.2, color = clone_a) +
    ggplot2::annotate("segment", x = 6.0, xend = 9.0, y = 4.80, yend = 5.90, linewidth = 1.2, color = clone_b) +
    ggplot2::annotate("point", x = c(6.0, 9.0), y = c(5.90, 4.80), size = 2.6, color = clone_a) +
    ggplot2::annotate("point", x = c(6.0, 9.0), y = c(4.80, 5.90), size = 2.6, color = clone_b) +
    ggplot2::annotate("text", x = 6.0, y = 4.25, label = "S1", family = font_family, size = 3.4) +
    ggplot2::annotate("text", x = 9.0, y = 4.25, label = "S2", family = font_family, size = 3.4) +
    ggplot2::annotate("text", x = 7.5, y = 6.25, label = "A > B in S1, but A < B in S2", family = font_family, size = 3.5) +
    ggplot2::annotate("text", x = 7.5, y = 3.95, label = "ordering flips, so they must branch", family = font_family, size = 3.2, color = red) +
    ggplot2::coord_cartesian(xlim = c(0, 10), ylim = c(0, 7.1), clip = "off") +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = grid::unit(c(8, 8, 8, 8), "pt")
    )

  ct_save_plot(p, filename, width, height)
}
