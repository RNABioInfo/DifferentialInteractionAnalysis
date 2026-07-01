plot_offset_diagnostics <- function(analysis) {
  diag <- analysis$normalization_diagnostics
  availability <- as.data.frame(table(diag$pair_background_available_fraction), stringsAsFactors = FALSE)
  names(availability) <- c("pair_background_available_fraction", "interactions")
  availability$pair_background_available_fraction <- factor(
    availability$pair_background_available_fraction,
    levels = availability$pair_background_available_fraction[order(as.numeric(availability$pair_background_available_fraction))]
  )
  p1 <- ggplot2::ggplot(availability, ggplot2::aes(x = pair_background_available_fraction, y = interactions)) +
    ggplot2::geom_col(fill = "#496A81") +
    ggplot2::geom_text(ggplot2::aes(label = interactions), vjust = -0.25, size = 3) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Fraction of model samples with pair background", y = "Interactions")
  p2 <- ggplot2::ggplot(diag, ggplot2::aes(x = max_log_offset - min_log_offset)) +
    ggplot2::geom_histogram(bins = 30, fill = "#B66D3A", color = "white") +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Log-offset range across model samples", y = "Interactions")
  print(p1)
  print(p2)
}
