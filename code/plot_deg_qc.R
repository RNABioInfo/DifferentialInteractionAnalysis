plot_concordance <- function(analysis) {
  if (is.null(analysis$concordance)) {
    plot.new()
    title("Method concordance unavailable")
    text(0.5, 0.5, "Run analysis_method: both to compute concordance")
    return(invisible(NULL))
  }
  ggplot2::ggplot(
    analysis$concordance,
    ggplot2::aes(x = edgeR_log2FoldChange, y = DESeq2_log2FoldChange, color = concordance_class)
  ) +
    ggplot2::geom_point(alpha = 0.5, size = 1.35) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::labs(x = "edgeR log2 fold-change", y = "DESeq2 log2 fold-change", title = "Method concordance")
}
