plot_method_volcano <- function(analysis, method) {
  res <- analysis$method_results[[method]]$annotated_results
  high_confidence_ids <- if (!is.null(analysis$high_confidence_results)) {
    analysis$high_confidence_results$interaction_id
  } else {
    character()
  }
  res$plot_class <- dplyr::case_when(
    res$interaction_id %in% high_confidence_ids ~ "high confidence",
    res$significant ~ "significant",
    TRUE ~ "not significant"
  )
  res$neg_log10_padj <- -log10(pmax(res$padj, .Machine$double.xmin))
  ggplot2::ggplot(res, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = plot_class)) +
    ggplot2::geom_point(alpha = 0.45, size = 1.2) +
    ggplot2::scale_color_manual(values = c(
      "high confidence" = "#C43B3B",
      "significant" = "#326C9C",
      "not significant" = "grey70"
    )) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "log2 fold-change", y = "-log10(FDR)", color = NULL, title = paste(method, "volcano plot"))
}
