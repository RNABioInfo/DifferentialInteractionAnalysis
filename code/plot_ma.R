plot_method_ma <- function(analysis, method) {
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
  abundance <- if ("logCPM" %in% names(res)) {
    res$logCPM
  } else {
    log2(res$baseMean + 1)
  }
  ggplot2::ggplot(res, ggplot2::aes(x = abundance, y = log2FoldChange, color = plot_class)) +
    ggplot2::geom_point(alpha = 0.45, size = 1.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c(
      "high confidence" = "#C43B3B",
      "significant" = "#326C9C",
      "not significant" = "grey70"
    )) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Average abundance", y = "log2 fold-change", color = NULL, title = paste(method, "MA plot"))
}
