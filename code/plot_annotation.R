plot_partner_type_confusion <- function(analysis, result_set = "high_confidence") {
  confusion <- analysis$postprocess$gff_annotation$partner_type_confusion
  if (is.null(confusion) || nrow(confusion) == 0) {
    return(empty_diagnostic_plot("Partner-type confusion", "No GFF partner-type summary available."))
  }
  plot_df <- confusion[confusion$result_set == result_set, , drop = FALSE]
  if (nrow(plot_df) == 0) {
    return(empty_diagnostic_plot(
      paste("Partner-type confusion:", gsub("_", " ", result_set)),
      "No interactions available for this result set."
    ))
  }
  levels <- partner_class_order()
  plot_df$partner_class_1 <- factor(plot_df$partner_class_1, levels = levels)
  plot_df$partner_class_2 <- factor(plot_df$partner_class_2, levels = rev(levels))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = partner_class_1, y = partner_class_2, fill = interactions)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = interactions), size = 3) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(low = "white", high = "#496A81") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = "Partner class",
      y = "Partner class",
      fill = "Interactions",
      title = paste("GFF partner-type matrix:", gsub("_", " ", result_set))
    )
}

plot_target_gene_normalized_counts <- function(analysis, top_n = 24) {
  targets <- analysis$postprocess$gff_annotation$target_gene_interactions
  if (is.null(targets) || nrow(targets) == 0) {
    return(empty_diagnostic_plot("Target-gene normalized counts", "No configured target IDs matched called interactions."))
  }
  normalized <- primary_normalized_counts(analysis)
  ids <- unique(targets$interaction_id)
  ids <- intersect(ids, rownames(normalized))
  if (length(ids) == 0) {
    return(empty_diagnostic_plot("Target-gene normalized counts", "No target interactions are present in normalized counts."))
  }
  target_order <- unique(targets$target_id[order(targets$padj, targets$target_id, na.last = TRUE)])
  ids <- head(ids, top_n)
  long_counts <- as.data.frame(as.table(normalized[ids, , drop = FALSE]), stringsAsFactors = FALSE)
  names(long_counts) <- c("interaction_id", "sample_id", "normalized_count")
  long_counts$log2_normalized_count <- log2(long_counts$normalized_count + 1)
  long_counts <- dplyr::left_join(long_counts, analysis$model_metadata, by = "sample_id")
  target_labels <- targets[, c("interaction_id", "target_id", "other_arm_partner_class"), drop = FALSE]
  target_labels <- target_labels[!duplicated(target_labels$interaction_id), , drop = FALSE]
  target_labels$facet_label <- paste0(
    target_labels$target_id,
    " | ",
    target_labels$interaction_id,
    " | partner: ",
    target_labels$other_arm_partner_class
  )
  long_counts <- dplyr::left_join(long_counts, target_labels, by = "interaction_id")
  long_counts$sample_id <- factor(long_counts$sample_id, levels = ordered_samples_by_role(analysis$model_metadata, long_counts$sample_id))
  long_counts$role <- factor(long_counts$role, levels = role_plot_levels())
  long_counts$target_id <- factor(long_counts$target_id, levels = target_order)
  ggplot2::ggplot(
    long_counts,
    ggplot2::aes(x = contrast_group, y = log2_normalized_count, color = role, shape = sample_id)
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      size = 3.3,
      alpha = 0.9
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = contrast_group),
      fun.data = function(y) {
        median_y <- stats::median(y, na.rm = TRUE)
        data.frame(y = median_y, ymin = median_y, ymax = median_y)
      },
      geom = "crossbar",
      width = 0.45,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::facet_wrap(~ facet_label, scales = "free_y") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 7)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "log2 offset-normalized count + 1",
      color = "Role",
      shape = "Sample",
      title = "Normalized counts for configured target-gene interactions"
    )
}
