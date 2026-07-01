empty_diagnostic_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title)
}

plot_pca_assessment <- function(analysis, pca_set) {
  scores <- analysis$postprocess$pca_scores[[pca_set]]
  if (is.null(scores) || nrow(scores) == 0) {
    return(empty_diagnostic_plot(
      paste("PCA:", pca_set),
      "PCA was skipped for this set."
    ))
  }
  scores$role <- factor(scores$role, levels = role_plot_levels())
  ggplot2::ggplot(
    scores,
    ggplot2::aes(x = PC1, y = PC2, color = role, shape = transform, label = sample_id)
  ) +
    ggplot2::geom_point(size = 4.5) +
    ggplot2::geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    ggplot2::facet_wrap(~ transform, scales = "free") +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = paste0("PC1 (", round(unique(scores$PC1_percent_variance)[1], 1), "%)"),
      y = paste0("PC2 (", round(unique(scores$PC2_percent_variance)[1], 1), "%)"),
      title = paste("Exploratory PCA:", gsub("_", " ", pca_set))
    )
}

plot_sample_correlation_heatmap <- function(analysis) {
  mat <- log2(primary_normalized_counts(analysis) + 1)
  sample_order <- ordered_samples_by_role(analysis$model_metadata, colnames(mat))
  mat <- mat[, sample_order, drop = FALSE]
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    return(empty_diagnostic_plot(
      "Sample correlation heatmap",
      "Need at least two samples and two interactions."
    ))
  }
  cor_mat <- stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("sample_1", "sample_2", "spearman_correlation")
  cor_df$sample_1 <- factor(cor_df$sample_1, levels = sample_order)
  cor_df$sample_2 <- factor(cor_df$sample_2, levels = rev(sample_order))
  ggplot2::ggplot(cor_df, ggplot2::aes(x = sample_1, y = sample_2, fill = spearman_correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), low = "#9A3A32", mid = "white", high = "#2D6F6D") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Spearman r", title = "Sample correlation on log2 offset-normalized counts")
}

plot_dia_rnanue_padj_comparison <- function(analysis) {
  comparison <- analysis$postprocess$dia_rnanue_padj_comparison
  if (is.null(comparison) || nrow(comparison) == 0) {
    return(empty_diagnostic_plot("DIA versus RNAnue padj", "No DIA/RNAnue padj comparison table available."))
  }
  plot_df <- comparison[
    is.finite(comparison$primary_dia_neg_log10_padj) &
      is.finite(comparison$rnanue_condition_neg_log10_padj_mean),
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0) {
    return(empty_diagnostic_plot(
      "DIA versus RNAnue padj",
      "No interactions had both primary DIA padj and RNAnue condition-sample padj values."
    ))
  }

  highlight_df <- plot_df[
    (!is.na(plot_df$result_class) & plot_df$result_class != "not_significant") |
      plot_df$high_confidence,
    ,
    drop = FALSE
  ]
  highlight_df$high_confidence_label <- ifelse(highlight_df$high_confidence, "yes", "no")
  method <- unique(plot_df$primary_method)[1]

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = primary_dia_neg_log10_padj,
      y = rnanue_condition_neg_log10_padj_mean
    )
  ) +
    ggplot2::geom_bin_2d(bins = 35) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65") +
    ggplot2::geom_point(
      data = highlight_df,
      ggplot2::aes(color = result_class, shape = high_confidence_label),
      size = 3,
      alpha = 0.9
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.2,
      label = paste0("n = ", nrow(plot_df)),
      size = 3.4
    ) +
    ggplot2::scale_fill_gradient(low = "#E8EEF2", high = "#305C7A") +
    ggplot2::theme_light() +
    ggplot2::scale_shape_manual(values = c(no = 16, yes = 17)) +
    ggplot2::labs(
      x = paste0("-log10(", method, " FDR)"),
      y = "Mean -log10 RNAnue padj across treatment samples",
      fill = "Interactions/bin",
      color = "Result class",
      shape = "High confidence",
      title = "Primary DIA FDR versus RNAnue condition-sample support",
      caption = "Density shows all complete rows; points highlight non-background result classes, including method-specific calls."
    )
}

plot_dia_rnanue_discordance <- function(analysis) {
  discordance <- analysis$postprocess$dia_rnanue_discordance
  if (is.null(discordance) || nrow(discordance) == 0) {
    return(empty_diagnostic_plot("DIA/RNAnue evidence classes", "No discordance table available."))
  }
  class_levels <- c("DIA+/RNAnue+", "DIA-only", "RNAnue-only", "neither")
  discordance$evidence_class <- factor(discordance$evidence_class, levels = class_levels)
  counts <- as.data.frame(table(discordance$evidence_class), stringsAsFactors = FALSE)
  names(counts) <- c("evidence_class", "interactions")
  counts$evidence_class <- factor(counts$evidence_class, levels = class_levels)

  ggplot2::ggplot(counts, ggplot2::aes(x = evidence_class, y = interactions, fill = evidence_class)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = interactions), vjust = -0.25, size = 3.4) +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = NULL,
      y = "Interactions",
      title = "DIA and RNAnue evidence classes"
    )
}

plot_partner_recurrence <- function(analysis, top_n = 20) {
  recurrence <- analysis$postprocess$partner_recurrence
  if (is.null(recurrence) || nrow(recurrence) == 0) {
    return(empty_diagnostic_plot("Partner recurrence", "No partner features available."))
  }
  plot_df <- head(recurrence[order(-recurrence$total_interactions), , drop = FALSE], top_n)
  plot_df$feature_id <- factor(plot_df$feature_id, levels = rev(plot_df$feature_id))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = feature_id, y = total_interactions, fill = high_confidence_interactions)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = NULL,
      y = "Interactions containing partner",
      fill = "High-confidence interactions",
      title = "Recurrent RNA partners"
    )
}

diagnostic_grid_message <- function(title, message) {
  grid::grid.newpage()
  grid::grid.text(
    paste(title, message, sep = "\n"),
    x = 0.5,
    y = 0.5,
    gp = grid::gpar(fontsize = 11)
  )
  invisible(NULL)
}

annotation_colors <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) {
    return(character())
  }
  stats::setNames(scales::hue_pal()(length(values)), values)
}

significant_set_ids <- function(analysis, set_name, top_n = Inf) {
  ids <- switch(
    set_name,
    edgeR_significant = {
      if ("edgeR" %in% names(analysis$method_results)) {
        res <- analysis$method_results$edgeR$annotated_results
        res$interaction_id[res$significant]
      } else {
        character()
      }
    },
    DESeq2_significant = {
      if ("DESeq2" %in% names(analysis$method_results)) {
        res <- analysis$method_results$DESeq2$annotated_results
        res$interaction_id[res$significant]
      } else {
        character()
      }
    },
    concordant_same_direction = {
      if (!is.null(analysis$concordance)) {
        analysis$concordance$interaction_id[
          analysis$concordance$concordance_class == "both_significant_same_direction"
        ]
      } else {
        character()
      }
    },
    high_confidence = {
      if (!is.null(analysis$high_confidence_results)) {
        analysis$high_confidence_results$interaction_id
      } else {
        character()
      }
    },
    called = called_interaction_ids(analysis),
    called_interaction_ids(analysis)
  )

  stability <- analysis$postprocess$candidate_stability
  order_df <- stability[match(ids, stability$interaction_id), , drop = FALSE]
  ids <- order_df$interaction_id[order(order_df$padj, order_df$interaction_id, na.last = TRUE)]
  ids <- ids[!is.na(ids)]
  if (is.finite(top_n)) {
    ids <- head(ids, top_n)
  }
  ids
}

draw_significant_complex_heatmap <- function(analysis, set_name, top_n = 80) {
  ids <- significant_set_ids(analysis, set_name, top_n = top_n)
  normalized <- primary_normalized_counts(analysis)
  sample_order <- ordered_samples_by_role(analysis$model_metadata, colnames(normalized))
  normalized <- normalized[, sample_order, drop = FALSE]
  ids <- intersect(ids, rownames(normalized))
  if (length(ids) < 2 || ncol(normalized) < 2) {
    return(diagnostic_grid_message(
      paste("Complex heatmap:", gsub("_", " ", set_name)),
      "Need at least two interactions and two samples."
    ))
  }

  mat <- log2(normalized[ids, , drop = FALSE] + 1)
  z_mat <- t(scale(t(mat)))
  z_mat[!is.finite(z_mat)] <- 0

  row_data <- analysis$postprocess$candidate_stability[
    match(rownames(z_mat), analysis$postprocess$candidate_stability$interaction_id),
    ,
    drop = FALSE
  ]
  row_labels <- paste0(
    rownames(z_mat),
    " | LFC=",
    round(row_data$log2FoldChange, 2)
  )

  column_data <- analysis$model_metadata[colnames(z_mat), , drop = FALSE]
  column_anno_df <- data.frame(
    role = column_data$role,
    condition = column_data$condition,
    stringsAsFactors = FALSE
  )
  for (batch_column in analysis$batch_columns_used) {
    column_anno_df[[batch_column]] <- as.character(column_data[[batch_column]])
  }
  column_colors <- lapply(column_anno_df, annotation_colors)

  row_anno_df <- data.frame(
    class = row_data$result_class,
    high_confidence = ifelse(row_data$high_confidence, "yes", "no"),
    stringsAsFactors = FALSE
  )
  row_colors <- list(
    class = annotation_colors(row_anno_df$class),
    high_confidence = c("yes" = "#C43B3B", "no" = "grey80")
  )

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = column_anno_df,
    col = column_colors,
    annotation_name_side = "left"
  )
  left_annotation <- ComplexHeatmap::rowAnnotation(
    df = row_anno_df,
    col = row_colors,
    annotation_name_side = "bottom"
  )

  heatmap <- ComplexHeatmap::Heatmap(
    z_mat,
    name = "row z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2D6F6D", "white", "#9A3A32")),
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    row_labels = row_labels,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = length(ids) <= 80,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = paste("Significant interaction heatmap:", gsub("_", " ", set_name)),
    heatmap_legend_param = list(title = "z")
  )
  grid::grid.newpage()
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(NULL)
}

plot_significant_normalized_counts <- function(analysis, set_name = "high_confidence", top_n = 20) {
  counts <- analysis$postprocess$significant_normalized_counts
  if (is.null(counts) || nrow(counts) == 0) {
    return(empty_diagnostic_plot("Normalized counts", "No significant normalized counts available."))
  }
  ids <- significant_set_ids(analysis, set_name, top_n = top_n)
  plot_df <- counts[counts$interaction_id %in% ids, , drop = FALSE]
  if (nrow(plot_df) == 0) {
    return(empty_diagnostic_plot(
      paste("Normalized counts:", gsub("_", " ", set_name)),
      "No interactions available for this set."
    ))
  }
  order_df <- unique(plot_df[, c("interaction_id", "padj"), drop = FALSE])
  order_df <- order_df[order(order_df$padj, order_df$interaction_id, na.last = TRUE), , drop = FALSE]
  plot_df$interaction_id <- factor(plot_df$interaction_id, levels = order_df$interaction_id)
  plot_df$contrast_group <- factor(
    plot_df$contrast_group,
    levels = levels(analysis$model_metadata$contrast_group)
  )
  plot_df$sample_id <- factor(plot_df$sample_id, levels = ordered_samples_by_role(analysis$model_metadata, plot_df$sample_id))
  plot_df$role <- factor(plot_df$role, levels = role_plot_levels())

  ggplot2::ggplot(
    plot_df,
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
    ggplot2::facet_wrap(~ interaction_id, scales = "free_y") +
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
      title = paste("Normalized counts for", gsub("_", " ", set_name), "interactions")
    )
}
