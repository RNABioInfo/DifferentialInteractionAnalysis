plot_sparsity_overview <- function(analysis) {
  sparsity <- data.frame(
    sample_id = colnames(analysis$counts_all),
    nonzero_interactions = colSums(analysis$counts_all > 0),
    zero_fraction = colMeans(analysis$counts_all == 0),
    stringsAsFactors = FALSE
  )
  sparsity <- dplyr::left_join(sparsity, analysis$metadata, by = "sample_id")
  sparsity$sample_id <- factor(sparsity$sample_id, levels = ordered_samples_by_role(analysis$metadata, sparsity$sample_id))
  sparsity$role <- factor(sparsity$role, levels = role_plot_levels())
  ggplot2::ggplot(sparsity, ggplot2::aes(x = sample_id, y = zero_fraction, fill = role)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Sample", y = "Fraction of zero-count interactions", title = "Interaction matrix sparsity")
}

plot_pca_counts <- function(analysis) {
  mat <- log2(as.matrix(analysis$counts_model) + 1)
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    plot.new()
    title("PCA unavailable")
    text(0.5, 0.5, "Need at least two samples and two interactions")
    return(invisible(NULL))
  }
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pca_df <- dplyr::left_join(pca_df, analysis$model_metadata, by = "sample_id")
  pca_df$role <- factor(pca_df$role, levels = role_plot_levels())
  ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = role, label = sample_id)) +
    ggplot2::geom_point(size = 4.5) +
    ggplot2::geom_text(vjust = -0.8, size = 3) +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      title = "Model-sample PCA on log2 raw interaction counts"
    )
}
