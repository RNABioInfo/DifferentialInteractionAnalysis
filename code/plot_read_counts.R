sample_qc_plots <- function(analysis) {
  qc <- analysis$sample_qc
  qc$sample_id <- factor(qc$sample_id, levels = ordered_samples_by_role(analysis$metadata, qc$sample_id))
  qc$role <- factor(qc$role, levels = role_plot_levels())
  mapped_read_count <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = mapped_read_count, fill = role)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = "Sample",
      y = "Mapped reads",
      title = "Mapped/read-classified reads by sample"
    )

  interaction_library_size <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = interaction_library_sum, fill = role)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Sample", y = "Interaction count sum", title = "Interaction library size")

  mapped_vs_interaction_sum <- ggplot2::ggplot(qc, ggplot2::aes(
    x = mapped_read_count,
    y = interaction_library_sum,
    color = role,
    label = sample_id
  )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    ggplot2::theme_light() +
    ggplot2::labs(
      x = "Mapped reads",
      y = "Interaction count sum",
      title = "Mapped reads vs interaction count sum"
    )

  split_fraction <- ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = split_fraction, fill = role)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Sample", y = "Split fraction", title = "RNAnue split-read fraction by sample")

  list(
    mapped_read_count = mapped_read_count,
    interaction_library_size = interaction_library_size,
    mapped_vs_interaction_sum = mapped_vs_interaction_sum,
    split_fraction = split_fraction
  )
}

plot_sample_qc <- function(analysis) {
  plots <- sample_qc_plots(analysis)
  print(plots$mapped_read_count)
  print(plots$interaction_library_size)
  print(plots$mapped_vs_interaction_sum)
  print(plots$split_fraction)
}
