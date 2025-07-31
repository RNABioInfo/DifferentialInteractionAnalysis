plot_read_counts = function(obj, ...) {
  UseMethod("plot_read_counts")
}

plot_read_counts.default = function(counts_df, title) {
  counts_plot_melt = reshape2::melt(counts_df, c("interaction", "sample"), value.name = "count")
  counts_plot_melt$sample = as.factor(counts_plot_melt$sample)
  
  counts_plot = ggplot(counts_plot_melt, aes(x = sample, y = count)) + 
    geom_violin() +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
    geom_boxplot(width=.1) + 
    labs(title = title, y = "Counts", x = "Sample") +
    theme_light()
  
  print(counts_plot)
}

plot_read_counts.DESeq2Analysis = function(obj) {
  counts = counts(obj$deseq_dataset, normalized = FALSE)
  
  plot_read_counts.default(counts, "Raw Gene Counts")
  
  plot_read_counts.default(obj$norm_counts, "Norm. Gene Counts")
}

plot_read_counts.EdgeRAnalysis = function(obj) {
  plot_read_counts.default(analysis$edgeR_dataset$counts, "Raw Gene Counts")
  
  counts_norm = cpm(obj$norm_counts, normalized.lib.sizes = TRUE)
  
  plot_read_counts.default(obj$norm_counts, "Norm. Gene Counts")
}
