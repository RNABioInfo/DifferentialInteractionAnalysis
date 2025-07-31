plot_volcano = function(obj, x, y, title) {
  UseMethod("plot_volcano")
}

plot_volcano.EdgeRAnalysis = function(obj) {
  sample_count = summary(obj$metadata_df$condition)
  condition_levels = names(sample_count)
  
  title = str_glue("{condition_levels[1]} vs {condition_levels[2]}")
  sample_string = sample_replicate_str(obj)
  
  print(
    EnhancedVolcano(obj$results,
                    lab = rownames(obj$results),
                    x = 'logFC',
                    y = 'FDR',
                    title = title,
                    subtitle = sample_string,
                    pCutoff = obj$params$padj_thresh)
  )
}

plot_volcano.DESeq2Analysis = function(obj) {
  sample_count = summary(obj$metadata_df$condition)
  condition_levels = names(sample_count)
  
  title = str_glue("{condition_levels[1]} vs {condition_levels[2]}")
  sample_string = sample_replicate_str(obj)
  
  print(
    EnhancedVolcano(obj$results,
                    lab = rownames(obj$results),
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = title,
                    subtitle = sample_string,
                    pCutoff = obj$params$padj_thresh)
  )
}