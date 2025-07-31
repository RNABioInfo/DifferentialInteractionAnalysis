plot_deg_qc = function(obj, ...) {
  UseMethod("plot_deg_qc")
}

plot_deg_qc.EdgeRAnalysis = function(obj) {
  sorted_norm_counts = obj$norm_counts[rownames(obj$results), ]
  
  print(degQC(sorted_norm_counts, obj$metadata_df$condition, pvalue = obj$results$PValue))
  
  print(degMV(counts = sorted_norm_counts, group = obj$metadata_df$condition, pvalues = obj$results$PValue))
}

plot_deg_qc.DESeq2Analysis = function(obj) {
  degQC(obj$norm_counts, obj$metadata_df$condition, pvalue = obj$results$pvalue)
}

#plotMD(analysis$edgeR_results)