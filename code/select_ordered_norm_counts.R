select_ordered_norm_counts = function(obj, n = 30, sort_by, zscore = TRUE) {
  UseMethod("select_ordered_norm_counts")
}

select_ordered_norm_counts.default = function(obj, n = 30, sorted_by, zscore = TRUE) {
  results = sig_results(obj)
  end = min(n, nrow(results))
  
  results = results[1:end,]
  
  results = results[rev(order(results[[sorted_by]])),]
  
  norm_counts = obj$norm_counts[rownames(results),]
  
  if (zscore) {
    names = colnames(norm_counts)
    
    norm_counts = t(apply(norm_counts, 1, scale))
    colnames(norm_counts) = names
  }
  
  return(norm_counts)
}

select_ordered_norm_counts.EdgeRAnalysis = function(obj, n = 30, sorted_by = "logFC", zscore = TRUE) {
  return(select_ordered_norm_counts.default(obj, n, sorted_by, zscore))
}

select_ordered_norm_counts.DESeq2Analysis = function(obj, n = 30, sorted_by = "log2FoldChange", zscore = TRUE) {
  return(select_ordered_norm_counts.default(obj, n, sorted_by, zscore))
}