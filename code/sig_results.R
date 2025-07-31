sig_results = function(obj, ...) {
  UseMethod("sig_results")
}

sig_results.EdgeRAnalysis = function(obj, sorted_by = "FDR") {
  res = as.data.frame(topTags(obj$edgeR_results, n = Inf, p.value = obj$params$padj_threshold))
  
  if (nrow(res) < 2) {
    return(res)
  }
  
  res = res[order(res[[sorted_by]]),]
  
  return(res)
}

sig_results.DESeq2Analysis = function(obj, sorted_by = "padj") {
  res = obj$deseq_results[obj$deseq_results$pvalue <= obj$params$padj_threshold,]
  res = res[order(res[[sorted_by]]),]
  
  
  return(as.data.frame(res))
}
