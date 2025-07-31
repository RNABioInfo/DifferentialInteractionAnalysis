select_ordered_padj = function(obj, ...) {
  UseMethod("select_ordered_padj")
}

select_ordered_padj.EdgeRAnalysis = function(obj, n = 30, sorted_by = "logFC") {
  results = sig_results(obj)
  end = min(n, nrow(results))
  
  results = results[1:end,]
  
  results = results[rev(order(results[[sorted_by]])),]
  
  padj = as.data.frame(results$FDR)
  
  rownames(padj) = rownames(results)
  colnames(padj) = c("FDR")
  return(padj)
}

select_ordered_padj.DESeq2Analysis = function(obj, n = 30, sorted_by = "log2FoldChange") {
  results = sig_results(obj)
  end = min(n, nrow(results))
  
  results = results[1:end,]
  
  results = results[rev(order(results[[sorted_by]])),]
  
  padj = as.data.frame(results$padj)
  
  rownames(padj) = rownames(results)
  colnames(padj) = c("padj")
  return(padj)
}