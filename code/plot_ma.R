plot_ma = function(obj, ...) {
  UseMethod("plot_ma")
}

plot_ma.EdgeRAnalysis = function (obj) {
  dataset = as.DEGSet(topTags(obj$edgeR_results, n = Inf, p.value = 1 ))
  
  colnames(dataset$raw) = c("log2FoldChange", "baseMean", "F", "pvalue", "padj")
  dataset$raw = subset(dataset$raw, select = -c(F))
  
  DEGreport::degMA(dataset, correlation = TRUE, diff = 0)
}

plot_ma.DESeq2Analysis = function(obj) {
  DEGreport::degMA(as.DEGSet(obj$deseq_results), correlation = TRUE)
}