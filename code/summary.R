summary.DESeq2Analysis = function(obj) {
  summary(obj$deseq_results)
}

summary.EdgeRAnalysis = function(obj) {
  summary(decideTests(obj$edgeR_results, p.value = obj$params$padj_thresh))
}