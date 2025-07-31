plot_deg = function(obj, ...) {
  UseMethod("plot_deg")
}

plot_deg.EdgeRAnalysis = function(obj) {}

plot_deg.DESeq2Analysis = function(obj) {
  degPlot(dds = obj$deseq_dataset, res = obj$deseq_results, n =20, xs = "condition")
}