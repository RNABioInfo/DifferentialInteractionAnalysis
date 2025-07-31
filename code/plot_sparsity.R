plot_sparsity = function(obj, ...) {
  UseMethod("plot_sparsity")
}

plot_sparsity.DESeq2Analysis = function(obj) {
  plotSparsity(obj$deseq_dataset)
}

plot_sparsity.EdgeRAnalysis = function(obj) {
  plotSparsity(obj$norm_counts)
}