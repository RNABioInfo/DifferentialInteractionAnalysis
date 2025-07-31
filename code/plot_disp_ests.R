plot_disp_ests = function(obj, ...) {
  UseMethod("plot_disp_ests")
}

plot_disp_ests.EdgeRAnalysis = function(obj, ...) {
  plotQLDisp(obj$edgeR_fit, main = "Dispersion Estimates")
}

plot_disp_ests.DESeq2Analysis = function(obj, ...) {
  plotDispEsts(obj$deseq_dataset, main = "Dispersion Estimates")
}