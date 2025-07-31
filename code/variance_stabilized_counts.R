variance_stabilized_counts = function(obj) {
  UseMethod("variance_stabilized_counts")
}

variance_stabilized_counts.EdgeRAnalysis = function(obj) {
  sig_res_df = sig_results(obj)
  
  vs_df = vst(obj$norm_counts, blind = FALSE)[rownames(sig_res_df),]
  
  return(vs_df)
}

variance_stabilized_counts.DESeq2Analysis = function(obj) {
  sig_res_df = sig_results(obj)
  
  vs_df = vst(obj$norm_counts, blind = FALSE)[rownames(sig_res_df),]
  
  return(vs_df)
}