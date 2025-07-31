annotate_sig_interactions = function(obj, ...) {
  UseMethod("annotate_sig_interactions")
}

annotate_sig_interactions.EdgeRAnalysis = function(obj, annotation_gff) {
  result_summary = annotate_interactions(obj$sig_interaction_bed, annotation_gff)
  
  result_summary$log2FoldChange = sig_results_df$logFC
  result_summary$padj = sig_results_df$FDR
  
  return(result_summary)
}

annotate_sig_interactions.DESeq2Analysis = function(obj, annotation_gff) {
  result_summary = annotate_interactions(obj$sig_interaction_bed, annotation_gff)
  
  result_summary$log2FoldChange = sig_results_df$log2FoldChange
  result_summary$padj = sig_results_df$padj
  
  return(result_summary)
}