EdgeRAnalysis <- function(params, metadata_df) {
  structure(
    .Data = list (
      params = AnalysisParameters(params),
      metadata_df = metadata_df,
      edgeR_dataset = NULL,
      edgeR_fit = NULL,
      edgeR_results = NULL,
      results = NULL,
      norm_counts = NULL,
      sig_interaction_bed = NULL
    ),
    class = "EdgeRAnalysis"
  )
}


get_filtered_dataset.EdgeRAnalysis <- function(obj) {
  counts_df = read_counts(obj)
  
  print(summary(counts_df))
  
  dge_list = DGEList(counts = counts_df, group = obj$metadata_df$condition)
  
  keep_min_count = rowSums(dge_list$counts >= obj$params$min_inter_contr) >= obj$params$min_n_samples_inter_contr
  
  keep_larger_zero = rowSums(dge_list$counts > 0) >= obj$params$min_n_samples_greater_zero
  
  combined_keep = keep_min_count & keep_larger_zero
  
  print("Number of removed and kept interactions")
  print(summary(combined_keep))
  
  obj$edgeR_dataset = dge_list[combined_keep,]
  
  return(obj)
  }