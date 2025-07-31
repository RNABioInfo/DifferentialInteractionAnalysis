DESeq2Analysis <- function(params, metadata_df) {
  structure(
    .Data = list (
      params = AnalysisParameters(params),
      metadata_df = metadata_df,
      deseq_dataset = NULL,
      deseq_results = NULL,
      results = NULL,
      norm_counts = NULL,
      sig_interaction_bed = NULL
    ),
    class = "DESeq2Analysis"
  )
}


get_filtered_dataset.DESeq2Analysis <- function(obj) {
  counts_df = read_counts(obj)
  
  print(summary(counts_df))
  
  # Round counts
  counts_df = counts_df %>%
    mutate(across(
      .cols = where(is.numeric),
      .fns  = ~ as.integer(round(.))
    ))
  
  dds = DESeqDataSetFromMatrix(counts_df, colData = obj$metadata_df, design = ~condition)
  
  keep_min_count = rowSums(counts(dds) >= obj$params$min_inter_contr) >= obj$params$min_n_samples_inter_contr
  
  keep_larger_zero = rowSums(counts(dds) > 0) >= obj$params$min_n_samples_greater_zero
  
  combined_keep = keep_min_count & keep_larger_zero
  
  print("Number of removed and kept interactions")
  print(summary(combined_keep))
  
  obj$deseq_dataset = dds[combined_keep,]
  
  return(obj)
}