AnalysisParameters <- function(params) {
  structure(
    .Data = list (
      out_dir = params$out_dir,
      analysis_method = params$analysis_method,
      metadata_path = params$metadata_file,
      interaction_counts_path = params$counts_file,
      interaction_regions_path = params$interactions_file,
      annotations_path = params$annotations_file,
      min_inter_contr = params$min_inter_contr,
      min_n_samples_inter_contr = params$min_n_samples_contr,
      min_n_samples_greater_zero = params$min_n_samples_greater_zero,
      features_ids_of_interest = params$feature_ids,
      padj_threshold = params$padj_thresh
    ),
    class = "AnalysisParameters"
  )
}