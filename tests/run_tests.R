source("code/load_libraries.R")
source_code_files <- setdiff(
  list.files("code", pattern = "\\.R$", full.names = TRUE),
  c("code/load_libraries.R", "code/check_environment.R")
)
invisible(lapply(source_code_files, source))

stopifnot(strand_policy_matches("+", "+", "same"))
stopifnot(!strand_policy_matches("+", "-", "same"))
stopifnot(strand_policy_matches("+", "-", "opposite"))
stopifnot(!strand_policy_matches("+", "+", "opposite"))
stopifnot(strand_policy_matches("+", "-", "both"))

sample_specific_features <- paste0(
  "c1:11111111-1111-1111-1111-111111111111;",
  "c2:22222222-2222-2222-2222-222222222222;"
)
stopifnot(identical(
  parse_feature_field_for_sample(sample_specific_features, "c1"),
  "11111111-1111-1111-1111-111111111111"
))
stopifnot(length(parse_feature_field_for_sample(sample_specific_features, "t1")) == 0)
stopifnot(identical(
  parse_feature_field_for_sample("c1:geneA;c2:geneA;", "t1", "geneA"),
  "geneA"
))

stable_uuid <- "aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa"
supplementary_uuid <- "bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb"
stable_id_gff <- tempfile(fileext = ".gff3")
writeLines(
  c(
    "##gff-version 3",
    paste("chr1", "test", "gene", 1, 10, ".", "+", ".", paste0("ID=", stable_uuid), sep = "\t"),
    paste("chr1", "test", "transcript", 20, 30, ".", "+", ".", "ID=child;Parent=parent1,parent2", sep = "\t"),
    paste("chr1", "test", "gene", 40, 50, ".", "+", ".", "ID=geneB", sep = "\t")
  ),
  stable_id_gff
)
stable_id_result <- read_rnanue_annotation_ids(stable_id_gff)
stopifnot(stable_id_result$available)
stopifnot(all(c(stable_uuid, "child", "parent1", "parent2", "geneB") %in% stable_id_result$ids))

mixed_feature_map <- parse_sample_feature_map(paste0(
  "c1:", stable_uuid, ",", supplementary_uuid, ";",
  "c2:", stable_uuid, ";"
))
stopifnot(setequal(
  feature_ids_for_sample(mixed_feature_map, "c1", stable_id_result$ids),
  c(stable_uuid, supplementary_uuid)
))
stopifnot(identical(
  feature_ids_for_sample(mixed_feature_map, "t1", stable_id_result$ids),
  stable_uuid
))
stopifnot(length(feature_ids_for_sample(
  parse_sample_feature_map(paste0("c1:", supplementary_uuid, ";")),
  "t1",
  stable_id_result$ids
)) == 0)

offset_test_counts <- matrix(
  1,
  nrow = 1,
  ncol = 4,
  dimnames = list("stable_uuid_interaction", c("c1", "c2", "t1", "t2"))
)
offset_test_bedpe <- data.frame(
  interaction_id = "stable_uuid_interaction",
  arm1_features = paste0("c1:", stable_uuid, ";c2:", stable_uuid, ";"),
  arm2_features = "c1:geneB;c2:geneB;",
  stringsAsFactors = FALSE,
  row.names = "stable_uuid_interaction"
)
offset_test_reads <- data.frame(
  sample_id = colnames(offset_test_counts),
  splits = rep(100, 4),
  stringsAsFactors = FALSE
)
offset_test_transcript_counts <- list(
  c1 = stats::setNames(c(10, 20), c(stable_uuid, "geneB")),
  c2 = stats::setNames(c(12, 22), c(stable_uuid, "geneB")),
  t1 = stats::setNames(c(25, 75), c("geneB", "other")),
  t2 = stats::setNames(c(30, 70), c("geneB", "other"))
)
stable_uuid_offsets <- build_offset_matrices(
  offset_test_counts,
  offset_test_bedpe,
  offset_test_reads,
  offset_test_transcript_counts,
  list(offset_pseudocount = 0.5, pair_background = TRUE, annotations_path = stable_id_gff)
)
stopifnot(stable_uuid_offsets$normalization_diagnostics$pair_background_available_fraction == 1)
stopifnot(stable_uuid_offsets$normalization_setup_diagnostics$stable_feature_id_count == 5)
stopifnot(stable_uuid_offsets$normalization_setup_diagnostics$pair_background_applied)

supplementary_only_gff <- tempfile(fileext = ".gff3")
writeLines(
  paste("chr1", "test", "gene", 40, 50, ".", "+", ".", "ID=geneB", sep = "\t"),
  supplementary_only_gff
)
supplementary_offsets <- build_offset_matrices(
  offset_test_counts,
  offset_test_bedpe,
  offset_test_reads,
  offset_test_transcript_counts,
  list(offset_pseudocount = 0.5, pair_background = TRUE, annotations_path = supplementary_only_gff)
)
stopifnot(supplementary_offsets$normalization_diagnostics$pair_background_available_fraction == 0.5)

missing_gff_warning <- NULL
missing_gff_offsets <- withCallingHandlers(
  build_offset_matrices(
    offset_test_counts,
    offset_test_bedpe,
    offset_test_reads,
    offset_test_transcript_counts,
    list(offset_pseudocount = 0.5, pair_background = TRUE, annotations_path = tempfile())
  ),
  warning = function(w) {
    missing_gff_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(grepl("Exposure-only offsets", missing_gff_warning))
stopifnot(!missing_gff_offsets$normalization_setup_diagnostics$pair_background_applied)
stopifnot(all(missing_gff_offsets$log_offset == log(100)))

malformed_gff <- tempfile(fileext = ".gff3")
writeLines("not\ta\tvalid\tgff", malformed_gff)
stopifnot(!read_rnanue_annotation_ids(malformed_gff)$available)
identifier_empty_gff <- tempfile(fileext = ".gff3")
writeLines(paste("chr1", "test", "gene", 1, 10, ".", "+", ".", "Name=no_id", sep = "\t"), identifier_empty_gff)
stopifnot(!read_rnanue_annotation_ids(identifier_empty_gff)$available)

multi_cluster_tokens <- parse_cluster_tokens("s1:c1,c2;s2:c3;")
stopifnot(identical(multi_cluster_tokens$sample_id, c("s1", "s1", "s2")))
stopifnot(identical(multi_cluster_tokens$cluster_ID, c("c1", "c2", "c3")))
single_cluster_tokens <- parse_cluster_tokens("s1:c1;s2:c2;")
stopifnot(identical(single_cluster_tokens$sample_id, c("s1", "s2")))
stopifnot(identical(single_cluster_tokens$cluster_ID, c("c1", "c2")))
blank_cluster_tokens <- parse_cluster_tokens("s1:c1, ,c2;;")
stopifnot(identical(blank_cluster_tokens$sample_id, c("s1", "s1")))
stopifnot(identical(blank_cluster_tokens$cluster_ID, c("c1", "c2")))

warning_capture_probe <- capture_analysis_warnings({
  warning("synthetic warning for capture test", call. = FALSE)
  42L
}, stage = "synthetic_stage", method = "synthetic_method")
stopifnot(identical(warning_capture_probe$value, 42L))
stopifnot(nrow(warning_capture_probe$warnings) == 1)
stopifnot(warning_capture_probe$warnings$message == "synthetic warning for capture test")

fixture_dir <- "tests/fixtures/synthetic_rnanue"
out_dir <- tempfile("dia_test_")

params <- list(
  out_dir = out_dir,
  analysis_method = "both",
  metadata_file = file.path(fixture_dir, "metadata.csv"),
  annotations_file = file.path(fixture_dir, "synthetic_annotations.gff3"),
  rnanue_results_dir = fixture_dir,
  primary_contrast = list(case_role = "treatment", control_role = "ligation_control"),
  qc_roles = c("no_ligation_control"),
  normalization = list(sample_exposure = "split_reads", pair_background = TRUE, offset_pseudocount = 0.5),
  count_handling = list(edgeR = "numeric", DESeq2 = "round"),
  annotation = list(
    enabled = TRUE,
    target_annotations_overwrite = file.path(fixture_dir, "synthetic_annotations.gff3"),
    feature_types = c("gene", "sRNA", "tRNA", "TU"),
    ignore_feature_types = c("region", "sequence_feature"),
    target_ids = c("geneA", "geneB", "targetA"),
    seqname_normalization = "auto_strip_version",
    strand_policy = "both",
    min_overlap_bp = 1
  ),
  min_inter_contr = 1,
  min_n_samples_contr = 2,
  min_n_samples_greater_zero = 2,
  feature_ids = "geneA,geneB",
  padj_thresh = 0.25
)

retired_param_error <- tryCatch(
  {
    AnalysisParameters(c(params, list(counts_file = file.path(fixture_dir, "old.gct"))))
    NULL
  },
  error = function(e) e$message
)
stopifnot(!is.null(retired_param_error))
stopifnot(grepl("Retired input parameter", retired_param_error))

retired_annotation_error <- tryCatch(
  {
    legacy_params <- params
    legacy_params$annotation$gff_file <- legacy_params$annotation$target_annotations_overwrite
    legacy_params$annotation$target_annotations_overwrite <- NULL
    AnalysisParameters(legacy_params)
    NULL
  },
  error = function(e) e$message
)
stopifnot(!is.null(retired_annotation_error))
stopifnot(grepl("annotation.target_annotations_overwrite", retired_annotation_error, fixed = TRUE))

retired_null_annotation_params <- params
retired_null_annotation_params$annotation["gff_file"] <- list(NULL)
retired_null_annotation_error <- tryCatch(
  {
    AnalysisParameters(retired_null_annotation_params)
    NULL
  },
  error = function(e) e$message
)
stopifnot(!is.null(retired_null_annotation_error))

default_annotation_params <- params
default_annotation_params$annotation$target_annotations_overwrite <- NULL
stopifnot(identical(
  AnalysisParameters(default_annotation_params)$target_annotations_path,
  params$annotations_file
))

make_temp_postprocess_dir <- function(counts_path, bedpe_path) {
  postprocess_dir <- tempfile("dia_test_postprocess_")
  dir.create(postprocess_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(
    counts_path,
    file.path(postprocess_dir, "complete_super_interaction_transcript_counts.gct")
  )
  file.copy(
    bedpe_path,
    file.path(postprocess_dir, "complete_super_interaction_regions.bedpe")
  )
  postprocess_dir
}

discovery_analysis <- make_analysis(params)
stopifnot(identical(
  normalizePath(discovery_analysis$params$interaction_counts_path),
  normalizePath(file.path(fixture_dir, "complete_super_interaction_transcript_counts.gct"))
))
stopifnot(identical(
  normalizePath(discovery_analysis$params$interaction_regions_path),
  normalizePath(file.path(fixture_dir, "complete_super_interaction_regions.bedpe"))
))

analysis <- make_analysis(params)
stopifnot(inherits(analysis, "InteractionAnalysis"))
stopifnot(setequal(analysis$model_metadata$sample_id, c("c1", "c2", "t1", "t2")))
stopifnot(!"n1" %in% analysis$model_metadata$sample_id)
stopifnot(nrow(analysis$no_ligation_qc) == 1)
stopifnot(analysis$design_diagnostics$batch_status == "no_usable_batch")
stopifnot(analysis$design_diagnostics$design_full_rank)
stopifnot(all(dim(analysis$log_offset) == dim(analysis$counts_model)))
stopifnot(analysis$normalization_diagnostics$pair_background_available_fraction[analysis$normalization_diagnostics$interaction_id == "si0"] == 1)
stopifnot(analysis$normalization_diagnostics$pair_background_available_fraction[analysis$normalization_diagnostics$interaction_id == "si3"] == 0.5)
stopifnot(analysis$normalization_diagnostics$unresolved_arm1[analysis$normalization_diagnostics$interaction_id == "si3"])
stopifnot(analysis$normalization_setup_diagnostics$pair_background_requested)
stopifnot(analysis$normalization_setup_diagnostics$pair_background_applied)
stopifnot(analysis$normalization_setup_diagnostics$stable_feature_id_count > 0)
stopifnot("11111111-1111-1111-1111-111111111111" %in% names(analysis$interaction_transcript_counts$c1))
stopifnot(is.null(analysis$contiguous_counts))
stopifnot(all(c("mapped_read_count", "mapped_read_source", "interaction_library_sum") %in% names(analysis$sample_qc)))
stopifnot(all(
  analysis$sample_qc$mapped_read_count ==
    analysis$sample_qc$splits + analysis$sample_qc$singletons
))
stopifnot(all(analysis$sample_qc$mapped_read_source == "rnanue_splits_plus_singletons"))
expected_interaction_sums <- colSums(read_counts(analysis$params$interaction_counts_path)[, analysis$sample_qc$sample_id, drop = FALSE])
stopifnot(all(analysis$sample_qc$interaction_library_sum == unname(expected_interaction_sums)))
stopifnot(!all(analysis$sample_qc$interaction_library_sum == analysis$sample_qc$mapped_read_count))

alternate_target_gff <- tempfile(fileext = ".gff3")
writeLines(
  paste("chr1", "test", "gene", 1, 10, ".", "+", ".", "ID=alternate_target", sep = "\t"),
  alternate_target_gff
)
alternate_annotation_params <- params
alternate_annotation_params$annotation$target_annotations_overwrite <- alternate_target_gff
alternate_annotation_analysis <- make_analysis(alternate_annotation_params)
stopifnot(identical(alternate_annotation_analysis$log_offset, analysis$log_offset))
stopifnot(identical(
  alternate_annotation_analysis$params$target_annotations_path,
  alternate_target_gff
))
stopifnot(identical(analysis$params$annotations_path, params$annotations_file))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines("stale", file.path(out_dir, "requested_feature_dia_results.tsv"))
writeLines("stale", file.path(out_dir, "requested_feature_normalized_counts.tsv"))
writeLines("stale", file.path(out_dir, "analysis_warnings.tsv"))
writeLines("stale", file.path(out_dir, "dia_rnanue_correlation_stats.tsv"))
writeLines("stale", file.path(out_dir, "dia_rnanue_threshold_enrichment.tsv"))
analysis <- run_analysis(analysis)
stopifnot(file.exists(file.path(out_dir, "design_diagnostics.tsv")))
stopifnot(file.exists(file.path(out_dir, "batch_balance.tsv")))
stopifnot(file.exists(file.path(out_dir, "analysis_warnings.tsv")))
stopifnot(file.exists(file.path(out_dir, "normalization_setup_diagnostics.tsv")))
stopifnot(readLines(file.path(out_dir, "analysis_warnings.tsv"), n = 1) == "warning_index\tstage\tmethod\tmessage")
stopifnot(all(c("edgeR", "DESeq2") %in% names(analysis$method_results)))
stopifnot(identical(
  analysis$method_results$edgeR$filter_keep,
  custom_count_filter(analysis$counts_model, analysis$params)
))
stopifnot(all(analysis$method_results$edgeR$dataset$samples$norm.factors == 1))
stopifnot(file.exists(file.path(out_dir, "edgeR_results.tsv")))
stopifnot(file.exists(file.path(out_dir, "DESeq2_results.tsv")))
stopifnot(file.exists(file.path(out_dir, "method_concordance.tsv")))
stopifnot(file.exists(file.path(out_dir, "candidate_stability.tsv")))
stopifnot(file.exists(file.path(out_dir, "significant_normalized_counts.tsv")))
stopifnot(file.exists(file.path(out_dir, "dia_rnanue_padj_comparison.tsv")))
stopifnot(file.exists(file.path(out_dir, "dia_rnanue_discordance.tsv")))
stopifnot(file.exists(file.path(out_dir, "rnanue_metric_group_summary.tsv")))
stopifnot(file.exists(file.path(out_dir, "partner_recurrence.tsv")))
stopifnot(file.exists(file.path(out_dir, "pca_diagnostics.tsv")))
stopifnot(file.exists(file.path(out_dir, "pca_scores_all_model_filtered.tsv")))
stopifnot(file.exists(file.path(out_dir, "high_confidence_interactions.tsv")))
stopifnot(file.exists(file.path(out_dir, "high_confidence_interactions.bedpe")))
stopifnot(file.exists(file.path(out_dir, "high_confidence_interaction_regions.bed")))
stopifnot(file.exists(file.path(out_dir, "high_confidence_interaction_regions.tsv")))
stopifnot(file.exists(file.path(out_dir, "DESeq2_rounding_diagnostics.tsv")))
stopifnot(file.exists(file.path(out_dir, "interaction_arm_gff_annotations.tsv")))
stopifnot(file.exists(file.path(out_dir, "interaction_partner_type_pairs.tsv")))
stopifnot(file.exists(file.path(out_dir, "partner_type_confusion.tsv")))
stopifnot(file.exists(file.path(out_dir, "target_gene_interactions.tsv")))
stopifnot(!file.exists(file.path(out_dir, "requested_feature_dia_results.tsv")))
stopifnot(!file.exists(file.path(out_dir, "requested_feature_normalized_counts.tsv")))
stopifnot(!file.exists(file.path(out_dir, "dia_rnanue_correlation_stats.tsv")))
stopifnot(!file.exists(file.path(out_dir, "dia_rnanue_threshold_enrichment.tsv")))
stopifnot(file.exists(file.path(out_dir, "target_gene_interactions.bedpe")))
stopifnot(file.exists(file.path(out_dir, "gff_annotation_diagnostics.tsv")))
stopifnot(nrow(analysis$concordance) > 0)
stopifnot(!is.null(analysis$high_confidence_results))
high_conf <- readr::read_tsv(file.path(out_dir, "high_confidence_interactions.tsv"), show_col_types = FALSE)
high_conf_bedpe <- readr::read_tsv(file.path(out_dir, "high_confidence_interactions.bedpe"), show_col_types = FALSE)
high_conf_region_bed <- readr::read_tsv(
  file.path(out_dir, "high_confidence_interaction_regions.bed"),
  col_names = c("chrom", "start", "end", "name", "score", "strand"),
  show_col_types = FALSE
)
high_conf_region_tsv <- readr::read_tsv(file.path(out_dir, "high_confidence_interaction_regions.tsv"), show_col_types = FALSE)
candidate_stability <- readr::read_tsv(file.path(out_dir, "candidate_stability.tsv"), show_col_types = FALSE)
significant_normalized_counts <- readr::read_tsv(file.path(out_dir, "significant_normalized_counts.tsv"), show_col_types = FALSE)
dia_rnanue_padj_comparison <- readr::read_tsv(file.path(out_dir, "dia_rnanue_padj_comparison.tsv"), show_col_types = FALSE)
dia_rnanue_discordance <- readr::read_tsv(file.path(out_dir, "dia_rnanue_discordance.tsv"), show_col_types = FALSE)
pca_diagnostics <- readr::read_tsv(file.path(out_dir, "pca_diagnostics.tsv"), show_col_types = FALSE)
sample_qc <- readr::read_tsv(file.path(out_dir, "sample_qc.tsv"), show_col_types = FALSE)
method_concordance <- readr::read_tsv(file.path(out_dir, "method_concordance.tsv"), show_col_types = FALSE)
arm_annotations <- readr::read_tsv(file.path(out_dir, "interaction_arm_gff_annotations.tsv"), show_col_types = FALSE)
partner_confusion <- readr::read_tsv(file.path(out_dir, "partner_type_confusion.tsv"), show_col_types = FALSE)
target_interactions <- readr::read_tsv(file.path(out_dir, "target_gene_interactions.tsv"), show_col_types = FALSE)
annotation_diagnostics <- readr::read_tsv(file.path(out_dir, "gff_annotation_diagnostics.tsv"), show_col_types = FALSE)
stopifnot(all(c("mapped_read_count", "mapped_read_source") %in% names(sample_qc)))
stopifnot(sample_qc$mapped_read_count[sample_qc$sample_id == "t2"] == 2700)
stopifnot(sample_qc$interaction_library_sum[sample_qc$sample_id == "t2"] == 102)
stopifnot(all(c("rnanue_condition_padj_values", "rnanue_control_padj_values") %in% names(method_concordance)))
si0_concordance <- method_concordance[method_concordance$interaction_id == "si0", , drop = FALSE]
stopifnot(nrow(si0_concordance) == 1)
stopifnot(grepl("^t1:.*;t2:", si0_concordance$rnanue_condition_padj_values))
stopifnot(grepl("^c1:.*;c2:", si0_concordance$rnanue_control_padj_values))
stopifnot(all(c(
  "interaction_id",
  "primary_method",
  "primary_dia_padj",
  "primary_dia_neg_log10_padj",
  "rnanue_condition_sample_padj_values",
  "rnanue_condition_padj_n",
  "rnanue_condition_neg_log10_padj_mean",
  "rnanue_condition_neg_log10_padj_sd",
  "result_class",
  "high_confidence"
) %in% names(dia_rnanue_padj_comparison)))
stopifnot(all(is.finite(dia_rnanue_padj_comparison$primary_dia_padj)))
stopifnot(all(c(
  "interaction_id",
  "evidence_class",
  "primary_dia_padj",
  "dia_positive",
  "rnanue_positive",
  "single_sample_dominance",
  "no_ligation_max_count"
) %in% names(dia_rnanue_discordance)))
stopifnot(nrow(dia_rnanue_discordance) == nrow(dia_rnanue_padj_comparison))
si0_padj_comparison <- dia_rnanue_padj_comparison[
  dia_rnanue_padj_comparison$interaction_id == "si0",
  ,
  drop = FALSE
]
stopifnot(nrow(si0_padj_comparison) == 1)
stopifnot(si0_padj_comparison$rnanue_condition_padj_n == 2)
stopifnot(grepl("t1:0.006;t2:0.006", si0_padj_comparison$rnanue_condition_sample_padj_values, fixed = TRUE))
stopifnot(is.finite(si0_padj_comparison$primary_dia_neg_log10_padj))
stopifnot(is.finite(si0_padj_comparison$rnanue_condition_neg_log10_padj_mean))

primary_method <- primary_result_method(analysis)
padj_mutation_id <- analysis$method_results[[primary_method]]$annotated_results$interaction_id[1]
na_padj_analysis <- analysis
na_padj_analysis$method_results[[primary_method]]$annotated_results$padj[1] <- NA_real_
na_padj_comparison <- make_dia_rnanue_padj_comparison(
  na_padj_analysis,
  na_padj_analysis$postprocess$result_classes
)
na_padj_discordance <- make_dia_rnanue_discordance(
  na_padj_comparison,
  na_padj_analysis$postprocess$candidate_stability,
  na_padj_analysis$params
)
stopifnot(!padj_mutation_id %in% na_padj_comparison$interaction_id)
stopifnot(!padj_mutation_id %in% na_padj_discordance$interaction_id)
stopifnot(all(is.finite(na_padj_comparison$primary_dia_padj)))

zero_padj_analysis <- analysis
zero_padj_analysis$method_results[[primary_method]]$annotated_results$padj[1] <- 0
zero_padj_comparison <- make_dia_rnanue_padj_comparison(
  zero_padj_analysis,
  zero_padj_analysis$postprocess$result_classes
)
zero_padj_row <- zero_padj_comparison[
  zero_padj_comparison$interaction_id == padj_mutation_id,
  ,
  drop = FALSE
]
stopifnot(nrow(zero_padj_row) == 1)
stopifnot(zero_padj_row$primary_dia_padj == 0)
stopifnot(is.finite(zero_padj_row$primary_dia_neg_log10_padj))

multi_cluster_bedpe_lines <- readLines(file.path(fixture_dir, "complete_super_interaction_regions.bedpe"), warn = FALSE)
multi_cluster_bedpe_lines <- gsub(
  "c1:c0;c2:c0;t1:c0;t2:c0;",
  "c1:c0;c2:c0;t1:c0,c1;t2:c0;",
  multi_cluster_bedpe_lines,
  fixed = TRUE
)
multi_cluster_bedpe_path <- tempfile("multi_cluster_rnanue_", fileext = ".bedpe")
writeLines(multi_cluster_bedpe_lines, multi_cluster_bedpe_path)
multi_cluster_params <- params
multi_cluster_params$out_dir <- tempfile("dia_test_multi_cluster_")
multi_cluster_params$postprocess_dir <- make_temp_postprocess_dir(
  analysis$params$interaction_counts_path,
  multi_cluster_bedpe_path
)
multi_cluster_analysis <- make_analysis(multi_cluster_params)
multi_cluster_analysis <- run_analysis(multi_cluster_analysis)
multi_cluster_concordance <- readr::read_tsv(
  file.path(multi_cluster_params$out_dir, "method_concordance.tsv"),
  show_col_types = FALSE
)
multi_cluster_padj_comparison <- readr::read_tsv(
  file.path(multi_cluster_params$out_dir, "dia_rnanue_padj_comparison.tsv"),
  show_col_types = FALSE
)
multi_cluster_si0 <- multi_cluster_concordance[
  multi_cluster_concordance$interaction_id == "si0",
  ,
  drop = FALSE
]
stopifnot(nrow(multi_cluster_si0) == 1)
stopifnot(!is.na(multi_cluster_si0$rnanue_condition_padj_values))
stopifnot(grepl("t1:0.006,0.04", multi_cluster_si0$rnanue_condition_padj_values, fixed = TRUE))
stopifnot(grepl(";t2:0.006", multi_cluster_si0$rnanue_condition_padj_values, fixed = TRUE))
multi_cluster_si0_padj_comparison <- multi_cluster_padj_comparison[
  multi_cluster_padj_comparison$interaction_id == "si0",
  ,
  drop = FALSE
]
stopifnot(nrow(multi_cluster_si0_padj_comparison) == 1)
stopifnot(grepl("t1:0.006;t2:0.006", multi_cluster_si0_padj_comparison$rnanue_condition_sample_padj_values, fixed = TRUE))
stopifnot(all(c("high_confidence", "DESeq2_log2FoldChange", "pair_background_available_fraction") %in% names(high_conf)))

relaxed_high_conf_analysis <- analysis
relaxed_high_conf_id <- relaxed_high_conf_analysis$method_results$edgeR$annotated_results$interaction_id[1]
edgeR_index <- match(
  relaxed_high_conf_id,
  relaxed_high_conf_analysis$method_results$edgeR$annotated_results$interaction_id
)
DESeq2_index <- match(
  relaxed_high_conf_id,
  relaxed_high_conf_analysis$method_results$DESeq2$annotated_results$interaction_id
)
normalization_index <- match(
  relaxed_high_conf_id,
  relaxed_high_conf_analysis$normalization_diagnostics$interaction_id
)
relaxed_high_conf_analysis$method_results$edgeR$annotated_results$significant[edgeR_index] <- TRUE
relaxed_high_conf_analysis$method_results$edgeR$annotated_results$log2FoldChange[edgeR_index] <- 1
relaxed_high_conf_analysis$method_results$edgeR$annotated_results$rnanue_min_padj_value[edgeR_index] <- 0.01
relaxed_high_conf_analysis$method_results$edgeR$annotated_results$no_ligation_max_count[edgeR_index] <- 1e6
relaxed_high_conf_analysis$method_results$edgeR$annotated_results$rnanue_coverage_profiles[edgeR_index] <- "broad_diffuse"
relaxed_high_conf_analysis$method_results$DESeq2$annotated_results$padj[DESeq2_index] <- 0.01
relaxed_high_conf_analysis$method_results$DESeq2$annotated_results$log2FoldChange[DESeq2_index] <- 1
relaxed_high_conf_analysis$normalization_diagnostics$pair_background_available_fraction[normalization_index] <- 0
relaxed_high_conf_results <- make_high_confidence_results(relaxed_high_conf_analysis)
stopifnot(relaxed_high_conf_id %in% relaxed_high_conf_results$interaction_id)

stopifnot(all(c("chr1", "start1", "end1", "chr2", "start2", "end2", "interaction_id") %in% names(high_conf_bedpe)))
stopifnot(all(c("chrom", "start", "end", "name", "score", "strand") %in% names(high_conf_region_tsv)))
stopifnot(nrow(high_conf_region_bed) == 2 * nrow(high_conf))
stopifnot(nrow(high_conf_region_tsv) == 2 * nrow(high_conf))
stopifnot(all(c("single_sample_dominance", "min_treatment_over_max_control", "result_class") %in% names(candidate_stability)))
stopifnot(all(c("interaction_id", "sample_id", "normalized_count", "log2_normalized_count", "result_class") %in% names(significant_normalized_counts)))
stopifnot(nrow(significant_normalized_counts) > 0)
stopifnot(all(c("pca_set", "transform", "status", "n_interactions") %in% names(pca_diagnostics)))
stopifnot(all(c(
  "selected_feature_id",
  "selected_feature_type",
  "selected_partner_class",
  "feature_annotation_source",
  "partner_class_source",
  "reciprocal_overlap_fraction"
) %in% names(arm_annotations)))
called_ids <- called_interaction_ids(analysis)
stopifnot(setequal(unique(arm_annotations$interaction_id), called_ids))
stopifnot(nrow(arm_annotations) == 2 * length(called_ids))
stopifnot(!any(is.na(arm_annotations$result_class) | arm_annotations$result_class == "not_significant"))
stopifnot(arm_annotations$selected_feature_id[
  arm_annotations$interaction_id == "si0" & arm_annotations$arm == "arm1"
] == "geneA")
stopifnot(arm_annotations$selected_partner_class[
  arm_annotations$interaction_id == "si0" & arm_annotations$arm == "arm2"
] == "sRNA")
stopifnot(arm_annotations$partner_class_source[
  arm_annotations$interaction_id == "si0" & arm_annotations$arm == "arm2"
] == "gff_overlap")
stopifnot(arm_annotations$selected_feature_id[
  arm_annotations$interaction_id == "si5" & arm_annotations$arm == "arm1"
] == "geneOpp")
stopifnot(arm_annotations$selected_partner_class[
  arm_annotations$interaction_id == "si5" & arm_annotations$arm == "arm1"
] == "rRNA")
stopifnot(arm_annotations$selected_partner_class[
  arm_annotations$interaction_id == "si5" & arm_annotations$arm == "arm2"
] == "unannotated")
stopifnot(arm_annotations$partner_class_source[
  arm_annotations$interaction_id == "si5" & arm_annotations$arm == "arm2"
] == "unannotated_no_gff_overlap")
stopifnot(any(arm_annotations$seqname_normalized[arm_annotations$selected_feature_id == "geneF"]))
stopifnot(all(c("result_set", "partner_class_1", "partner_class_2", "interactions") %in% names(partner_confusion)))
stopifnot(all(c("target_id", "interaction_id", "target_arm", "other_arm_partner_class") %in% names(target_interactions)))
stopifnot(any(target_interactions$target_id == "targetA"))
stopifnot(annotation_diagnostics$value[annotation_diagnostics$metric == "status"] == "computed")

mutated_bedpe_lines <- readLines(file.path(fixture_dir, "complete_super_interaction_regions.bedpe"), warn = FALSE)
mutated_bedpe_lines <- gsub("c1:geneB;c2:geneB;t1:geneB;t2:geneB;", "c1:proteinCodingLike;c2:proteinCodingLike;t1:proteinCodingLike;t2:proteinCodingLike;", mutated_bedpe_lines, fixed = TRUE)
mutated_bedpe_lines <- gsub("c1:11111111-1111-1111-1111-111111111111;c2:11111111-1111-1111-1111-111111111111;", "c1:proteinCodingLike;c2:proteinCodingLike;", mutated_bedpe_lines, fixed = TRUE)
mutated_bedpe_lines <- gsub("c1:geneG;c2:geneG;t1:geneG;t2:geneG;", "c1:proteinCodingLike;c2:proteinCodingLike;t1:proteinCodingLike;t2:proteinCodingLike;", mutated_bedpe_lines, fixed = TRUE)
mutated_bedpe_path <- tempfile("mutated_rnanue_ids_", fileext = ".bedpe")
writeLines(mutated_bedpe_lines, mutated_bedpe_path)
mutated_params <- params
mutated_params$out_dir <- tempfile("dia_test_mutated_rnanue_ids_")
mutated_params$postprocess_dir <- make_temp_postprocess_dir(
  analysis$params$interaction_counts_path,
  mutated_bedpe_path
)
mutated_analysis <- make_analysis(mutated_params)
mutated_arms <- interaction_arm_table(mutated_analysis, analysis$postprocess$result_classes)
mutated_arms <- mutated_arms[mutated_arms$interaction_id %in% called_ids, , drop = FALSE]
mutated_features <- prepare_gff_features(mutated_analysis$params, mutated_arms$chrom)
mutated_arm_annotations <- make_arm_annotations(mutated_arms, mutated_features$features, mutated_analysis$params)
mutated_partner_pairs <- make_partner_type_pairs(mutated_arm_annotations)
mutated_partner_confusion <- make_partner_type_confusion(analysis, mutated_partner_pairs)
classification_cols <- c(
  "interaction_id",
  "arm",
  "selected_feature_id",
  "selected_feature_type",
  "selected_partner_class",
  "feature_annotation_source",
  "partner_class_source"
)
normalize_missing_strings <- function(x) {
  x <- as.data.frame(x)
  for (column in names(x)) {
    if (is.character(x[[column]])) {
      x[[column]][is.na(x[[column]])] <- ""
    }
  }
  x
}
base_classification <- arm_annotations[, classification_cols, drop = FALSE]
mutated_classification <- mutated_arm_annotations[, classification_cols, drop = FALSE]
base_classification <- base_classification[order(base_classification$interaction_id, base_classification$arm), , drop = FALSE]
mutated_classification <- mutated_classification[order(mutated_classification$interaction_id, mutated_classification$arm), , drop = FALSE]
rownames(base_classification) <- NULL
rownames(mutated_classification) <- NULL
stopifnot(isTRUE(all.equal(
  normalize_missing_strings(base_classification),
  normalize_missing_strings(mutated_classification),
  check.attributes = FALSE
)))
confusion_cols <- c("result_set", "partner_class_1", "partner_class_2", "interactions")
base_confusion <- partner_confusion[, confusion_cols, drop = FALSE]
mutated_confusion <- mutated_partner_confusion[, confusion_cols, drop = FALSE]
base_confusion <- base_confusion[order(base_confusion$result_set, base_confusion$partner_class_1, base_confusion$partner_class_2), , drop = FALSE]
mutated_confusion <- mutated_confusion[order(mutated_confusion$result_set, mutated_confusion$partner_class_1, mutated_confusion$partner_class_2), , drop = FALSE]
rownames(base_confusion) <- NULL
rownames(mutated_confusion) <- NULL
stopifnot(isTRUE(all.equal(
  as.data.frame(base_confusion),
  as.data.frame(mutated_confusion),
  check.attributes = FALSE
)))

edgeR_rows_first <- nrow(readr::read_tsv(file.path(out_dir, "edgeR_results.tsv"), show_col_types = FALSE))
analysis <- run_analysis(analysis)
edgeR_rows_second <- nrow(readr::read_tsv(file.path(out_dir, "edgeR_results.tsv"), show_col_types = FALSE))
stopifnot(edgeR_rows_first == edgeR_rows_second)
stopifnot(file.exists(file.path(out_dir, "high_confidence_interactions.bedpe")))
stopifnot(file.exists(file.path(out_dir, "high_confidence_interaction_regions.bed")))

metadata_no_batch <- readr::read_csv(file.path(fixture_dir, "metadata.csv"), show_col_types = FALSE)
metadata_balanced_batch <- metadata_no_batch
metadata_balanced_batch$batch <- c("b1", "b2", "b1", "b2", "b1")
balanced_metadata_path <- tempfile("metadata_balanced_batch_", fileext = ".csv")
readr::write_csv(metadata_balanced_batch, balanced_metadata_path)
balanced_params <- params
balanced_params$out_dir <- tempfile("dia_test_batch_")
balanced_params$metadata_file <- balanced_metadata_path
balanced_params$batch_correction <- list(
  mode = "auto",
  columns = c("batch"),
  on_confounded = "fail",
  visualization_remove_batch = TRUE
)
balanced_analysis <- make_analysis(balanced_params)
stopifnot(grepl("batch", paste(deparse(balanced_analysis$design_formula), collapse = "")))
stopifnot(balanced_analysis$design_diagnostics$batch_status == "batch_used")
stopifnot(balanced_analysis$design_diagnostics$design_full_rank)
balanced_analysis <- run_analysis(balanced_analysis)
stopifnot(file.exists(file.path(balanced_params$out_dir, "pca_scores_all_model_filtered.tsv")))

metadata_confounded_batch <- metadata_no_batch
metadata_confounded_batch$batch <- ifelse(
  metadata_confounded_batch$role == "treatment",
  "b2",
  "b1"
)
confounded_metadata_path <- tempfile("metadata_confounded_batch_", fileext = ".csv")
readr::write_csv(metadata_confounded_batch, confounded_metadata_path)
confounded_params <- params
confounded_params$out_dir <- tempfile("dia_test_confounded_")
confounded_params$metadata_file <- confounded_metadata_path
confounded_params$batch_correction <- list(
  mode = "auto",
  columns = c("batch"),
  on_confounded = "fail",
  visualization_remove_batch = TRUE
)
confounded_error <- tryCatch(
  {
    make_analysis(confounded_params)
    NULL
  },
  error = function(e) e$message
)
stopifnot(!is.null(confounded_error))
stopifnot(grepl("Design is not full rank", confounded_error))
stopifnot(file.exists(file.path(confounded_params$out_dir, "design_diagnostics.tsv")))
stopifnot(file.exists(file.path(confounded_params$out_dir, "batch_balance.tsv")))

off_params <- confounded_params
off_params$out_dir <- tempfile("dia_test_batch_off_")
off_params$batch_correction$mode <- "off"
off_analysis <- make_analysis(off_params)
stopifnot(off_analysis$design_diagnostics$batch_status == "batch_off")
stopifnot(!grepl("batch", paste(deparse(off_analysis$design_formula), collapse = "")))

message("Tests passed. Output directory: ", out_dir)
