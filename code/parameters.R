`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) {
    return(y)
  }
  x
}

as_character_vector <- function(x, default = character()) {
  if (is.null(x)) {
    return(default)
  }
  if (is.list(x)) {
    return(unlist(x, use.names = FALSE))
  }
  x
}

sanitize_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

is_present <- function(x) {
  !is.null(x) && length(x) > 0 && !is.na(x) && nzchar(x)
}

resolve_stage_dir <- function(root_dir, stage_name, explicit_dir = NULL) {
  if (is_present(explicit_dir)) {
    return(explicit_dir)
  }
  if (!is_present(root_dir)) {
    return(NA_character_)
  }

  unlabelled <- file.path(root_dir, stage_name)
  if (dir.exists(unlabelled)) {
    return(unlabelled)
  }

  stage_matches <- list.files(
    root_dir,
    pattern = paste0("^", stage_name, "($|_)"),
    full.names = TRUE
  )
  stage_matches <- stage_matches[dir.exists(stage_matches)]
  if (length(stage_matches) == 1) {
    return(stage_matches[1])
  }

  unlabelled
}

resolve_rnanue_result_file <- function(rnanue_results_dir, postprocess_dir, filename) {
  candidates <- character()
  if (is_present(postprocess_dir)) {
    candidates <- c(candidates, file.path(postprocess_dir, filename))
  }
  if (is_present(rnanue_results_dir)) {
    candidates <- c(candidates, file.path(rnanue_results_dir, filename))
  }

  existing_candidates <- unique(candidates[file.exists(candidates)])
  if (length(existing_candidates) > 0) {
    return(existing_candidates[1])
  }

  if (is_present(rnanue_results_dir) && dir.exists(rnanue_results_dir)) {
    matches <- list.files(
      rnanue_results_dir,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    matches <- unique(matches[file.exists(matches)])
    postprocess_matches <- matches[grepl("^05_postprocess($|_)", basename(dirname(matches)))]

    if (length(postprocess_matches) == 1) {
      return(postprocess_matches[1])
    }
    if (length(matches) == 1) {
      return(matches[1])
    }
    if (length(matches) > 1) {
      stop(
        "Could not infer ", filename, " from rnanue_results_dir because multiple matches were found. ",
        "Set postprocess_dir explicitly. Matches: ",
        paste(matches, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (is_present(postprocess_dir)) {
    return(file.path(postprocess_dir, filename))
  }
  stop("rnanue_results_dir or postprocess_dir must be set to discover ", filename, ".", call. = FALSE)
}

reject_retired_input_params <- function(params) {
  retired_params <- c("counts_file", "interactions_file")
  present_params <- retired_params[vapply(retired_params, function(name) !is.null(params[[name]]), logical(1))]
  if (length(present_params) > 0) {
    stop(
      "Retired input parameter(s) are no longer supported: ",
      paste(present_params, collapse = ", "),
      ". Use rnanue_results_dir and optional stage-level overrides instead.",
      call. = FALSE
    )
  }

  annotation <- params$annotation %||% list()
  if ("gff_file" %in% names(annotation)) {
    stop(
      "Retired parameter annotation.gff_file is no longer supported. ",
      "Use annotation.target_annotations_overwrite for an optional final annotation override; ",
      "annotations_file must point to the exact GFF used by RNAnue.",
      call. = FALSE
    )
  }
}

AnalysisParameters <- function(params) {
  reject_retired_input_params(params)

  rnanue_results_dir <- params$rnanue_results_dir %||% NA_character_

  postprocess_dir <- resolve_stage_dir(
    rnanue_results_dir,
    "05_postprocess",
    explicit_dir = params$postprocess_dir
  )

  detect_dir <- resolve_stage_dir(
    rnanue_results_dir,
    "03_detect",
    explicit_dir = params$detect_dir
  )

  analyze_dir <- resolve_stage_dir(
    rnanue_results_dir,
    "04_analyze",
    explicit_dir = params$analyze_dir
  )

  primary_contrast <- params$primary_contrast %||% list()
  normalization <- params$normalization %||% list()
  count_handling <- params$count_handling %||% list()
  high_confidence <- params$high_confidence %||% list()
  batch_correction <- params$batch_correction %||% list()
  annotation <- params$annotation %||% list()
  legacy_feature_ids <- params$feature_ids %||% ""
  target_annotations_path <-
    annotation$target_annotations_overwrite %||% params$annotations_file %||% NA_character_
  annotation_strand_policy <- tolower(annotation$strand_policy %||% "")
  if (!nzchar(annotation_strand_policy)) {
    annotation_strand_policy <- if (!is.null(annotation$same_strand)) {
      if (isTRUE(annotation$same_strand)) "same" else "both"
    } else {
      "both"
    }
  }
  batch_mode <- tolower(batch_correction$mode %||% "auto")
  batch_on_confounded <- tolower(batch_correction$on_confounded %||% "fail")
  if (!annotation_strand_policy %in% c("same", "opposite", "both")) {
    stop("annotation.strand_policy must be one of same, opposite, or both.", call. = FALSE)
  }
  if (!batch_mode %in% c("auto", "off", "required")) {
    stop("batch_correction.mode must be one of auto, off, or required.", call. = FALSE)
  }
  if (!batch_on_confounded %in% c("fail", "warn_skip")) {
    stop("batch_correction.on_confounded must be one of fail or warn_skip.", call. = FALSE)
  }

  structure(
    list(
      out_dir = params$out_dir %||% "results",
      analysis_method = params$analysis_method %||% "both",
      metadata_path = params$metadata_file %||% "data/metadata.csv",
      interaction_counts_path = resolve_rnanue_result_file(
        rnanue_results_dir,
        postprocess_dir,
        "complete_super_interaction_transcript_counts.gct"
      ),
      interaction_regions_path = resolve_rnanue_result_file(
        rnanue_results_dir,
        postprocess_dir,
        "complete_super_interaction_regions.bedpe"
      ),
      annotations_path = params$annotations_file %||% NA_character_,
      rnanue_results_dir = rnanue_results_dir,
      postprocess_dir = postprocess_dir,
      detect_dir = detect_dir,
      analyze_dir = analyze_dir,
      case_role = primary_contrast$case_role %||% "treatment",
      control_role = primary_contrast$control_role %||% "ligation_control",
      qc_roles = as_character_vector(params$qc_roles, "no_ligation_control"),
      sample_exposure = normalization$sample_exposure %||% "split_reads",
      pair_background = isTRUE(normalization$pair_background %||% TRUE),
      offset_pseudocount = as.numeric(normalization$offset_pseudocount %||% 0.5),
      batch_correction_mode = batch_mode,
      batch_columns = as_character_vector(batch_correction$columns, "batch"),
      batch_on_confounded = batch_on_confounded,
      batch_visualization_remove_batch = isTRUE(batch_correction$visualization_remove_batch %||% TRUE),
      edgeR_count_handling = count_handling$edgeR %||% "numeric",
      DESeq2_count_handling = count_handling$DESeq2 %||% "round",
      high_confidence_rnanue_padj_max = as.numeric(high_confidence$rnanue_padj_max %||% 0.1),
      annotation_enabled = isTRUE(annotation$enabled %||% TRUE),
      target_annotations_path = target_annotations_path,
      annotation_feature_types = as_character_vector(
        annotation$feature_types,
        c("gene", "sRNA", "tRNA", "rRNA", "ncRNA", "transcript", "TU", "CDS", "mobile_genetic_element")
      ),
      annotation_ignore_feature_types = as_character_vector(
        annotation$ignore_feature_types,
        c("region", "sequence_feature")
      ),
      annotation_match_attributes = as_character_vector(
        annotation$match_attributes,
        c("ID", "Name", "gene", "locus_tag", "Alias", "Parent")
      ),
      annotation_biotype_attributes = as_character_vector(
        annotation$biotype_attributes,
        c("gene_biotype", "biotype", "gbkey", "type")
      ),
      annotation_seqname_normalization = annotation$seqname_normalization %||% "auto_strip_version",
      annotation_strand_policy = annotation_strand_policy,
      annotation_same_strand = identical(annotation_strand_policy, "same"),
      annotation_min_overlap_bp = as.integer(annotation$min_overlap_bp %||% 1),
      annotation_target_ids = as_character_vector(annotation$target_ids, legacy_feature_ids),
      annotation_target_result_sets = as_character_vector(
        annotation$target_result_sets,
        c("high_confidence", "concordant_same_direction", "edgeR_significant", "DESeq2_significant")
      ),
      min_inter_contr = as.numeric(params$min_inter_contr %||% 10),
      min_n_samples_inter_contr = as.integer(params$min_n_samples_contr %||% 2),
      min_n_samples_greater_zero = as.integer(params$min_n_samples_greater_zero %||% 3),
      features_ids_of_interest = legacy_feature_ids,
      padj_threshold = as.numeric(params$padj_thresh %||% params$padj_threshold %||% 0.1)
    ),
    class = "AnalysisParameters"
  )
}
