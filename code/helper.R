role_plot_levels <- function() {
  c("ligation_control", "treatment", "no_ligation_control", "other_qc")
}

ordered_samples_by_role <- function(metadata, sample_ids = NULL) {
  if (is.null(sample_ids)) {
    sample_ids <- metadata$sample_id
  }
  sample_ids <- as.character(sample_ids)
  sample_metadata <- metadata[metadata$sample_id %in% sample_ids, , drop = FALSE]
  sample_metadata$role <- factor(sample_metadata$role, levels = role_plot_levels())
  sample_metadata <- sample_metadata[order(sample_metadata$role, sample_metadata$condition, sample_metadata$sample_id), , drop = FALSE]
  c(sample_metadata$sample_id, setdiff(sample_ids, sample_metadata$sample_id))
}

batch_value_is_present <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

build_batch_balance <- function(model_metadata, requested_columns, used_columns = character()) {
  empty <- data.frame(
    sample_id = character(),
    condition = character(),
    role = character(),
    contrast_group = character(),
    batch_column = character(),
    batch_value = character(),
    used_in_model = logical(),
    stringsAsFactors = FALSE
  )
  available <- requested_columns[requested_columns %in% names(model_metadata)]
  if (length(available) == 0) {
    return(empty)
  }

  rows <- lapply(available, function(column) {
    data.frame(
      sample_id = model_metadata$sample_id,
      condition = model_metadata$condition,
      role = model_metadata$role,
      contrast_group = as.character(model_metadata$contrast_group),
      batch_column = column,
      batch_value = as.character(model_metadata[[column]]),
      used_in_model = column %in% used_columns,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

design_diagnostics_table <- function(params, design_formula, design, contrast_coef,
                                     requested_columns, used_columns, status, message) {
  data.frame(
    design_formula = paste(deparse(design_formula), collapse = ""),
    batch_mode = params$batch_correction_mode,
    requested_batch_columns = paste(requested_columns, collapse = ","),
    used_batch_columns = paste(used_columns, collapse = ","),
    batch_status = status,
    batch_message = message,
    visualization_remove_batch = params$batch_visualization_remove_batch,
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    design_full_rank = qr(design)$rank == ncol(design),
    contrast_coef = paste(contrast_coef, collapse = ","),
    stringsAsFactors = FALSE
  )
}

write_failed_design_diagnostics <- function(params, diagnostics, batch_balance) {
  dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(diagnostics, file.path(params$out_dir, "design_diagnostics.tsv"))
  readr::write_tsv(batch_balance, file.path(params$out_dir, "batch_balance.tsv"))
}

make_analysis_design <- function(model_metadata, params) {
  requested_columns <- unique(params$batch_columns)
  requested_columns <- requested_columns[!is.na(requested_columns) & nzchar(requested_columns)]
  used_columns <- character()
  dropped_columns <- character()
  status <- "batch_off"
  message <- "Batch correction disabled."

  if (params$batch_correction_mode != "off") {
    for (column in requested_columns) {
      if (!column %in% names(model_metadata)) {
        dropped_columns <- c(dropped_columns, paste0(column, ":missing"))
        next
      }
      values <- as.character(model_metadata[[column]])
      values[!batch_value_is_present(values)] <- NA_character_
      if (any(is.na(values))) {
        dropped_columns <- c(dropped_columns, paste0(column, ":missing_values"))
        next
      }
      if (length(unique(values)) < 2) {
        dropped_columns <- c(dropped_columns, paste0(column, ":single_level"))
        next
      }
      model_metadata[[column]] <- factor(values)
      used_columns <- c(used_columns, column)
    }

    if (length(used_columns) == 0) {
      if (params$batch_correction_mode == "required") {
        stop(
          "batch_correction.mode is 'required', but no requested batch column is usable. Dropped columns: ",
          paste(dropped_columns, collapse = ", "),
          call. = FALSE
        )
      }
      status <- "no_usable_batch"
      message <- if (length(dropped_columns) > 0) {
        paste("No usable batch columns; dropped", paste(dropped_columns, collapse = ", "))
      } else {
        "No requested batch columns were present."
      }
    } else {
      status <- "batch_used"
      message <- paste("Using batch columns:", paste(used_columns, collapse = ", "))
    }
  }

  terms <- c(used_columns, "contrast_group")
  design_formula <- stats::reformulate(terms)
  design <- stats::model.matrix(design_formula, data = model_metadata)
  contrast_coef <- grep("^contrast_group", colnames(design), value = TRUE)
  diagnostics <- design_diagnostics_table(
    params,
    design_formula,
    design,
    contrast_coef,
    requested_columns,
    used_columns,
    status,
    message
  )
  batch_balance <- build_batch_balance(model_metadata, requested_columns, used_columns)

  if (length(contrast_coef) != 1) {
    stop(
      "Could not identify a single treatment-vs-control contrast coefficient in design matrix.",
      call. = FALSE
    )
  }

  if (qr(design)$rank < ncol(design)) {
    diagnostics$batch_status <- "confounded_or_rank_deficient"
    diagnostics$batch_message <- paste0(
      "Design is not full rank. Batch columns are confounded with condition or another covariate. ",
      "Requested columns: ", paste(requested_columns, collapse = ", "),
      "; used columns: ", paste(used_columns, collapse = ", ")
    )
    if (params$batch_on_confounded == "warn_skip") {
      warning(diagnostics$batch_message, " Falling back to ~ contrast_group.", call. = FALSE)
      used_columns <- character()
      design_formula <- stats::reformulate("contrast_group")
      design <- stats::model.matrix(design_formula, data = model_metadata)
      contrast_coef <- grep("^contrast_group", colnames(design), value = TRUE)
      diagnostics <- design_diagnostics_table(
        params,
        design_formula,
        design,
        contrast_coef,
        requested_columns,
        used_columns,
        "batch_skipped_confounded",
        "Batch design was rank deficient; fallback design uses contrast_group only."
      )
      batch_balance <- build_batch_balance(model_metadata, requested_columns, used_columns)
    } else {
      write_failed_design_diagnostics(params, diagnostics, batch_balance)
      stop(
        diagnostics$batch_message,
        " See ",
        file.path(params$out_dir, "design_diagnostics.tsv"),
        " and ",
        file.path(params$out_dir, "batch_balance.tsv"),
        ".",
        call. = FALSE
      )
    }
  }

  list(
    model_metadata = model_metadata,
    design_formula = design_formula,
    design = design,
    contrast_coef = contrast_coef,
    batch_columns_used = used_columns,
    design_diagnostics = diagnostics,
    batch_balance = batch_balance
  )
}

make_analysis <- function(params) {
  params <- AnalysisParameters(params)
  metadata <- read_metadata(params$metadata_path, params)
  counts_all <- read_counts(params$interaction_counts_path)
  bedpe <- read_interactions(params$interaction_regions_path)

  model_metadata <- metadata[
    metadata$include != "exclude" &
      metadata$role %in% c(params$control_role, params$case_role),
    ,
    drop = FALSE
  ]

  if (nrow(model_metadata) == 0) {
    stop("No model samples remain after applying primary roles and include flags.", call. = FALSE)
  }

  role_counts <- table(model_metadata$role)
  for (role in c(params$control_role, params$case_role)) {
    if (!role %in% names(role_counts) || role_counts[[role]] < 2) {
      stop(
        "Primary contrast role '", role, "' has fewer than two model samples. ",
        "Replicated negative-binomial inference needs biological replication.",
        call. = FALSE
      )
    }
  }

  missing_model_samples <- setdiff(model_metadata$sample_id, colnames(counts_all))
  if (length(missing_model_samples) > 0) {
    stop(
      "Model samples missing from count table: ",
      paste(missing_model_samples, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- metadata[metadata$sample_id %in% colnames(counts_all), , drop = FALSE]

  case_condition <- unique(model_metadata$condition[model_metadata$role == params$case_role])
  control_condition <- unique(model_metadata$condition[model_metadata$role == params$control_role])
  if (length(case_condition) != 1 || length(control_condition) != 1) {
    stop("Each primary role must map to exactly one condition.", call. = FALSE)
  }

  model_metadata$contrast_group <- ifelse(
    model_metadata$role == params$case_role,
    case_condition,
    control_condition
  )
  model_metadata$contrast_group <- factor(
    model_metadata$contrast_group,
    levels = c(control_condition, case_condition)
  )
  rownames(model_metadata) <- model_metadata$sample_id
  design_info <- make_analysis_design(model_metadata, params)
  model_metadata <- design_info$model_metadata

  common_ids <- intersect(rownames(counts_all), bedpe$interaction_id)
  if (length(common_ids) == 0) {
    stop("No overlapping interaction IDs between GCT counts and BEDPE regions.", call. = FALSE)
  }

  counts_all <- counts_all[common_ids, , drop = FALSE]
  bedpe <- bedpe[common_ids, , drop = FALSE]
  counts_model <- counts_all[, model_metadata$sample_id, drop = FALSE]

  all_metadata_sample_ids <- metadata$sample_id
  read_summaries <- read_read_count_summaries(all_metadata_sample_ids, params$detect_dir)
  contiguous_counts <- read_contiguous_counts(model_metadata$sample_id, params$detect_dir)
  native_interactions <- read_native_interactions(all_metadata_sample_ids, params$analyze_dir)
  native_metrics <- build_rnanue_metric_summary(bedpe, native_interactions)
  native_padj_by_role <- build_rnanue_padj_by_role(bedpe, native_interactions, metadata, params)
  sample_qc <- build_sample_qc(metadata, counts_all, read_summaries)
  no_ligation_qc <- build_no_ligation_qc(sample_qc, params)
  offsets <- build_offset_matrices(counts_model, bedpe, read_summaries, contiguous_counts, params)

  structure(
    list(
      params = params,
      metadata = metadata,
      model_metadata = model_metadata,
      counts_all = counts_all,
      counts_model = counts_model,
      bedpe = bedpe,
      sample_qc = sample_qc,
      no_ligation_qc = no_ligation_qc,
      read_summaries = read_summaries,
      contiguous_counts = contiguous_counts,
      native_interactions = native_interactions,
      native_metrics = native_metrics,
      native_padj_by_role = native_padj_by_role,
      log_offset = offsets$log_offset,
      pair_background_available = offsets$pair_background_available,
      sample_exposure = offsets$sample_exposure,
      normalization_diagnostics = offsets$normalization_diagnostics,
      design_formula = design_info$design_formula,
      design = design_info$design,
      contrast_coef = design_info$contrast_coef,
      batch_columns_used = design_info$batch_columns_used,
      design_diagnostics = design_info$design_diagnostics,
      batch_balance = design_info$batch_balance,
      case_condition = case_condition,
      control_condition = control_condition,
      contrast_id = paste0(sanitize_id(case_condition), "_vs_", sanitize_id(control_condition)),
      method_results = list(),
      concordance = NULL,
      high_confidence_results = NULL,
      postprocess = list(),
      analysis_warnings = empty_analysis_warnings()
    ),
    class = "InteractionAnalysis"
  )
}

requested_methods <- function(analysis_method) {
  method <- tolower(analysis_method)
  if (method == "both") {
    return(c("edgeR", "DESeq2"))
  }
  if (method == "edger") {
    return("edgeR")
  }
  if (method == "deseq2") {
    return("DESeq2")
  }
  stop("analysis_method must be one of edgeR, DESeq2, or both.", call. = FALSE)
}

custom_count_filter <- function(counts, params) {
  keep_min_count <- rowSums(counts >= params$min_inter_contr) >=
    params$min_n_samples_inter_contr
  keep_larger_zero <- rowSums(counts > 0) >= params$min_n_samples_greater_zero
  keep_min_count & keep_larger_zero
}

offset_normalized_counts <- function(counts, log_offset) {
  offset_scale <- exp(log_offset)
  row_center <- exp(rowMeans(log(offset_scale)))
  sweep(counts / offset_scale, 1, row_center, `*`)
}

rounding_diagnostics <- function(raw_counts, rounded_counts) {
  data.frame(
    sample_id = colnames(raw_counts),
    total_raw_count = colSums(raw_counts),
    total_rounded_count = colSums(rounded_counts),
    total_absolute_rounding_delta = colSums(abs(rounded_counts - raw_counts)),
    max_absolute_rounding_delta = apply(abs(rounded_counts - raw_counts), 2, max),
    fractional_row_count = colSums(abs(raw_counts - round(raw_counts)) > 1e-8),
    stringsAsFactors = FALSE
  )
}

empty_analysis_warnings <- function() {
  data.frame(
    warning_index = integer(),
    stage = character(),
    method = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

capture_analysis_warnings <- function(expr, stage, method = NA_character_) {
  warning_rows <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
        stage = stage,
        method = method,
        message = conditionMessage(w),
        stringsAsFactors = FALSE
      )
      invokeRestart("muffleWarning")
    }
  )

  warnings <- if (length(warning_rows) > 0) {
    captured <- dplyr::bind_rows(warning_rows)
    captured$warning_index <- seq_len(nrow(captured))
    captured[, c("warning_index", "stage", "method", "message"), drop = FALSE]
  } else {
    empty_analysis_warnings()
  }

  list(value = value, warnings = warnings)
}

append_analysis_warnings <- function(existing, warnings) {
  if (is.null(existing)) {
    existing <- empty_analysis_warnings()
  }
  if (is.null(warnings) || nrow(warnings) == 0) {
    return(existing)
  }

  warnings$warning_index <- seq.int(nrow(existing) + 1L, length.out = nrow(warnings))
  warnings <- warnings[, c("warning_index", "stage", "method", "message"), drop = FALSE]
  dplyr::bind_rows(existing, warnings)
}

annotation_table <- function(analysis) {
  bedpe <- analysis$bedpe
  data.frame(
    interaction_id = bedpe$interaction_id,
    chr1 = bedpe$chr1,
    start1 = bedpe$start1,
    end1 = bedpe$end1,
    chr2 = bedpe$chr2,
    start2 = bedpe$start2,
    end2 = bedpe$end2,
    strand1 = bedpe$strand1,
    strand2 = bedpe$strand2,
    cluster_ids = bedpe$cluster_ids,
    arm1_features = bedpe$arm1_features,
    arm2_features = bedpe$arm2_features,
    stringsAsFactors = FALSE
  )
}

augment_method_results <- function(analysis, method_result) {
  results <- method_result$results
  anno <- annotation_table(analysis)
  metrics <- analysis$native_metrics

  augmented <- dplyr::left_join(results, anno, by = "interaction_id")
  augmented <- dplyr::left_join(augmented, metrics, by = "interaction_id")

  raw_counts <- analysis$counts_model[augmented$interaction_id, , drop = FALSE]
  names(raw_counts) <- paste0("raw_count_", colnames(raw_counts))
  normalized_counts <- method_result$normalized_counts[augmented$interaction_id, , drop = FALSE]
  names(normalized_counts) <- paste0("normalized_count_", colnames(normalized_counts))

  control_samples <- analysis$model_metadata$sample_id[analysis$model_metadata$role == analysis$params$control_role]
  case_samples <- analysis$model_metadata$sample_id[analysis$model_metadata$role == analysis$params$case_role]
  qc_samples <- analysis$metadata$sample_id[analysis$metadata$role %in% analysis$params$qc_roles]
  qc_samples <- intersect(qc_samples, colnames(analysis$counts_all))

  augmented$control_raw_mean <- rowMeans(analysis$counts_model[augmented$interaction_id, control_samples, drop = FALSE])
  augmented$case_raw_mean <- rowMeans(analysis$counts_model[augmented$interaction_id, case_samples, drop = FALSE])
  if (length(qc_samples) > 0) {
    qc_counts <- analysis$counts_all[augmented$interaction_id, qc_samples, drop = FALSE]
    augmented$no_ligation_max_count <- apply(qc_counts, 1, max, na.rm = TRUE)
    augmented$no_ligation_nonzero_samples <- rowSums(qc_counts > 0)
  } else {
    augmented$no_ligation_max_count <- NA_real_
    augmented$no_ligation_nonzero_samples <- NA_integer_
  }
  augmented$significant <- !is.na(augmented$padj) &
    augmented$padj <= analysis$params$padj_threshold

  dplyr::bind_cols(augmented, as.data.frame(raw_counts), as.data.frame(normalized_counts))
}

make_concordance <- function(method_results, padj_threshold, native_padj_by_role = NULL) {
  if (!all(c("edgeR", "DESeq2") %in% names(method_results))) {
    return(NULL)
  }
  edgeR_res <- method_results$edgeR$annotated_results[, c("interaction_id", "log2FoldChange", "pvalue", "padj")]
  names(edgeR_res) <- c("interaction_id", "edgeR_log2FoldChange", "edgeR_pvalue", "edgeR_padj")
  deseq_res <- method_results$DESeq2$annotated_results[, c("interaction_id", "log2FoldChange", "pvalue", "padj")]
  names(deseq_res) <- c("interaction_id", "DESeq2_log2FoldChange", "DESeq2_pvalue", "DESeq2_padj")

  concordance <- dplyr::full_join(edgeR_res, deseq_res, by = "interaction_id")
  if (!is.null(native_padj_by_role) && nrow(native_padj_by_role) > 0) {
    concordance <- dplyr::left_join(concordance, native_padj_by_role, by = "interaction_id")
  }
  concordance$edgeR_significant <- !is.na(concordance$edgeR_padj) &
    concordance$edgeR_padj <= padj_threshold
  concordance$DESeq2_significant <- !is.na(concordance$DESeq2_padj) &
    concordance$DESeq2_padj <= padj_threshold
  concordance$same_direction <- sign(concordance$edgeR_log2FoldChange) ==
    sign(concordance$DESeq2_log2FoldChange)
  concordance$concordance_class <- dplyr::case_when(
    concordance$edgeR_significant & concordance$DESeq2_significant & concordance$same_direction ~ "both_significant_same_direction",
    concordance$edgeR_significant & concordance$DESeq2_significant & !concordance$same_direction ~ "both_significant_opposite_direction",
    concordance$edgeR_significant & !concordance$DESeq2_significant ~ "edgeR_only",
    !concordance$edgeR_significant & concordance$DESeq2_significant ~ "DESeq2_only",
    TRUE ~ "not_significant"
  )
  concordance[order(concordance$edgeR_padj, concordance$DESeq2_padj), , drop = FALSE]
}

coverage_not_only_excluded <- function(profile_string, excluded_profiles) {
  if (is.na(profile_string) || !nzchar(profile_string)) {
    return(FALSE)
  }
  profiles <- unique(trimws(unlist(strsplit(profile_string, ";", fixed = TRUE), use.names = FALSE)))
  profiles <- profiles[nzchar(profiles)]
  length(profiles) > 0 && !all(profiles %in% excluded_profiles)
}

empty_high_confidence_results <- function() {
  data.frame(
    interaction_id = character(),
    high_confidence = logical(),
    stringsAsFactors = FALSE
  )
}

make_high_confidence_results <- function(analysis) {
  if (!all(c("edgeR", "DESeq2") %in% names(analysis$method_results))) {
    return(empty_high_confidence_results())
  }

  edgeR_results <- analysis$method_results$edgeR$annotated_results
  DESeq2_results <- analysis$method_results$DESeq2$annotated_results
  DESeq2_summary <- DESeq2_results[, intersect(
    c(
      "interaction_id",
      "log2FoldChange",
      "pvalue",
      "padj",
      "log2FoldChange_shrunken",
      "lfc_shrinkage_ratio"
    ),
    names(DESeq2_results)
  ), drop = FALSE]
  names(DESeq2_summary) <- sub("^log2FoldChange$", "DESeq2_log2FoldChange", names(DESeq2_summary))
  names(DESeq2_summary) <- sub("^pvalue$", "DESeq2_pvalue", names(DESeq2_summary))
  names(DESeq2_summary) <- sub("^padj$", "DESeq2_padj", names(DESeq2_summary))
  names(DESeq2_summary) <- sub("^log2FoldChange_shrunken$", "DESeq2_log2FoldChange_shrunken", names(DESeq2_summary))
  names(DESeq2_summary) <- sub("^lfc_shrinkage_ratio$", "DESeq2_lfc_shrinkage_ratio", names(DESeq2_summary))

  confidence <- dplyr::left_join(edgeR_results, DESeq2_summary, by = "interaction_id")
  confidence <- dplyr::left_join(
    confidence,
    analysis$normalization_diagnostics[, c("interaction_id", "pair_background_available_fraction"), drop = FALSE],
    by = "interaction_id"
  )

  confidence$high_conf_edgeR_significant <- confidence$significant &
    confidence$log2FoldChange > 0
  confidence$high_conf_DESeq2_significant_positive <- !is.na(confidence$DESeq2_padj) &
    confidence$DESeq2_padj <= analysis$params$padj_threshold &
    confidence$DESeq2_log2FoldChange > 0
  confidence$high_conf_same_positive_direction <- confidence$log2FoldChange > 0 &
    confidence$DESeq2_log2FoldChange > 0
  confidence$high_conf_rnanue_native <- !is.na(confidence$rnanue_min_padj_value) &
    confidence$rnanue_min_padj_value <= analysis$params$high_confidence_rnanue_padj_max
  confidence$high_conf_no_ligation_clean <- is.na(confidence$no_ligation_max_count) |
    confidence$no_ligation_max_count <= analysis$params$high_confidence_max_no_ligation_count
  confidence$high_conf_pair_background <- !is.na(confidence$pair_background_available_fraction) &
    confidence$pair_background_available_fraction >= analysis$params$high_confidence_min_pair_background_fraction
  confidence$high_conf_coverage_profile <- vapply(
    confidence$rnanue_coverage_profiles,
    coverage_not_only_excluded,
    logical(1),
    excluded_profiles = analysis$params$high_confidence_excluded_only_coverage_profiles
  )

  confidence$high_confidence <- confidence$high_conf_edgeR_significant &
    confidence$high_conf_DESeq2_significant_positive &
    confidence$high_conf_same_positive_direction &
    confidence$high_conf_rnanue_native &
    confidence$high_conf_no_ligation_clean &
    confidence$high_conf_pair_background &
    confidence$high_conf_coverage_profile

  confidence <- confidence[confidence$high_confidence, , drop = FALSE]
  confidence[order(confidence$padj, confidence$DESeq2_padj), , drop = FALSE]
}

high_confidence_bedpe <- function(high_confidence_results) {
  if (nrow(high_confidence_results) == 0) {
    return(data.frame(
      chr1 = character(),
      start1 = integer(),
      end1 = integer(),
      chr2 = character(),
      start2 = integer(),
      end2 = integer(),
      interaction_id = character(),
      score = numeric(),
      strand1 = character(),
      strand2 = character(),
      edgeR_log2FoldChange = numeric(),
      edgeR_padj = numeric(),
      DESeq2_log2FoldChange = numeric(),
      DESeq2_padj = numeric(),
      rnanue_min_padj_value = numeric(),
      no_ligation_max_count = numeric(),
      arm1_features = character(),
      arm2_features = character(),
      cluster_ids = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    chr1 = high_confidence_results$chr1,
    start1 = high_confidence_results$start1,
    end1 = high_confidence_results$end1,
    chr2 = high_confidence_results$chr2,
    start2 = high_confidence_results$start2,
    end2 = high_confidence_results$end2,
    interaction_id = high_confidence_results$interaction_id,
    score = -log10(pmax(high_confidence_results$padj, .Machine$double.xmin)),
    strand1 = high_confidence_results$strand1,
    strand2 = high_confidence_results$strand2,
    edgeR_log2FoldChange = high_confidence_results$log2FoldChange,
    edgeR_padj = high_confidence_results$padj,
    DESeq2_log2FoldChange = high_confidence_results$DESeq2_log2FoldChange,
    DESeq2_padj = high_confidence_results$DESeq2_padj,
    rnanue_min_padj_value = high_confidence_results$rnanue_min_padj_value,
    no_ligation_max_count = high_confidence_results$no_ligation_max_count,
    arm1_features = high_confidence_results$arm1_features,
    arm2_features = high_confidence_results$arm2_features,
    cluster_ids = high_confidence_results$cluster_ids,
    stringsAsFactors = FALSE
  )
}

high_confidence_region_records <- function(high_confidence_results) {
  empty_regions <- data.frame(
    chrom = character(),
    start = integer(),
    end = integer(),
    name = character(),
    score = integer(),
    strand = character(),
    interaction_id = character(),
    arm = character(),
    partner_chrom = character(),
    partner_start = integer(),
    partner_end = integer(),
    partner_strand = character(),
    edgeR_log2FoldChange = numeric(),
    edgeR_padj = numeric(),
    DESeq2_log2FoldChange = numeric(),
    DESeq2_padj = numeric(),
    rnanue_min_padj_value = numeric(),
    no_ligation_max_count = numeric(),
    arm_features = character(),
    partner_features = character(),
    cluster_ids = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(high_confidence_results) == 0) {
    return(empty_regions)
  }

  score <- pmin(
    1000L,
    as.integer(round(-log10(pmax(high_confidence_results$padj, .Machine$double.xmin)) * 100))
  )
  common <- data.frame(
    score = score,
    interaction_id = high_confidence_results$interaction_id,
    edgeR_log2FoldChange = high_confidence_results$log2FoldChange,
    edgeR_padj = high_confidence_results$padj,
    DESeq2_log2FoldChange = high_confidence_results$DESeq2_log2FoldChange,
    DESeq2_padj = high_confidence_results$DESeq2_padj,
    rnanue_min_padj_value = high_confidence_results$rnanue_min_padj_value,
    no_ligation_max_count = high_confidence_results$no_ligation_max_count,
    cluster_ids = high_confidence_results$cluster_ids,
    stringsAsFactors = FALSE
  )

  arm1 <- data.frame(
    chrom = high_confidence_results$chr1,
    start = high_confidence_results$start1,
    end = high_confidence_results$end1,
    name = paste0(high_confidence_results$interaction_id, "|arm1"),
    strand = high_confidence_results$strand1,
    arm = "arm1",
    partner_chrom = high_confidence_results$chr2,
    partner_start = high_confidence_results$start2,
    partner_end = high_confidence_results$end2,
    partner_strand = high_confidence_results$strand2,
    arm_features = high_confidence_results$arm1_features,
    partner_features = high_confidence_results$arm2_features,
    stringsAsFactors = FALSE
  )
  arm2 <- data.frame(
    chrom = high_confidence_results$chr2,
    start = high_confidence_results$start2,
    end = high_confidence_results$end2,
    name = paste0(high_confidence_results$interaction_id, "|arm2"),
    strand = high_confidence_results$strand2,
    arm = "arm2",
    partner_chrom = high_confidence_results$chr1,
    partner_start = high_confidence_results$start1,
    partner_end = high_confidence_results$end1,
    partner_strand = high_confidence_results$strand1,
    arm_features = high_confidence_results$arm2_features,
    partner_features = high_confidence_results$arm1_features,
    stringsAsFactors = FALSE
  )

  records <- dplyr::bind_rows(
    dplyr::bind_cols(arm1, common),
    dplyr::bind_cols(arm2, common)
  )
  records[, names(empty_regions), drop = FALSE]
}

high_confidence_regions_bed <- function(high_confidence_results) {
  records <- high_confidence_region_records(high_confidence_results)
  records[, c("chrom", "start", "end", "name", "score", "strand"), drop = FALSE]
}

write_interaction_outputs <- function(analysis) {
  dir.create(analysis$params$out_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(analysis$sample_qc, file.path(analysis$params$out_dir, "sample_qc.tsv"))
  readr::write_tsv(analysis$normalization_diagnostics, file.path(analysis$params$out_dir, "normalization_diagnostics.tsv"))
  readr::write_tsv(analysis$design_diagnostics, file.path(analysis$params$out_dir, "design_diagnostics.tsv"))
  readr::write_tsv(analysis$batch_balance, file.path(analysis$params$out_dir, "batch_balance.tsv"))
  readr::write_tsv(analysis$analysis_warnings, file.path(analysis$params$out_dir, "analysis_warnings.tsv"))
  if (nrow(analysis$no_ligation_qc) > 0) {
    readr::write_tsv(analysis$no_ligation_qc, file.path(analysis$params$out_dir, "no_ligation_qc.tsv"))
  }

  for (method in names(analysis$method_results)) {
    result <- analysis$method_results[[method]]
    readr::write_tsv(
      result$annotated_results,
      file.path(analysis$params$out_dir, paste0(method, "_results.tsv"))
    )

    sig <- result$annotated_results[result$annotated_results$significant, , drop = FALSE]
    readr::write_tsv(
      sig,
      file.path(analysis$params$out_dir, paste0("significant_interactions_", method, ".bedpe"))
    )

    if (!is.null(result$rounding)) {
      readr::write_tsv(
        result$rounding,
        file.path(analysis$params$out_dir, paste0(method, "_rounding_diagnostics.tsv"))
      )
    }
  }

  if (!is.null(analysis$concordance)) {
    readr::write_tsv(
      analysis$concordance,
      file.path(analysis$params$out_dir, "method_concordance.tsv")
    )
  }

  if (!is.null(analysis$high_confidence_results)) {
    readr::write_tsv(
      analysis$high_confidence_results,
      file.path(analysis$params$out_dir, "high_confidence_interactions.tsv")
    )
    readr::write_tsv(
      high_confidence_bedpe(analysis$high_confidence_results),
      file.path(analysis$params$out_dir, "high_confidence_interactions.bedpe")
    )
    readr::write_tsv(
      high_confidence_regions_bed(analysis$high_confidence_results),
      file.path(analysis$params$out_dir, "high_confidence_interaction_regions.bed"),
      col_names = FALSE
    )
    readr::write_tsv(
      high_confidence_region_records(analysis$high_confidence_results),
      file.path(analysis$params$out_dir, "high_confidence_interaction_regions.tsv")
    )
  }

  if (exists("write_postprocess_outputs", mode = "function")) {
    write_postprocess_outputs(analysis)
  }

  invisible(analysis)
}
