primary_result_method <- function(analysis) {
  if ("edgeR" %in% names(analysis$method_results)) {
    return("edgeR")
  }
  names(analysis$method_results)[1]
}

primary_normalized_counts <- function(analysis) {
  method <- primary_result_method(analysis)
  as.matrix(analysis$method_results[[method]]$normalized_counts)
}

interaction_result_classes <- function(analysis) {
  method <- primary_result_method(analysis)
  ids <- analysis$method_results[[method]]$annotated_results$interaction_id
  classes <- data.frame(
    interaction_id = ids,
    result_class = "not_significant",
    high_confidence = FALSE,
    stringsAsFactors = FALSE
  )

  if (!is.null(analysis$concordance)) {
    classes <- dplyr::left_join(
      classes,
      analysis$concordance[, c("interaction_id", "concordance_class"), drop = FALSE],
      by = "interaction_id"
    )
    classes$result_class <- ifelse(
      is.na(classes$concordance_class),
      classes$result_class,
      classes$concordance_class
    )
  } else {
    for (method_name in names(analysis$method_results)) {
      sig_ids <- analysis$method_results[[method_name]]$annotated_results$interaction_id[
        analysis$method_results[[method_name]]$annotated_results$significant
      ]
      classes$result_class[classes$interaction_id %in% sig_ids] <- paste0(method_name, "_significant")
    }
  }

  high_confidence_ids <- if (!is.null(analysis$high_confidence_results)) {
    analysis$high_confidence_results$interaction_id
  } else {
    character()
  }
  classes$high_confidence <- classes$interaction_id %in% high_confidence_ids
  classes$result_class[classes$high_confidence] <- "high_confidence"
  classes
}

pca_score_empty <- function() {
  data.frame(
    pca_set = character(),
    transform = character(),
    sample_id = character(),
    PC1 = numeric(),
    PC2 = numeric(),
    PC1_percent_variance = numeric(),
    PC2_percent_variance = numeric(),
    n_interactions = integer(),
    stringsAsFactors = FALSE
  )
}

pca_diag_row <- function(pca_set, transform, status, message, n_interactions) {
  data.frame(
    pca_set = pca_set,
    transform = transform,
    status = status,
    message = message,
    n_interactions = n_interactions,
    stringsAsFactors = FALSE
  )
}

visualization_batch_removed_matrix <- function(mat, analysis) {
  if (!analysis$params$batch_visualization_remove_batch ||
    length(analysis$batch_columns_used) == 0) {
    return(NULL)
  }
  keep_design <- stats::model.matrix(~ contrast_group, data = analysis$model_metadata)
  batch_cols <- setdiff(colnames(analysis$design), colnames(keep_design))
  if (length(batch_cols) == 0) {
    return(NULL)
  }
  fit <- limma::lmFit(mat, analysis$design)
  batch_coef <- fit$coefficients[, batch_cols, drop = FALSE]
  batch_coef[!is.finite(batch_coef)] <- 0
  batch_effect <- batch_coef %*% t(analysis$design[, batch_cols, drop = FALSE])
  mat - batch_effect
}

compute_pca_scores <- function(mat, analysis, pca_set, transform) {
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(list(
      scores = pca_score_empty(),
      diagnostics = pca_diag_row(
        pca_set,
        transform,
        "skipped",
        "Need at least two interactions and two samples for PCA.",
        nrow(mat)
      )
    ))
  }
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  if (ncol(pca$x) < 2 || sum(pca$sdev^2) == 0) {
    return(list(
      scores = pca_score_empty(),
      diagnostics = pca_diag_row(
        pca_set,
        transform,
        "skipped",
        "PCA returned fewer than two usable components.",
        nrow(mat)
      )
    ))
  }

  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(
    pca_set = pca_set,
    transform = transform,
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_percent_variance = var_exp[1],
    PC2_percent_variance = var_exp[2],
    n_interactions = nrow(mat),
    stringsAsFactors = FALSE
  )
  scores <- dplyr::left_join(scores, analysis$model_metadata, by = "sample_id")
  sample_order <- ordered_samples_by_role(analysis$model_metadata, scores$sample_id)
  scores$sample_id <- factor(scores$sample_id, levels = sample_order)
  scores$role <- factor(scores$role, levels = role_plot_levels())
  scores <- scores[order(scores$sample_id), , drop = FALSE]
  scores$sample_id <- as.character(scores$sample_id)
  list(
    scores = scores,
    diagnostics = pca_diag_row(pca_set, transform, "computed", "PCA computed.", nrow(mat))
  )
}

make_pca_assessments <- function(analysis, result_classes) {
  normalized <- primary_normalized_counts(analysis)
  normalized <- normalized[rownames(normalized) %in% result_classes$interaction_id, , drop = FALSE]

  sets <- list(
    all_model_filtered = rownames(normalized)
  )
  if ("edgeR" %in% names(analysis$method_results)) {
    edgeR_res <- analysis$method_results$edgeR$annotated_results
    sets$edgeR_significant <- edgeR_res$interaction_id[edgeR_res$significant]
  }
  if ("DESeq2" %in% names(analysis$method_results)) {
    deseq_res <- analysis$method_results$DESeq2$annotated_results
    sets$DESeq2_significant <- deseq_res$interaction_id[deseq_res$significant]
  }
  if (!is.null(analysis$concordance)) {
    sets$concordant_same_direction <- analysis$concordance$interaction_id[
      analysis$concordance$concordance_class == "both_significant_same_direction"
    ]
  }
  if (!is.null(analysis$high_confidence_results)) {
    sets$high_confidence <- analysis$high_confidence_results$interaction_id
  }

  scores <- list()
  diagnostics <- list()
  for (set_name in names(sets)) {
    ids <- intersect(sets[[set_name]], rownames(normalized))
    mat <- log2(normalized[ids, , drop = FALSE] + 1)
    raw <- compute_pca_scores(mat, analysis, set_name, "log2_offset_normalized")
    scores[[set_name]] <- raw$scores
    diagnostics[[paste(set_name, "raw", sep = "_")]] <- raw$diagnostics

    batch_removed <- visualization_batch_removed_matrix(mat, analysis)
    if (!is.null(batch_removed)) {
      adjusted <- compute_pca_scores(batch_removed, analysis, set_name, "batch_removed_visualization")
      scores[[set_name]] <- dplyr::bind_rows(scores[[set_name]], adjusted$scores)
      diagnostics[[paste(set_name, "batch_removed", sep = "_")]] <- adjusted$diagnostics
    }
  }

  list(
    scores = scores,
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}

safe_cv <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  if (!is.finite(mean_x) || mean_x <= 0) {
    return(NA_real_)
  }
  stats::sd(x, na.rm = TRUE) / mean_x
}

make_candidate_stability <- function(analysis, result_classes) {
  method <- primary_result_method(analysis)
  primary_results <- analysis$method_results[[method]]$annotated_results
  normalized <- primary_normalized_counts(analysis)
  normalized <- normalized[primary_results$interaction_id, , drop = FALSE]
  raw_counts <- analysis$counts_model[primary_results$interaction_id, , drop = FALSE]

  case_samples <- analysis$model_metadata$sample_id[analysis$model_metadata$role == analysis$params$case_role]
  control_samples <- analysis$model_metadata$sample_id[analysis$model_metadata$role == analysis$params$control_role]

  stability <- data.frame(
    interaction_id = primary_results$interaction_id,
    primary_method = method,
    treatment_norm_mean = rowMeans(normalized[, case_samples, drop = FALSE], na.rm = TRUE),
    control_norm_mean = rowMeans(normalized[, control_samples, drop = FALSE], na.rm = TRUE),
    treatment_norm_cv = apply(normalized[, case_samples, drop = FALSE], 1, safe_cv),
    control_norm_cv = apply(normalized[, control_samples, drop = FALSE], 1, safe_cv),
    treatment_detection_prevalence = rowMeans(raw_counts[, case_samples, drop = FALSE] > 0),
    control_detection_prevalence = rowMeans(raw_counts[, control_samples, drop = FALSE] > 0),
    single_sample_dominance = apply(
      normalized,
      1,
      function(x) if (sum(x, na.rm = TRUE) > 0) max(x, na.rm = TRUE) / sum(x, na.rm = TRUE) else NA_real_
    ),
    min_treatment_over_max_control = (
      apply(normalized[, case_samples, drop = FALSE], 1, min, na.rm = TRUE) + 0.5
    ) / (
      apply(normalized[, control_samples, drop = FALSE], 1, max, na.rm = TRUE) + 0.5
    ),
    stringsAsFactors = FALSE
  )

  stability <- dplyr::left_join(
    stability,
    primary_results[, intersect(
      c(
        "interaction_id", "log2FoldChange", "pvalue", "padj", "significant",
        "no_ligation_max_count", "no_ligation_nonzero_samples",
        "rnanue_min_padj_value", "rnanue_mean_support_per_effective_bp",
        "rnanue_median_arm_balance", "rnanue_max_coverage_components",
        "rnanue_coverage_profiles"
      ),
      names(primary_results)
    ), drop = FALSE],
    by = "interaction_id"
  )
  dplyr::left_join(stability, result_classes, by = "interaction_id")
}

safe_median <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

make_rnanue_metric_group_summary <- function(candidate_stability) {
  groups <- split(candidate_stability, candidate_stability$result_class)
  rows <- lapply(names(groups), function(group_name) {
    group <- groups[[group_name]]
    data.frame(
      result_class = group_name,
      interactions = nrow(group),
      median_rnanue_min_padj_value = safe_median(group$rnanue_min_padj_value),
      median_support_per_effective_bp = safe_median(group$rnanue_mean_support_per_effective_bp),
      median_arm_balance = safe_median(group$rnanue_median_arm_balance),
      median_max_coverage_components = safe_median(group$rnanue_max_coverage_components),
      median_no_ligation_max_count = safe_median(group$no_ligation_max_count),
      median_single_sample_dominance = safe_median(group$single_sample_dominance),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

make_partner_recurrence <- function(analysis, result_classes) {
  method <- primary_result_method(analysis)
  results <- analysis$method_results[[method]]$annotated_results
  arm_feature_rows <- function(interaction_id, arm, features) {
    if (length(features) == 0) {
      return(data.frame(
        interaction_id = character(),
        arm = character(),
        feature_id = character(),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      interaction_id = interaction_id,
      arm = arm,
      feature_id = features,
      stringsAsFactors = FALSE
    )
  }
  rows <- lapply(seq_len(nrow(results)), function(i) {
    interaction_id <- results$interaction_id[i]
    dplyr::bind_rows(
      arm_feature_rows(
        interaction_id = interaction_id,
        arm = "arm1",
        features = parse_feature_field(results$arm1_features[i])
      ),
      arm_feature_rows(
        interaction_id = interaction_id,
        arm = "arm2",
        features = parse_feature_field(results$arm2_features[i])
      )
    )
  })
  feature_rows <- dplyr::bind_rows(rows)
  if (nrow(feature_rows) == 0) {
    return(data.frame(
      feature_id = character(),
      total_interactions = integer(),
      significant_interactions = integer(),
      high_confidence_interactions = integer(),
      concordant_interactions = integer(),
      arm1_interactions = integer(),
      arm2_interactions = integer(),
      stringsAsFactors = FALSE
    ))
  }
  feature_rows <- dplyr::left_join(feature_rows, result_classes, by = "interaction_id")
  split_features <- split(feature_rows, feature_rows$feature_id)
  recurrence <- lapply(names(split_features), function(feature_id) {
    feature <- split_features[[feature_id]]
    data.frame(
      feature_id = feature_id,
      total_interactions = length(unique(feature$interaction_id)),
      significant_interactions = length(unique(feature$interaction_id[
        feature$result_class != "not_significant"
      ])),
      high_confidence_interactions = length(unique(feature$interaction_id[feature$high_confidence])),
      concordant_interactions = length(unique(feature$interaction_id[
        feature$result_class %in% c("both_significant_same_direction", "high_confidence")
      ])),
      arm1_interactions = length(unique(feature$interaction_id[feature$arm == "arm1"])),
      arm2_interactions = length(unique(feature$interaction_id[feature$arm == "arm2"])),
      stringsAsFactors = FALSE
    )
  })
  recurrence <- dplyr::bind_rows(recurrence)
  recurrence[order(-recurrence$total_interactions, -recurrence$high_confidence_interactions), , drop = FALSE]
}

called_interaction_ids <- function(analysis) {
  ids <- character()
  for (method in names(analysis$method_results)) {
    results <- analysis$method_results[[method]]$annotated_results
    ids <- c(ids, results$interaction_id[results$significant])
  }
  if (!is.null(analysis$high_confidence_results)) {
    ids <- c(ids, analysis$high_confidence_results$interaction_id)
  }
  unique(ids)
}

make_significant_normalized_counts <- function(analysis, result_classes) {
  ids <- called_interaction_ids(analysis)
  normalized <- primary_normalized_counts(analysis)
  ids <- intersect(ids, rownames(normalized))
  if (length(ids) == 0) {
    return(data.frame(
      interaction_id = character(),
      sample_id = character(),
      normalized_count = numeric(),
      log2_normalized_count = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  long_counts <- as.data.frame(as.table(normalized[ids, , drop = FALSE]), stringsAsFactors = FALSE)
  names(long_counts) <- c("interaction_id", "sample_id", "normalized_count")
  long_counts$log2_normalized_count <- log2(long_counts$normalized_count + 1)
  long_counts <- dplyr::left_join(long_counts, analysis$model_metadata, by = "sample_id")
  long_counts <- dplyr::left_join(long_counts, result_classes, by = "interaction_id")

  method <- primary_result_method(analysis)
  primary_results <- analysis$method_results[[method]]$annotated_results
  stat_cols <- intersect(
    c("interaction_id", "log2FoldChange", "pvalue", "padj", "significant", "no_ligation_max_count"),
    names(primary_results)
  )
  long_counts <- dplyr::left_join(
    long_counts,
    primary_results[, stat_cols, drop = FALSE],
    by = "interaction_id"
  )
  sample_order <- ordered_samples_by_role(analysis$model_metadata, long_counts$sample_id)
  long_counts$sample_id <- factor(long_counts$sample_id, levels = sample_order)
  long_counts <- long_counts[order(long_counts$padj, long_counts$interaction_id, long_counts$sample_id), , drop = FALSE]
  long_counts$sample_id <- as.character(long_counts$sample_id)
  long_counts
}

empty_dia_rnanue_padj_comparison <- function() {
  data.frame(
    interaction_id = character(),
    primary_method = character(),
    primary_dia_log2FoldChange = numeric(),
    primary_dia_pvalue = numeric(),
    primary_dia_padj = numeric(),
    primary_dia_padj_floor = numeric(),
    primary_dia_neg_log10_padj = numeric(),
    significant = logical(),
    result_class = character(),
    high_confidence = logical(),
    rnanue_condition_sample_padj_values = character(),
    rnanue_condition_padj_n = integer(),
    rnanue_condition_padj_zero_count = integer(),
    rnanue_condition_padj_min = numeric(),
    rnanue_condition_padj_max = numeric(),
    rnanue_condition_padj_arithmetic_mean = numeric(),
    rnanue_condition_padj_geomean = numeric(),
    rnanue_condition_padj_floor = numeric(),
    rnanue_condition_neg_log10_padj_mean = numeric(),
    rnanue_condition_neg_log10_padj_sd = numeric(),
    stringsAsFactors = FALSE
  )
}

padj_plot_floor <- function(values) {
  positive <- values[is.finite(values) & values > 0]
  if (length(positive) == 0) {
    return(.Machine$double.xmin)
  }
  max(min(positive, na.rm = TRUE) / 10, .Machine$double.xmin)
}

format_named_padj_values <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA_character_)
  }
  formatted <- trimws(formatC(values, digits = 6, format = "g"))
  paste(paste(names(values), formatted, sep = ":"), collapse = ";")
}

condition_sample_rnanue_padj_values <- function(analysis, interaction_id, condition_samples) {
  out <- rep(NA_real_, length(condition_samples))
  names(out) <- condition_samples
  bedpe_i <- match(interaction_id, analysis$bedpe$interaction_id)
  if (is.na(bedpe_i)) {
    return(out)
  }
  clusters <- parse_cluster_tokens(analysis$bedpe$cluster_ids[bedpe_i])
  if (nrow(clusters) == 0) {
    return(out)
  }

  for (sample_id in condition_samples) {
    sample_clusters <- clusters$cluster_ID[clusters$sample_id == sample_id]
    sample_df <- analysis$native_interactions[[sample_id]]
    if (length(sample_clusters) == 0 ||
      is.null(sample_df) ||
      !"padj_value" %in% names(sample_df)) {
      next
    }
    sample_clusters <- sample_clusters[sample_clusters %in% rownames(sample_df)]
    if (length(sample_clusters) == 0) {
      next
    }
    values <- suppressWarnings(as.numeric(sample_df[sample_clusters, "padj_value", drop = TRUE]))
    values <- values[is.finite(values) & values >= 0]
    if (length(values) > 0) {
      out[sample_id] <- min(values, na.rm = TRUE)
    }
  }
  out
}

geometric_mean_padj <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  if (any(values == 0)) {
    return(0)
  }
  exp(mean(log(values)))
}

make_dia_rnanue_padj_comparison <- function(analysis, result_classes) {
  method <- primary_result_method(analysis)
  primary_results <- analysis$method_results[[method]]$annotated_results
  primary_results <- primary_results[is.finite(primary_results$padj), , drop = FALSE]
  if (nrow(primary_results) == 0) {
    return(empty_dia_rnanue_padj_comparison())
  }

  condition_samples <- analysis$model_metadata$sample_id[
    analysis$model_metadata$role == analysis$params$case_role
  ]
  sample_values <- lapply(
    primary_results$interaction_id,
    condition_sample_rnanue_padj_values,
    analysis = analysis,
    condition_samples = condition_samples
  )
  rnanue_floor <- padj_plot_floor(unlist(sample_values, use.names = FALSE))
  dia_floor <- padj_plot_floor(primary_results$padj)

  rows <- lapply(seq_len(nrow(primary_results)), function(i) {
    values <- sample_values[[i]]
    present_values <- values[!is.na(values)]
    neg_log_values <- -log10(pmax(present_values, rnanue_floor))
    data.frame(
      interaction_id = primary_results$interaction_id[i],
      primary_method = method,
      primary_dia_log2FoldChange = primary_results$log2FoldChange[i],
      primary_dia_pvalue = primary_results$pvalue[i],
      primary_dia_padj = primary_results$padj[i],
      primary_dia_padj_floor = dia_floor,
      primary_dia_neg_log10_padj = -log10(pmax(primary_results$padj[i], dia_floor)),
      significant = primary_results$significant[i],
      rnanue_condition_sample_padj_values = format_named_padj_values(values),
      rnanue_condition_padj_n = length(present_values),
      rnanue_condition_padj_zero_count = sum(present_values == 0),
      rnanue_condition_padj_min = if (length(present_values) == 0) NA_real_ else min(present_values),
      rnanue_condition_padj_max = if (length(present_values) == 0) NA_real_ else max(present_values),
      rnanue_condition_padj_arithmetic_mean = if (length(present_values) == 0) NA_real_ else mean(present_values),
      rnanue_condition_padj_geomean = geometric_mean_padj(present_values),
      rnanue_condition_padj_floor = rnanue_floor,
      rnanue_condition_neg_log10_padj_mean = if (length(neg_log_values) == 0) NA_real_ else mean(neg_log_values),
      rnanue_condition_neg_log10_padj_sd = if (length(neg_log_values) < 2) NA_real_ else stats::sd(neg_log_values),
      stringsAsFactors = FALSE
    )
  })

  comparison <- dplyr::bind_rows(rows)
  comparison <- dplyr::left_join(comparison, result_classes, by = "interaction_id")
  comparison <- comparison[order(comparison$primary_dia_padj, comparison$interaction_id, na.last = TRUE), , drop = FALSE]
  comparison[, names(empty_dia_rnanue_padj_comparison()), drop = FALSE]
}

empty_dia_rnanue_discordance <- function() {
  data.frame(
    interaction_id = character(),
    evidence_class = character(),
    primary_method = character(),
    primary_dia_log2FoldChange = numeric(),
    primary_dia_pvalue = numeric(),
    primary_dia_padj = numeric(),
    primary_dia_neg_log10_padj = numeric(),
    dia_positive = logical(),
    rnanue_positive = logical(),
    rnanue_support_cutoff = numeric(),
    rnanue_condition_padj_n = integer(),
    rnanue_condition_padj_geomean = numeric(),
    rnanue_condition_neg_log10_padj_mean = numeric(),
    rnanue_condition_neg_log10_padj_sd = numeric(),
    treatment_norm_mean = numeric(),
    control_norm_mean = numeric(),
    treatment_norm_cv = numeric(),
    control_norm_cv = numeric(),
    single_sample_dominance = numeric(),
    min_treatment_over_max_control = numeric(),
    no_ligation_max_count = numeric(),
    no_ligation_nonzero_samples = numeric(),
    result_class = character(),
    high_confidence = logical(),
    stringsAsFactors = FALSE
  )
}

make_dia_rnanue_discordance <- function(comparison, candidate_stability, params) {
  if (is.null(comparison) || nrow(comparison) == 0) {
    return(empty_dia_rnanue_discordance())
  }

  rnanue_cutoff <- params$high_confidence_rnanue_padj_max %||% 0.1
  discordance <- comparison
  discordance$dia_positive <- !is.na(discordance$significant) & discordance$significant
  discordance$rnanue_positive <- is.finite(discordance$rnanue_condition_padj_geomean) &
    discordance$rnanue_condition_padj_geomean <= rnanue_cutoff
  discordance$rnanue_support_cutoff <- rnanue_cutoff
  discordance$evidence_class <- dplyr::case_when(
    discordance$dia_positive & discordance$rnanue_positive ~ "DIA+/RNAnue+",
    discordance$dia_positive & !discordance$rnanue_positive ~ "DIA-only",
    !discordance$dia_positive & discordance$rnanue_positive ~ "RNAnue-only",
    TRUE ~ "neither"
  )

  stability_cols <- intersect(
    c(
      "interaction_id", "treatment_norm_mean", "control_norm_mean",
      "treatment_norm_cv", "control_norm_cv", "single_sample_dominance",
      "min_treatment_over_max_control", "no_ligation_max_count",
      "no_ligation_nonzero_samples"
    ),
    names(candidate_stability)
  )
  if (length(stability_cols) > 1) {
    discordance <- dplyr::left_join(
      discordance,
      candidate_stability[, stability_cols, drop = FALSE],
      by = "interaction_id"
    )
  }

  for (col in setdiff(names(empty_dia_rnanue_discordance()), names(discordance))) {
    discordance[[col]] <- NA
  }
  discordance <- discordance[order(discordance$evidence_class, discordance$primary_dia_padj), , drop = FALSE]
  discordance[, names(empty_dia_rnanue_discordance()), drop = FALSE]
}

make_postprocess_assessments <- function(analysis) {
  result_classes <- interaction_result_classes(analysis)
  pca <- make_pca_assessments(analysis, result_classes)
  candidate_stability <- make_candidate_stability(analysis, result_classes)
  gff_annotation <- make_gff_annotation_assessments(analysis, result_classes)
  dia_rnanue_padj_comparison <- make_dia_rnanue_padj_comparison(analysis, result_classes)

  list(
    result_classes = result_classes,
    pca_scores = pca$scores,
    pca_diagnostics = pca$diagnostics,
    candidate_stability = candidate_stability,
    rnanue_metric_group_summary = make_rnanue_metric_group_summary(candidate_stability),
    partner_recurrence = make_partner_recurrence(analysis, result_classes),
    significant_normalized_counts = make_significant_normalized_counts(analysis, result_classes),
    dia_rnanue_padj_comparison = dia_rnanue_padj_comparison,
    dia_rnanue_discordance = make_dia_rnanue_discordance(
      dia_rnanue_padj_comparison,
      candidate_stability,
      analysis$params
    ),
    gff_annotation = gff_annotation
  )
}

write_postprocess_outputs <- function(analysis) {
  if (length(analysis$postprocess) == 0) {
    return(invisible(NULL))
  }
  unlink(file.path(
    analysis$params$out_dir,
    c(
      "requested_feature_dia_results.tsv",
      "requested_feature_normalized_counts.tsv",
      "dia_rnanue_correlation_stats.tsv",
      "dia_rnanue_threshold_enrichment.tsv"
    )
  ))

  for (set_name in names(analysis$postprocess$pca_scores)) {
    readr::write_tsv(
      analysis$postprocess$pca_scores[[set_name]],
      file.path(analysis$params$out_dir, paste0("pca_scores_", sanitize_id(set_name), ".tsv"))
    )
  }
  readr::write_tsv(
    analysis$postprocess$pca_diagnostics,
    file.path(analysis$params$out_dir, "pca_diagnostics.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$candidate_stability,
    file.path(analysis$params$out_dir, "candidate_stability.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$rnanue_metric_group_summary,
    file.path(analysis$params$out_dir, "rnanue_metric_group_summary.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$partner_recurrence,
    file.path(analysis$params$out_dir, "partner_recurrence.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$significant_normalized_counts,
    file.path(analysis$params$out_dir, "significant_normalized_counts.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$dia_rnanue_padj_comparison,
    file.path(analysis$params$out_dir, "dia_rnanue_padj_comparison.tsv")
  )
  readr::write_tsv(
    analysis$postprocess$dia_rnanue_discordance,
    file.path(analysis$params$out_dir, "dia_rnanue_discordance.tsv")
  )
  if (!is.null(analysis$postprocess$gff_annotation)) {
    readr::write_tsv(
      analysis$postprocess$gff_annotation$arm_annotations,
      file.path(analysis$params$out_dir, "interaction_arm_gff_annotations.tsv")
    )
    readr::write_tsv(
      analysis$postprocess$gff_annotation$partner_type_pairs,
      file.path(analysis$params$out_dir, "interaction_partner_type_pairs.tsv")
    )
    readr::write_tsv(
      analysis$postprocess$gff_annotation$partner_type_confusion,
      file.path(analysis$params$out_dir, "partner_type_confusion.tsv")
    )
    readr::write_tsv(
      analysis$postprocess$gff_annotation$target_gene_interactions,
      file.path(analysis$params$out_dir, "target_gene_interactions.tsv")
    )
    readr::write_tsv(
      analysis$postprocess$gff_annotation$target_gene_bedpe,
      file.path(analysis$params$out_dir, "target_gene_interactions.bedpe")
    )
    readr::write_tsv(
      analysis$postprocess$gff_annotation$diagnostics,
      file.path(analysis$params$out_dir, "gff_annotation_diagnostics.tsv")
    )
  }

  invisible(NULL)
}
