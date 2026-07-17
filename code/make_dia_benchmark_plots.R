#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(jsonlite)
  library(readr)
  library(scales)
  library(tidyr)
})

call_set_labels <- c(
  RNAnue_native_padj = "RNAnue native",
  edgeR = "edgeR",
  DESeq2 = "DESeq2",
  high_confidence_DIA = "High-confidence DIA"
)

call_set_levels <- unname(call_set_labels[c(
  "DESeq2",
  "edgeR",
  "RNAnue_native_padj",
  "high_confidence_DIA"
)])

method_colors <- c(
  "RNAnue native" = "#4C78A8",
  "edgeR" = "#238B35",
  "DESeq2" = "#CC6677",
  "High-confidence DIA" = "#B12A7A"
)

truth_colors <- c(
  "True Chimera positive" = "#238B35",
  "Background candidate" = "#9B9B9B"
)

metric_labels <- c(
  precision = "Precision",
  recall = "Recall",
  f1 = "F1",
  auprc = "AUPRC"
)

parse_args <- function(args) {
  get_value <- function(flag, default = NULL) {
    index <- match(flag, args)
    if (is.na(index)) {
      return(default)
    }
    if (index == length(args)) {
      stop("Missing value after ", flag, call. = FALSE)
    }
    args[[index + 1]]
  }

  benchmark_dir <- get_value("--benchmark-dir")
  if (is.null(benchmark_dir)) {
    stop("--benchmark-dir is required", call. = FALSE)
  }

  list(
    benchmark_dir = normalizePath(benchmark_dir, mustWork = TRUE),
    output_dir = get_value(
      "--output-dir",
      file.path("results", "benchmarks", paste0("dia_pipeline_v", Sys.Date()))
    ),
    formats = strsplit(get_value("--formats", "png,pdf,svg"), ",", fixed = TRUE)[[1]],
    dpi = as.integer(get_value("--dpi", "300")),
    padj_threshold = as.numeric(get_value("--padj-threshold", "0.1"))
  )
}

find_one <- function(directory, pattern) {
  matches <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(matches) != 1) {
    stop(
      "Expected exactly one file matching ",
      pattern,
      " in ",
      directory,
      "; found ",
      length(matches),
      call. = FALSE
    )
  }
  matches[[1]]
}

label_call_set <- function(x) {
  out <- as.character(x)
  matched <- out %in% names(call_set_labels)
  out[matched] <- unname(call_set_labels[out[matched]])
  out
}

condition_display <- function(condition_id, treatment_fraction) {
  ifelse(
    is.na(treatment_fraction) | treatment_fraction == 0,
    "Null",
    paste0(scales::percent(treatment_fraction, accuracy = 1), " True Chimera")
  )
}

format_cutoff <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return("NA")
  }
  format(x, scientific = FALSE, trim = TRUE)
}

benchmark_subtitle <- function(metrics, padj_threshold = NULL) {
  non_null <- metrics %>%
    filter(treatment_target_true_chimera_fraction > 0)
  treatment_fractions <- sort(unique(non_null$treatment_target_true_chimera_fraction))

  parts <- character()
  if (!is.null(padj_threshold)) {
    parts <- c(parts, paste0("Cutoff ", format_cutoff(padj_threshold)))
  } else if ("padj_threshold" %in% names(metrics)) {
    parts <- c(parts, paste0("Cutoff ", format_cutoff(unique(metrics$padj_threshold)[[1]])))
  }
  if ("total_reads" %in% names(metrics)) {
    parts <- c(parts, paste0(scales::comma(unique(metrics$total_reads)[[1]]), " reads"))
  }
  if ("read_length" %in% names(metrics)) {
    parts <- c(parts, paste0(unique(metrics$read_length)[[1]], " nt"))
  }
  if ("background_fraction" %in% names(metrics)) {
    parts <- c(parts, paste0("background ", scales::percent(unique(metrics$background_fraction)[[1]], accuracy = 1)))
  }
  if (length(treatment_fractions) > 0) {
    parts <- c(
      parts,
      paste0(
        "treatment True Chimera ",
        paste(scales::percent(treatment_fractions, accuracy = 1), collapse = "/")
      )
    )
  }

  paste(parts, collapse = " | ")
}

benchmark_theme <- function(base_size = 11) {
  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
      strip.text = ggplot2::element_text(face = "bold", color = "grey15")
    )
}

save_plot <- function(plot, output_dir, stem, width = 10, height = 6, formats = c("png", "pdf"), dpi = 300) {
  paths <- character()
  for (format in formats) {
    format <- trimws(tolower(format))
    if (!nzchar(format)) {
      next
    }
    path <- file.path(output_dir, paste0(stem, ".", format))
    ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
    paths <- c(paths, path)
  }
  paths
}

read_benchmark_tables <- function(benchmark_dir) {
  metrics_path <- find_one(benchmark_dir, "\\.metrics\\.tsv$")
  candidate_path <- find_one(benchmark_dir, "\\.candidate_labels\\.tsv$")
  summary_path <- find_one(benchmark_dir, "\\.summary\\.json$")
  config_path <- find_one(benchmark_dir, "\\.config\\.json$")

  metrics <- readr::read_tsv(metrics_path, na = c("", "NA", "NaN"), show_col_types = FALSE)
  candidates <- readr::read_tsv(candidate_path, na = c("", "NA", "NaN"), show_col_types = FALSE)
  summary <- jsonlite::fromJSON(summary_path, flatten = TRUE)
  config <- jsonlite::fromJSON(config_path, flatten = TRUE)

  metrics <- metrics %>%
    mutate(
      call_set_label = factor(label_call_set(call_set), levels = call_set_levels),
      condition_label = condition_display(condition_id, treatment_target_true_chimera_fraction),
      condition_label = factor(condition_label, levels = unique(condition_label))
    )

  candidates <- candidates %>%
    mutate(truth_positive = as.logical(truth_positive))

  list(
    metrics = metrics,
    candidates = candidates,
    summary = summary,
    config = config,
    paths = list(
      metrics = metrics_path,
      candidates = candidate_path,
      summary = summary_path,
      config = config_path
    )
  )
}

read_method_scores <- function(benchmark_dir, candidates) {
  condition_root <- file.path(benchmark_dir, "conditions")
  condition_dirs <- list.dirs(condition_root, recursive = FALSE, full.names = TRUE)
  rows <- list()
  index <- 1

  for (condition_dir in condition_dirs) {
    condition_id <- basename(condition_dir)
    condition_candidates <- candidates %>%
      filter(.data$condition_id == .env$condition_id) %>%
      select(condition_id, interaction_id, truth_positive)

    for (method in c("edgeR", "DESeq2")) {
      result_path <- file.path(condition_dir, "dia", "results", paste0(method, "_results.tsv"))
      if (!file.exists(result_path)) {
        next
      }
      result <- readr::read_tsv(result_path, na = c("", "NA", "NaN"), show_col_types = FALSE) %>%
        select(interaction_id, log2FoldChange, pvalue, padj) %>%
        mutate(
          condition_id = .env$condition_id,
          method = .env$method,
          padj = as.numeric(padj),
          log2FoldChange = as.numeric(log2FoldChange),
          score_padj = padj
        ) %>%
        left_join(condition_candidates, by = c("condition_id", "interaction_id"))
      rows[[index]] <- result
      index <- index + 1
    }

    high_confidence_path <- file.path(condition_dir, "dia", "results", "high_confidence_interactions.tsv")
    if (file.exists(high_confidence_path)) {
      high_confidence <- readr::read_tsv(high_confidence_path, na = c("", "NA", "NaN"), show_col_types = FALSE) %>%
        filter(is.na(high_confidence) | high_confidence %in% c(TRUE, "TRUE", "true", "1", 1)) %>%
        mutate(
          edgeR_padj_value = as.numeric(padj),
          DESeq2_padj_value = as.numeric(DESeq2_padj),
          edgeR_log2FoldChange_value = as.numeric(log2FoldChange),
          DESeq2_log2FoldChange_value = as.numeric(DESeq2_log2FoldChange)
        ) %>%
        transmute(
          interaction_id,
          log2FoldChange = pmin(edgeR_log2FoldChange_value, DESeq2_log2FoldChange_value),
          pvalue = NA_real_,
          padj = pmax(edgeR_padj_value, DESeq2_padj_value),
          condition_id = .env$condition_id,
          method = "high_confidence_DIA",
          score_padj = edgeR_padj_value
        ) %>%
        mutate(
          log2FoldChange = ifelse(is.infinite(log2FoldChange), NA_real_, log2FoldChange),
          padj = ifelse(is.infinite(padj), NA_real_, padj),
          score_padj = ifelse(is.infinite(score_padj), NA_real_, score_padj)
        ) %>%
        left_join(condition_candidates, by = c("condition_id", "interaction_id"))
      rows[[index]] <- high_confidence
      index <- index + 1
    }
  }

  rnanue_scores <- candidates %>%
    transmute(
      condition_id,
      interaction_id,
      log2FoldChange = NA_real_,
      pvalue = NA_real_,
      padj = as.numeric(native_treatment_min_padj),
      method = "RNAnue_native_padj",
      score_padj = padj,
      truth_positive
    )

  bind_rows(rows, list(rnanue_scores)) %>%
    mutate(
      method_label = label_call_set(method),
      method_label = factor(method_label, levels = call_set_levels)
    )
}

score_from_padj <- function(padj) {
  padj <- as.numeric(padj)
  score <- -log10(pmax(padj, .Machine$double.xmin))
  score[is.na(padj) | !is.finite(padj)] <- NA_real_
  score
}

average_precision <- function(labels, scores) {
  positive_count <- sum(labels, na.rm = TRUE)
  if (positive_count == 0) {
    return(NA_real_)
  }

  score_values <- scores[names(labels)]
  score_values[is.na(score_values)] <- -Inf
  ranked <- tibble::tibble(
    interaction_id = names(labels),
    truth_positive = as.logical(labels),
    score = as.numeric(score_values)
  ) %>%
    arrange(desc(score), desc(interaction_id))

  true_indices <- which(ranked$truth_positive)
  sum(seq_along(true_indices) / true_indices) / positive_count
}

recompute_operating_metrics <- function(candidates, scores, source_metrics, padj_threshold) {
  conditions <- source_metrics %>%
    distinct(
      condition_id,
      treatment_target_true_chimera_fraction,
      control_target_true_chimera_fraction,
      background_fraction,
      total_reads,
      read_length
    ) %>%
    arrange(treatment_target_true_chimera_fraction, condition_id)

  call_sets <- c("RNAnue_native_padj", "edgeR", "DESeq2", "high_confidence_DIA")
  rows <- list()
  row_index <- 1

  for (condition in conditions$condition_id) {
    condition_meta <- conditions %>%
      filter(.data$condition_id == .env$condition) %>%
      slice(1)
    condition_candidates <- candidates %>%
      filter(.data$condition_id == .env$condition) %>%
      mutate(truth_positive = as.logical(truth_positive))
    labels <- stats::setNames(condition_candidates$truth_positive, condition_candidates$interaction_id)
    candidate_ids <- names(labels)
    truth_ids <- candidate_ids[labels]
    candidate_count <- length(candidate_ids)
    truth_count <- length(truth_ids)

    ambiguous_overlap_count <- source_metrics %>%
      filter(.data$condition_id == .env$condition) %>%
      summarise(value = suppressWarnings(max(ambiguous_overlap_count, na.rm = TRUE)), .groups = "drop") %>%
      pull(value)
    if (length(ambiguous_overlap_count) == 0 || !is.finite(ambiguous_overlap_count)) {
      ambiguous_overlap_count <- sum(condition_candidates$ambiguous_truth_overlap %in% TRUE, na.rm = TRUE)
    }

    for (method in call_sets) {
      method_scores <- scores %>%
        filter(.data$condition_id == .env$condition, .data$method == .env$method) %>%
        mutate(
          padj = as.numeric(padj),
          log2FoldChange = as.numeric(log2FoldChange),
          score_value = score_from_padj(score_padj)
        ) %>%
        arrange(padj, interaction_id) %>%
        distinct(interaction_id, .keep_all = TRUE)

      callable <- method_scores %>%
        filter(!is.na(padj), is.finite(padj), padj >= 0, padj <= 1) %>%
        filter(.env$method == "RNAnue_native_padj" | (!is.na(log2FoldChange) & log2FoldChange > 0))
      called_ids <- callable %>%
        filter(padj <= .env$padj_threshold) %>%
        pull(interaction_id) %>%
        intersect(candidate_ids)

      tp <- length(intersect(called_ids, truth_ids))
      fp <- length(setdiff(called_ids, truth_ids))
      fn <- truth_count - tp
      called_count <- length(called_ids)
      precision <- if (called_count > 0) tp / called_count else NA_real_
      recall <- if (truth_count > 0) tp / truth_count else NA_real_
      fdr <- if (called_count > 0) fp / called_count else NA_real_
      f1 <- if (!is.na(precision + recall) && (precision + recall) > 0) {
        2 * precision * recall / (precision + recall)
      } else {
        NA_real_
      }

      ranking_scores <- method_scores %>%
        filter(.env$method == "RNAnue_native_padj" | (!is.na(log2FoldChange) & log2FoldChange > 0))
      score_values <- ranking_scores$score_value
      names(score_values) <- ranking_scores$interaction_id
      auprc <- average_precision(labels, score_values)

      result_candidate_count <- if (method == "RNAnue_native_padj") {
        candidate_count
      } else if (method == "high_confidence_DIA") {
        NA_integer_
      } else {
        length(intersect(unique(method_scores$interaction_id), candidate_ids))
      }
      count_filter_loss <- if (is.na(result_candidate_count)) {
        NA_integer_
      } else {
        candidate_count - result_candidate_count
      }

      rows[[row_index]] <- data.frame(
        condition_id = condition,
        condition_label = condition_display(
          condition,
          condition_meta$treatment_target_true_chimera_fraction[[1]]
        ),
        treatment_target_true_chimera_fraction = condition_meta$treatment_target_true_chimera_fraction[[1]],
        control_target_true_chimera_fraction = condition_meta$control_target_true_chimera_fraction[[1]],
        background_fraction = condition_meta$background_fraction[[1]],
        total_reads = condition_meta$total_reads[[1]],
        read_length = condition_meta$read_length[[1]],
        call_set = method,
        call_set_label = label_call_set(method),
        status = "ok",
        candidate_count = candidate_count,
        true_chimera_positive_candidate_count = truth_count,
        truth_positive_candidate_count = truth_count,
        called_count = called_count,
        tp = tp,
        fp = fp,
        fn = fn,
        precision = precision,
        recall = recall,
        fdr = fdr,
        f1 = f1,
        auprc = auprc,
        count_filter_loss = count_filter_loss,
        ambiguous_overlap_count = ambiguous_overlap_count,
        padj_threshold = padj_threshold,
        rnanue_native_padj_max = padj_threshold,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1
    }
  }

  bind_rows(rows) %>%
    mutate(
      call_set_label = factor(call_set_label, levels = call_set_levels),
      condition_label = factor(condition_label, levels = unique(condition_label))
    )
}

make_threshold_curve <- function(candidates, scores, metrics, padj_threshold) {
  curve_conditions <- metrics %>%
    distinct(condition_id, condition_label, treatment_target_true_chimera_fraction)

  rows <- list()
  row_index <- 1

  for (condition in curve_conditions$condition_id) {
    condition_candidates <- candidates %>%
      filter(.data$condition_id == .env$condition) %>%
      select(condition_id, interaction_id, truth_positive)
    total_truth <- sum(condition_candidates$truth_positive, na.rm = TRUE)
    candidate_count <- nrow(condition_candidates)

    for (method in c("DESeq2", "edgeR", "RNAnue_native_padj", "high_confidence_DIA")) {
      method_scores <- scores %>%
        filter(.data$condition_id == .env$condition, .data$method == .env$method) %>%
        select(interaction_id, padj, log2FoldChange)

      candidate_scores <- condition_candidates %>%
        left_join(method_scores, by = "interaction_id")

      callable <- candidate_scores %>%
        filter(!is.na(padj), is.finite(padj), padj >= 0, padj <= 1) %>%
        filter(.env$method == "RNAnue_native_padj" | (!is.na(log2FoldChange) & log2FoldChange > 0)) %>%
        arrange(padj, interaction_id)

      threshold_tolerance <- max(sqrt(.Machine$double.eps), abs(padj_threshold) * sqrt(.Machine$double.eps))
      score_thresholds <- callable$padj[abs(callable$padj - padj_threshold) > threshold_tolerance]
      thresholds <- sort(unique(c(0, padj_threshold, 0.075, 0.1, 0.2, 0.5, 1, score_thresholds)))
      thresholds <- thresholds[is.finite(thresholds) & thresholds >= 0 & thresholds <= 1]

      if (nrow(callable) == 0) {
        curve <- data.frame(
          threshold = thresholds,
          tp = 0,
          fp = 0,
          called_count = 0
        )
      } else {
        curve <- data.frame(threshold = thresholds)
        curve$tp <- vapply(
          thresholds,
          function(threshold) sum(callable$truth_positive[callable$padj <= threshold], na.rm = TRUE),
          numeric(1)
        )
        curve$fp <- vapply(
          thresholds,
          function(threshold) sum(!callable$truth_positive[callable$padj <= threshold], na.rm = TRUE),
          numeric(1)
        )
        curve$called_count <- curve$tp + curve$fp
      }

      curve <- curve %>%
        mutate(
          condition_id = .env$condition,
          method = .env$method,
          method_label = label_call_set(method),
          candidate_count = candidate_count,
          truth_positive_candidate_count = total_truth,
          fn = truth_positive_candidate_count - tp,
          precision = ifelse(called_count > 0, tp / called_count, NA_real_),
          recall = ifelse(truth_positive_candidate_count > 0, tp / truth_positive_candidate_count, NA_real_),
          fdr = ifelse(called_count > 0, fp / called_count, NA_real_),
          f1 = ifelse(!is.na(precision + recall) & (precision + recall) > 0, 2 * precision * recall / (precision + recall), NA_real_),
          is_operating_threshold = abs(threshold - padj_threshold) <= threshold_tolerance
        )

      metric_operating <- metrics %>%
        filter(.data$condition_id == .env$condition, .data$call_set == .env$method) %>%
        slice(1)
      operating_index <- which(curve$is_operating_threshold)
      if (nrow(metric_operating) == 1 && length(operating_index) == 1) {
        curve$tp[operating_index] <- metric_operating$tp
        curve$fp[operating_index] <- metric_operating$fp
        curve$fn[operating_index] <- metric_operating$fn
        curve$called_count[operating_index] <- metric_operating$called_count
        curve$precision[operating_index] <- metric_operating$precision
        curve$recall[operating_index] <- metric_operating$recall
        curve$fdr[operating_index] <- metric_operating$fdr
        curve$f1[operating_index] <- metric_operating$f1
      }
      rows[[row_index]] <- curve
      row_index <- row_index + 1
    }
  }

  bind_rows(rows) %>%
    left_join(curve_conditions, by = "condition_id") %>%
    mutate(
      method_label = factor(method_label, levels = call_set_levels),
      condition_label = factor(condition_label, levels = unique(curve_conditions$condition_label)),
      true_chimera_axis_label = paste0(
        scales::percent(treatment_target_true_chimera_fraction, accuracy = 1),
        " True Chimera"
      )
    )
}

make_score_distribution <- function(scores, metrics) {
  condition_lookup <- metrics %>%
    distinct(condition_id, condition_label, treatment_target_true_chimera_fraction, padj_threshold)
  score_floor <- 1e-300

  scores %>%
    filter(!is.na(padj), is.finite(padj), padj >= 0, condition_id %in% condition_lookup$condition_id) %>%
    left_join(condition_lookup, by = "condition_id") %>%
    filter(treatment_target_true_chimera_fraction > 0) %>%
    mutate(
      neg_log10_padj = -log10(pmax(padj, score_floor)),
      neg_log10_padj_capped = pmin(neg_log10_padj, 50),
      truth_class = ifelse(truth_positive, "True Chimera positive", "Background candidate"),
      truth_class = factor(truth_class, levels = c("True Chimera positive", "Background candidate")),
      condition_label = factor(condition_label, levels = unique(condition_lookup$condition_label)),
      method_label = factor(method_label, levels = call_set_levels)
    )
}

make_truth_funnel <- function(metrics, summary_conditions) {
  summary_non_null <- summary_conditions %>%
    filter(treatment_target_true_chimera_fraction > 0) %>%
    select(
      condition_id,
      supported_truth_total,
      supported_truth_found_by_rnanue,
      generated_but_rnanue_missed,
      truth_positive_candidate_count
    )

  metric_counts <- metrics %>%
    filter(condition_id %in% summary_non_null$condition_id) %>%
    select(condition_id, call_set, tp)

  funnel_rows <- summary_non_null %>%
    transmute(condition_id, step = "Supported True Chimera", count = supported_truth_total) %>%
    bind_rows(
      summary_non_null %>%
        transmute(condition_id, step = "RNAnue-found True Chimera", count = supported_truth_found_by_rnanue),
      summary_non_null %>%
        transmute(condition_id, step = "Candidate-labeled True Chimera", count = truth_positive_candidate_count),
      metric_counts %>%
        filter(call_set %in% c("RNAnue_native_padj", "edgeR", "DESeq2", "high_confidence_DIA")) %>%
        transmute(
          condition_id,
          step = paste0(label_call_set(call_set), " TP"),
          count = tp
        )
    )

  step_levels <- c(
    "Supported True Chimera",
    "RNAnue-found True Chimera",
    "Candidate-labeled True Chimera",
    "RNAnue native TP",
    "edgeR TP",
    "DESeq2 TP",
    "High-confidence DIA TP"
  )

  metrics %>%
    distinct(
      condition_id,
      condition_label,
      treatment_target_true_chimera_fraction,
      background_fraction,
      total_reads,
      read_length,
      padj_threshold
    ) %>%
    right_join(funnel_rows, by = "condition_id") %>%
    mutate(
      step = factor(step, levels = step_levels),
      condition_label = factor(condition_label, levels = unique(metrics$condition_label))
    )
}

write_summary_readme <- function(output_dir, tables, figures, args, summary_lines) {
  figure_descriptions <- figures %>%
    distinct(stem, description) %>%
    mutate(line = paste0("- `", stem, "`: ", description)) %>%
    pull(line)

  readme <- c(
    "# Differential interaction benchmark",
    "",
    paste0("Source run: `", args$benchmark_dir, "`"),
    paste0("Operating adjusted-p-value cutoff: `", format_cutoff(args$padj_threshold), "`."),
    "",
    "Generated artifacts:",
    "",
    figure_descriptions,
    "",
    "Call-set semantics:",
    "",
    "- `RNAnue_native_padj`: RNAnue treatment native adjusted p-value at or below the cutoff.",
    "- `edgeR`: edgeR adjusted p-value at or below the cutoff with positive log2 fold change.",
    "- `DESeq2`: DESeq2 adjusted p-value at or below the cutoff with positive log2 fold change.",
    "- `high_confidence_DIA`: the report high-confidence gate, additionally evaluated at the folder cutoff.",
    "",
    summary_lines
  )
  writeLines(readme, file.path(output_dir, "README.md"))
}

summarise_best_recall <- function(metrics) {
  boosted <- metrics %>%
    filter(treatment_target_true_chimera_fraction > 0, !is.na(recall), !is.na(precision)) %>%
    arrange(desc(recall), desc(precision)) %>%
    slice(1)

  if (nrow(boosted) == 0) {
    return(character())
  }

  paste0(
    "Best recall at this cutoff is `",
    boosted$call_set_label[[1]],
    "` in `",
    boosted$condition_id[[1]],
    "` with recall ",
    scales::percent(boosted$recall[[1]], accuracy = 0.01),
    " and precision ",
    scales::percent(boosted$precision[[1]], accuracy = 0.01),
    "."
  )
}

summarise_key_metrics <- function(metrics) {
  boosted <- metrics %>%
    filter(treatment_target_true_chimera_fraction > 0) %>%
    mutate(label = label_call_set(call_set))

  method_lines <- boosted %>%
    filter(call_set %in% names(call_set_labels)) %>%
    group_by(label) %>%
    summarise(
      precision = paste(range_percent(precision), collapse = " to "),
      recall = paste(range_percent(recall), collapse = " to "),
      .groups = "drop"
    ) %>%
    arrange(match(label, call_set_levels)) %>%
    transmute(line = paste0("- ", label, ": precision ", precision, "; recall ", recall, "."))

  null_lines <- metrics %>%
    filter(treatment_target_true_chimera_fraction == 0) %>%
    mutate(label = label_call_set(call_set)) %>%
    arrange(match(label, call_set_levels)) %>%
    transmute(line = paste0("- Null false positives at the operating threshold, ", label, ": ", called_count, "."))

  c(summarise_best_recall(metrics), "", method_lines$line, null_lines$line)
}

range_percent <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) {
    return("NA")
  }
  scales::percent(range(x), accuracy = 0.1)
}

plot_operating_metrics <- function(metrics) {
  plot_df <- metrics %>%
    filter(treatment_target_true_chimera_fraction > 0) %>%
    select(condition_label, call_set_label, precision, recall, f1, auprc) %>%
    pivot_longer(c(precision, recall, f1, auprc), names_to = "metric", values_to = "value") %>%
    mutate(metric_label = factor(metric_labels[metric], levels = unname(metric_labels)))

  ggplot(plot_df, aes(x = condition_label, y = value, color = call_set_label, group = call_set_label)) +
    geom_line(linewidth = 0.65, alpha = 0.85) +
    geom_point(size = 2.4) +
    facet_wrap(~ metric_label, ncol = 2) +
    scale_color_manual(values = method_colors) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    benchmark_theme() +
    labs(
      x = "Treatment True Chimera fraction",
      y = "Value",
      color = "Call set",
      title = "Differential interaction benchmark metrics",
      subtitle = benchmark_subtitle(metrics)
    )
}

plot_confusion_counts <- function(metrics) {
  plot_df <- metrics %>%
    filter(treatment_target_true_chimera_fraction > 0) %>%
    select(condition_label, call_set_label, tp, fp, fn) %>%
    pivot_longer(c(tp, fp, fn), names_to = "count_type", values_to = "count") %>%
    mutate(
      count_type = factor(
        count_type,
        levels = c("tp", "fp", "fn"),
        labels = c("True positives", "False positives", "False negatives")
      )
    )

  ggplot(plot_df, aes(x = call_set_label, y = count, fill = call_set_label)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = scales::comma(count)), vjust = -0.2, size = 2.8) +
    facet_grid(count_type ~ condition_label, scales = "free_y") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.18))) +
    benchmark_theme(base_size = 10) +
    labs(
      x = NULL,
      y = "Interactions",
      title = "True-positive, false-positive, and false-negative counts",
      subtitle = "Panels use independent y scales so rare false positives remain visible"
    )
}

plot_null_false_positives <- function(metrics) {
  plot_df <- metrics %>%
    filter(treatment_target_true_chimera_fraction == 0) %>%
    mutate(false_positive_calls = called_count)

  ggplot(plot_df, aes(x = call_set_label, y = false_positive_calls, fill = call_set_label)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = scales::comma(false_positive_calls)), vjust = -0.2, size = 3.2) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    benchmark_theme() +
    labs(
      x = NULL,
      y = "False-positive calls",
      title = "Null-condition false-positive burden",
      subtitle = "The null simulation contains no supported True Chimera interactions"
    )
}

plot_precision_recall_curve <- function(curve_df, padj_threshold) {
  curve_plot_df <- curve_df %>%
    filter(called_count > 0, !is.na(precision), !is.na(recall)) %>%
    distinct(condition_id, method_label, condition_label, recall, precision, .keep_all = TRUE)

  operating <- curve_df %>%
    filter(is_operating_threshold)

  ggplot(curve_plot_df, aes(x = recall, y = precision, color = method_label)) +
    geom_step(linewidth = 0.7, alpha = 0.9, direction = "vh") +
    geom_point(data = operating, aes(x = recall, y = precision), size = 2.6, stroke = 0.6) +
    facet_wrap(~ condition_label, ncol = 2) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.02))) +
    benchmark_theme() +
    labs(
      x = "Recall",
      y = "Precision",
      color = "Call set",
      title = "Threshold precision-recall tradeoff",
      subtitle = paste0("Points mark the configured padj threshold of ", padj_threshold, "; DIA call sets require positive log2 fold change")
    )
}

plot_score_separation <- function(score_distribution) {
  plot_df <- score_distribution %>%
    filter(method %in% c("edgeR", "DESeq2"), !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    arrange(truth_class == "True Chimera positive")
  threshold <- unique(plot_df$padj_threshold)
  threshold <- threshold[!is.na(threshold) & is.finite(threshold) & threshold > 0]
  threshold <- if (length(threshold) == 0) NA_real_ else threshold[[1]]

  ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj_capped, color = truth_class)) +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_point(alpha = 0.42, size = 0.38, stroke = 0) +
    facet_grid(method_label ~ condition_label) +
    scale_color_manual(values = truth_colors, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 50, 10), limits = c(0, 50), expand = expansion(mult = c(0.01, 0.03))) +
    benchmark_theme(base_size = 10) +
    labs(
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Candidate class",
      title = "DIA score direction separates True Chimera-positive and background candidates",
      subtitle = "padj < 1e-50 floored."
    )
}

plot_truth_funnel <- function(funnel_df) {
  totals <- funnel_df %>%
    filter(step == "Candidate-labeled True Chimera") %>%
    transmute(condition_id, total_true_chimera_candidates = count)

  plot_df <- funnel_df %>%
    filter(grepl(" TP$", as.character(step))) %>%
    left_join(totals, by = "condition_id") %>%
    filter(!is.na(total_true_chimera_candidates), total_true_chimera_candidates > 0) %>%
    mutate(
      call_set_label = sub(" TP$", "", as.character(step)),
      call_set_label = factor(call_set_label, levels = call_set_levels),
      recovered = pmin(count, total_true_chimera_candidates),
      missed = pmax(total_true_chimera_candidates - recovered, 0)
    ) %>%
    select(
      condition_id,
      condition_label,
      treatment_target_true_chimera_fraction,
      background_fraction,
      total_reads,
      read_length,
      padj_threshold,
      call_set_label,
      total_true_chimera_candidates,
      recovered,
      missed
    ) %>%
    pivot_longer(c(recovered, missed), names_to = "status", values_to = "count") %>%
    mutate(
      status = factor(ifelse(status == "recovered", "Recovered", "Missed"), levels = c("Missed", "Recovered")),
      fraction = count / total_true_chimera_candidates
    )

  ggplot(plot_df, aes(x = call_set_label, y = fraction, fill = status)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    facet_wrap(~ condition_label, nrow = 1) +
    scale_fill_manual(
      values = c("Recovered" = "#238B35", "Missed" = "#B8B8B8"),
      breaks = c("Recovered", "Missed")
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    benchmark_theme(base_size = 10) +
    labs(
      x = "Call set",
      y = "True Chimera-positive candidates",
      fill = "Status",
      title = "True Chimera recovery by call set",
      subtitle = benchmark_subtitle(plot_df)
    )
}

write_input_copies <- function(paths, output_dir) {
  input_dir <- file.path(output_dir, "inputs")
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  copied <- vapply(paths, function(path) {
    file.copy(path, file.path(input_dir, basename(path)), overwrite = TRUE)
  }, logical(1))
  if (any(!copied)) {
    warning("Some input files could not be copied into ", input_dir, call. = FALSE)
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  tables <- read_benchmark_tables(args$benchmark_dir)
  scores <- read_method_scores(args$benchmark_dir, tables$candidates)
  tables$metrics <- recompute_operating_metrics(
    tables$candidates,
    scores,
    tables$metrics,
    args$padj_threshold
  )
  threshold_curve <- make_threshold_curve(
    tables$candidates,
    scores,
    tables$metrics,
    args$padj_threshold
  )
  score_distribution <- make_score_distribution(scores, tables$metrics)
  score_summary <- score_distribution %>%
    filter(method %in% c("edgeR", "DESeq2", "RNAnue_native_padj")) %>%
    mutate(
      true_chimera_axis_label = as.character(condition_label),
      candidate_class = ifelse(truth_class == "True Chimera positive", "true_chimera_positive", "background_candidate"),
      score_source = ifelse(method == "RNAnue_native_padj", "RNAnue_native_padj", paste0(method, "_padj"))
    ) %>%
    group_by(condition_id, true_chimera_axis_label, candidate_class, score_source) %>%
    summarise(
      n = n(),
      median_padj = median(padj, na.rm = TRUE),
      q25_padj = quantile(padj, 0.25, na.rm = TRUE),
      q75_padj = quantile(padj, 0.75, na.rm = TRUE),
      median_neg_log10_padj = median(neg_log10_padj, na.rm = TRUE),
      q25_neg_log10_padj = quantile(neg_log10_padj, 0.25, na.rm = TRUE),
      q75_neg_log10_padj = quantile(neg_log10_padj, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  truth_funnel <- make_truth_funnel(tables$metrics, tables$summary$conditions)

  metric_scoreboard <- tables$metrics %>%
    select(
      condition_id,
      treatment_target_true_chimera_fraction,
      control_target_true_chimera_fraction,
      background_fraction,
      total_reads,
      read_length,
      call_set,
      call_set_label,
      status,
      candidate_count,
      true_chimera_positive_candidate_count,
      called_count,
      tp,
      fp,
      fn,
      precision,
      recall,
      fdr,
      f1,
      auprc,
      count_filter_loss,
      ambiguous_overlap_count,
      padj_threshold,
      rnanue_native_padj_max
    )

  threshold_curve_output <- threshold_curve %>%
    transmute(
      condition_id,
      treatment_target_true_chimera_fraction,
      true_chimera_axis_label,
      call_set = method,
      call_set_label = as.character(method_label),
      threshold,
      called_count,
      tp,
      fp,
      precision,
      recall,
      is_operating_threshold
    )

  readr::write_tsv(metric_scoreboard, file.path(args$output_dir, "metric_scoreboard.tsv"), na = "")
  readr::write_tsv(threshold_curve_output, file.path(args$output_dir, "threshold_precision_recall_curve.tsv"), na = "")
  readr::write_tsv(score_summary, file.path(args$output_dir, "score_distribution_summary.tsv"), na = "")

  figure_specs <- list(
    list(
      stem = "01_operating_metrics",
      description = "Precision, recall, F1, and AUPRC on a fixed 0-100% scale.",
      plot = plot_operating_metrics(tables$metrics),
      width = 11,
      height = 6.4
    ),
    list(
      stem = "02_tp_fp_fn_counts",
      description = "True-positive, false-positive, and false-negative counts for boosted conditions.",
      plot = plot_confusion_counts(tables$metrics),
      width = 11,
      height = 8
    ),
    list(
      stem = "03_null_false_positives",
      description = "False-positive calls in the null condition.",
      plot = plot_null_false_positives(tables$metrics),
      width = 8.5,
      height = 5.4
    ),
    list(
      stem = "04_score_separation",
      description = "edgeR and DESeq2 scatter plots of log2 fold change versus `-log10(adjusted p-value)`.",
      plot = plot_score_separation(score_distribution),
      width = 11,
      height = 7.2
    ),
    list(
      stem = "05_threshold_precision_recall",
      description = "Precision-recall curves with the selected cutoff marked.",
      plot = plot_precision_recall_curve(threshold_curve, args$padj_threshold),
      width = 10.5,
      height = 5.6
    ),
    list(
      stem = "06_true_chimera_recovery_funnel",
      description = "Recovered versus missed True Chimera-positive candidates on a fixed 0-100% scale.",
      plot = plot_truth_funnel(truth_funnel),
      width = 11,
      height = 5.8
    )
  )

  figure_rows <- lapply(figure_specs, function(spec) {
    paths <- save_plot(
      spec$plot,
      args$output_dir,
      spec$stem,
      width = spec$width,
      height = spec$height,
      formats = args$formats,
      dpi = args$dpi
    )
    data.frame(
      stem = spec$stem,
      file = basename(paths),
      description = spec$description,
      stringsAsFactors = FALSE
    )
  })
  figure_index <- bind_rows(figure_rows)

  write_summary_readme(
    args$output_dir,
    tables,
    figure_index,
    args,
    summarise_key_metrics(tables$metrics)
  )

  cat(args$output_dir, "\n")
  invisible(args$output_dir)
}

if (sys.nframe() == 0) {
  main()
}
