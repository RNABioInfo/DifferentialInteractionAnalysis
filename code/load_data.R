read_metadata <- function(metadata_file_path, params) {
  metadata_df <- as.data.frame(
    readr::read_csv(
      metadata_file_path,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    ),
    stringsAsFactors = FALSE
  )

  if (!"sample_id" %in% names(metadata_df) && "name" %in% names(metadata_df)) {
    metadata_df$sample_id <- metadata_df$name
  }

  if (!"role" %in% names(metadata_df)) {
    condition_levels <- unique(metadata_df$condition)
    if (length(condition_levels) != 2) {
      stop(
        "Legacy metadata without a role column can only be inferred for exactly two conditions.",
        call. = FALSE
      )
    }
    warning(
      "Metadata has no role column. Inferring first condition as ligation_control ",
      "and second condition as treatment. Add explicit roles for production use.",
      call. = FALSE
    )
    metadata_df$role <- ifelse(
      metadata_df$condition == condition_levels[1],
      params$control_role,
      params$case_role
    )
  }

  required <- c("sample_id", "condition", "role")
  missing_required <- setdiff(required, names(metadata_df))
  if (length(missing_required) > 0) {
    stop(
      "Metadata is missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"include" %in% names(metadata_df)) {
    metadata_df$include <- ifelse(
      metadata_df$role %in% c(params$case_role, params$control_role),
      "model",
      "qc"
    )
  }

  metadata_df$sample_id <- as.character(metadata_df$sample_id)
  metadata_df$condition <- as.character(metadata_df$condition)
  metadata_df$role <- as.character(metadata_df$role)
  metadata_df$include <- tolower(as.character(metadata_df$include))
  metadata_df$include[is.na(metadata_df$include) | metadata_df$include == ""] <- "model"
  metadata_df$include <- dplyr::case_when(
    metadata_df$include %in% c("true", "yes", "include") ~ "model",
    metadata_df$include %in% c("false", "no", "exclude") ~ "exclude",
    metadata_df$include %in% c("model", "qc") ~ metadata_df$include,
    TRUE ~ metadata_df$include
  )

  valid_roles <- c("treatment", "ligation_control", "no_ligation_control", "other_qc")
  invalid_roles <- setdiff(unique(metadata_df$role), valid_roles)
  if (length(invalid_roles) > 0) {
    stop(
      "Metadata contains invalid roles: ",
      paste(invalid_roles, collapse = ", "),
      ". Valid roles are: ",
      paste(valid_roles, collapse = ", "),
      call. = FALSE
    )
  }

  duplicate_samples <- metadata_df$sample_id[duplicated(metadata_df$sample_id)]
  if (length(duplicate_samples) > 0) {
    stop(
      "Metadata contains duplicate sample_id values: ",
      paste(unique(duplicate_samples), collapse = ", "),
      call. = FALSE
    )
  }

  rownames(metadata_df) <- metadata_df$sample_id
  metadata_df
}

read_counts <- function(counts_path) {
  if (!file.exists(counts_path)) {
    stop("Interaction count file does not exist: ", counts_path, call. = FALSE)
  }

  counts_raw <- readr::read_delim(
    counts_path,
    delim = "\t",
    skip = 2,
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )

  if (!all(c("Name", "Description") %in% names(counts_raw))) {
    stop("GCT count table must contain Name and Description columns.", call. = FALSE)
  }

  sample_cols <- setdiff(names(counts_raw), c("Name", "Description"))
  counts_df <- as.data.frame(
    lapply(counts_raw[sample_cols], as.numeric),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(counts_df) <- counts_raw$Name
  counts_df
}

read_interactions <- function(interactions_path) {
  if (!file.exists(interactions_path)) {
    stop("Interaction BEDPE file does not exist: ", interactions_path, call. = FALSE)
  }

  bedpe <- as.data.frame(
    readr::read_delim(
      interactions_path,
      delim = "\t",
      skip = 2,
      col_names = FALSE,
      name_repair = "minimal",
      show_col_types = FALSE,
      progress = FALSE
    ),
    stringsAsFactors = FALSE
  )

  expected_names <- c(
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "interaction_id", "score", "strand1", "strand2", "cluster_ids",
    "arm1_features", "arm2_features", "arm1_count", "arm2_count", "color"
  )
  names(bedpe)[seq_len(min(length(expected_names), ncol(bedpe)))] <-
    expected_names[seq_len(min(length(expected_names), ncol(bedpe)))]

  if (!"interaction_id" %in% names(bedpe)) {
    stop("BEDPE table must contain an interaction ID column in position 7.", call. = FALSE)
  }

  bedpe$interaction_id <- as.character(bedpe$interaction_id)
  rownames(bedpe) <- bedpe$interaction_id
  bedpe
}

parse_feature_field <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "" || x == ".") {
    return(character())
  }
  tokens <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  ids <- sub("^[^:]+:", "", tokens)
  ids <- unlist(strsplit(ids, ",", fixed = TRUE), use.names = FALSE)
  ids <- trimws(ids)
  unique(ids[nzchar(ids) & !ids %in% c(".", "*")])
}

split_gff_identifier_values <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) {
    return(character())
  }
  values <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  values <- trimws(values)
  unique(values[nzchar(values) & !values %in% c(".", "*")])
}

read_rnanue_annotation_ids <- function(gff_path) {
  unavailable <- function(reason) {
    list(ids = character(), available = FALSE, reason = reason)
  }

  if (!is_present(gff_path)) {
    return(unavailable("annotations_file is missing or not configured"))
  }
  if (!file.exists(gff_path)) {
    return(unavailable(paste("annotations_file does not exist:", gff_path)))
  }

  lines <- tryCatch(
    {
      connection <- if (grepl("\\.gz$", gff_path, ignore.case = TRUE)) {
        gzfile(gff_path, open = "rt")
      } else {
        file(gff_path, open = "rt")
      }
      on.exit(close(connection), add = TRUE)
      readLines(connection, warn = FALSE)
    },
    error = function(e) e
  )
  if (inherits(lines, "error")) {
    return(unavailable(paste("annotations_file could not be read:", conditionMessage(lines))))
  }

  lines <- lines[nzchar(trimws(lines)) & !grepl("^#", lines)]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 9]
  if (length(fields) == 0) {
    return(unavailable("annotations_file contained no usable nine-column GFF rows"))
  }

  identifiers <- character()
  for (field in fields) {
    attributes <- unlist(strsplit(field[9], ";", fixed = TRUE), use.names = FALSE)
    attributes <- trimws(attributes)
    for (attribute in attributes) {
      key_value <- strsplit(attribute, "=", fixed = TRUE)[[1]]
      if (length(key_value) < 2) {
        next
      }
      key <- trimws(key_value[1])
      if (!key %in% c("ID", "Parent")) {
        next
      }
      value <- paste(key_value[-1], collapse = "=")
      value <- tryCatch(utils::URLdecode(value), error = function(e) value)
      identifiers <- c(identifiers, split_gff_identifier_values(value))
    }
  }
  identifiers <- unique(identifiers)
  if (length(identifiers) == 0) {
    return(unavailable("annotations_file contained no ID or Parent attributes"))
  }

  list(ids = identifiers, available = TRUE, reason = "")
}

parse_sample_feature_map <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "" || x == ".") {
    return(list(by_sample = list(), all_ids = character()))
  }

  tokens <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  by_sample <- list()
  all_ids <- character()

  for (token in tokens) {
    separator <- regexpr(":", token, fixed = TRUE)[1]
    if (separator > 0) {
      sample_id <- substr(token, 1, separator - 1)
      id_text <- substr(token, separator + 1, nchar(token))
    } else {
      sample_id <- NA_character_
      id_text <- token
    }

    ids <- unlist(strsplit(id_text, ",", fixed = TRUE), use.names = FALSE)
    ids <- trimws(ids)
    ids <- unique(ids[nzchar(ids) & !ids %in% c(".", "*")])
    all_ids <- c(all_ids, ids)

    if (!is.na(sample_id) && length(ids) > 0) {
      by_sample[[sample_id]] <- unique(c(by_sample[[sample_id]], ids))
    }
  }

  list(by_sample = by_sample, all_ids = unique(all_ids))
}

feature_ids_for_sample <- function(feature_map, sample_id, stable_feature_ids = character()) {
  sample_ids <- feature_map$by_sample[[sample_id]]

  if (length(sample_ids) > 0) {
    return(sample_ids)
  }

  intersect(feature_map$all_ids, stable_feature_ids)
}

parse_feature_field_for_sample <- function(x, sample_id, stable_feature_ids = character()) {
  feature_ids_for_sample(parse_sample_feature_map(x), sample_id, stable_feature_ids)
}

parse_cluster_tokens <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "" || x == ".") {
    return(data.frame(sample_id = character(), cluster_ID = character()))
  }
  tokens <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  rows <- lapply(tokens, function(token) {
    sample_id <- trimws(sub(":.*$", "", token))
    cluster_ids <- sub("^[^:]+:", "", token)
    cluster_ids <- unlist(strsplit(cluster_ids, ",", fixed = TRUE), use.names = FALSE)
    cluster_ids <- trimws(cluster_ids)
    cluster_ids <- cluster_ids[nzchar(cluster_ids)]
    if (!nzchar(sample_id) || length(cluster_ids) == 0) {
      return(NULL)
    }
    data.frame(
      sample_id = rep(sample_id, length(cluster_ids)),
      cluster_ID = cluster_ids,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(sample_id = character(), cluster_ID = character()))
  }
  parsed <- do.call(rbind, rows)
  rownames(parsed) <- NULL
  parsed
}

find_sample_file <- function(base_dir, sample_id, suffix) {
  if (is.na(base_dir) || !dir.exists(base_dir)) {
    return(NA_character_)
  }
  pattern <- paste0("^", sample_id, suffix, "$")
  files <- list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    return(NA_character_)
  }
  files[1]
}

read_read_count_summaries <- function(sample_ids, detect_dir) {
  rows <- lapply(sample_ids, function(sample_id) {
    file_path <- find_sample_file(detect_dir, sample_id, "_read_counts_summary.tsv")
    if (is.na(file_path)) {
      return(data.frame(
        sample_id = sample_id,
        splits = NA_real_,
        singletons = NA_real_,
        mapped_reads = NA_real_,
        read_count_summary_path = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    summary_df <- readr::read_tsv(file_path, show_col_types = FALSE, progress = FALSE)
    mapped_reads <- if ("mapped_reads" %in% names(summary_df)) {
      as.numeric(summary_df$mapped_reads[1])
    } else {
      NA_real_
    }
    data.frame(
      sample_id = sample_id,
      splits = as.numeric(summary_df$splits[1]),
      singletons = as.numeric(summary_df$singletons[1]),
      mapped_reads = mapped_reads,
      read_count_summary_path = file_path,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

read_interaction_transcript_counts <- function(sample_ids, analyze_dir) {
  counts <- list()
  for (sample_id in sample_ids) {
    file_path <- find_sample_file(analyze_dir, sample_id, "_interaction_transcript_counts.tsv")
    if (is.na(file_path)) {
      counts[[sample_id]] <- NULL
      next
    }
    transcript_counts <- readr::read_delim(
      file_path,
      delim = "\t",
      col_names = c("feature_id", "count"),
      show_col_types = FALSE,
      progress = FALSE
    )
    values <- as.numeric(transcript_counts$count)
    names(values) <- as.character(transcript_counts$feature_id)
    counts[[sample_id]] <- values
  }
  counts
}

read_native_interactions <- function(sample_ids, analyze_dir) {
  native <- list()
  for (sample_id in sample_ids) {
    file_path <- find_sample_file(analyze_dir, sample_id, "_interactions.tsv")
    if (is.na(file_path)) {
      native[[sample_id]] <- NULL
      next
    }
    native_df <- as.data.frame(
      readr::read_tsv(file_path, show_col_types = FALSE, progress = FALSE),
      stringsAsFactors = FALSE
    )
    if ("cluster_ID" %in% names(native_df)) {
      rownames(native_df) <- native_df$cluster_ID
    }
    native[[sample_id]] <- native_df
  }
  native
}

build_rnanue_metric_summary <- function(bedpe, native_interactions) {
  rows <- lapply(seq_len(nrow(bedpe)), function(i) {
    clusters <- parse_cluster_tokens(bedpe$cluster_ids[i])
    hits <- list()
    if (nrow(clusters) > 0) {
      for (j in seq_len(nrow(clusters))) {
        sample_df <- native_interactions[[clusters$sample_id[j]]]
        if (!is.null(sample_df) && clusters$cluster_ID[j] %in% rownames(sample_df)) {
          hits[[length(hits) + 1]] <- sample_df[clusters$cluster_ID[j], , drop = FALSE]
        }
      }
    }

    if (length(hits) == 0) {
      return(data.frame(
        interaction_id = bedpe$interaction_id[i],
        rnanue_min_p_value = NA_real_,
        rnanue_min_padj_value = NA_real_,
        rnanue_mean_support_per_effective_bp = NA_real_,
        rnanue_median_arm_balance = NA_real_,
        rnanue_max_coverage_components = NA_real_,
        rnanue_coverage_profiles = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    hit_df <- dplyr::bind_rows(hits)
    p_value <- suppressWarnings(as.numeric(hit_df$p_value))
    padj_value <- suppressWarnings(as.numeric(hit_df$padj_value))
    support <- suppressWarnings(as.numeric(hit_df$support_per_effective_bp))
    arm_balance <- suppressWarnings(as.numeric(hit_df$arm_balance))
    components <- suppressWarnings(as.numeric(hit_df$coverage_components))
    profiles <- unique(as.character(hit_df$coverage_profile))
    profiles <- profiles[!is.na(profiles) & nzchar(profiles)]

    data.frame(
      interaction_id = bedpe$interaction_id[i],
      rnanue_min_p_value = suppressWarnings(min(p_value, na.rm = TRUE)),
      rnanue_min_padj_value = suppressWarnings(min(padj_value, na.rm = TRUE)),
      rnanue_mean_support_per_effective_bp = suppressWarnings(mean(support, na.rm = TRUE)),
      rnanue_median_arm_balance = suppressWarnings(stats::median(arm_balance, na.rm = TRUE)),
      rnanue_max_coverage_components = suppressWarnings(max(components, na.rm = TRUE)),
      rnanue_coverage_profiles = if (length(profiles) == 0) NA_character_ else paste(profiles, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  metric_df <- dplyr::bind_rows(rows)
  numeric_cols <- setdiff(names(metric_df), c("interaction_id", "rnanue_coverage_profiles"))
  for (col in numeric_cols) {
    metric_df[[col]][is.infinite(metric_df[[col]])] <- NA_real_
  }
  metric_df
}

format_rnanue_padj_values <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(trimws(formatC(values, digits = 6, format = "g")), collapse = ",")
}

build_rnanue_padj_by_role <- function(bedpe, native_interactions, metadata, params) {
  case_samples <- metadata$sample_id[metadata$role == params$case_role]
  control_samples <- metadata$sample_id[metadata$role == params$control_role]

  collapse_role_values <- function(cluster_rows, role_samples) {
    entries <- vapply(role_samples, function(sample_id) {
      sample_clusters <- cluster_rows$cluster_ID[cluster_rows$sample_id == sample_id]
      sample_df <- native_interactions[[sample_id]]
      if (length(sample_clusters) == 0 || is.null(sample_df) || !"padj_value" %in% names(sample_df)) {
        return(NA_character_)
      }
      sample_clusters <- sample_clusters[sample_clusters %in% rownames(sample_df)]
      if (length(sample_clusters) == 0) {
        return(NA_character_)
      }
      values <- format_rnanue_padj_values(sample_df[sample_clusters, "padj_value", drop = TRUE])
      if (is.na(values)) {
        return(NA_character_)
      }
      paste0(sample_id, ":", values)
    }, character(1))
    entries <- entries[!is.na(entries) & nzchar(entries)]
    if (length(entries) == 0) {
      return(NA_character_)
    }
    paste(entries, collapse = ";")
  }

  rows <- lapply(seq_len(nrow(bedpe)), function(i) {
    clusters <- parse_cluster_tokens(bedpe$cluster_ids[i])
    data.frame(
      interaction_id = bedpe$interaction_id[i],
      rnanue_condition_padj_values = collapse_role_values(clusters, case_samples),
      rnanue_control_padj_values = collapse_role_values(clusters, control_samples),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

build_sample_qc <- function(metadata, counts_all, read_summaries) {
  sample_ids <- metadata$sample_id
  available_samples <- intersect(sample_ids, colnames(counts_all))
  interaction_library_sum <- rep(NA_real_, length(sample_ids))
  detected_interactions <- rep(NA_integer_, length(sample_ids))
  names(interaction_library_sum) <- sample_ids
  names(detected_interactions) <- sample_ids
  interaction_library_sum[available_samples] <- colSums(counts_all[, available_samples, drop = FALSE])
  detected_interactions[available_samples] <- colSums(counts_all[, available_samples, drop = FALSE] > 0)

  qc <- dplyr::left_join(metadata, read_summaries, by = "sample_id")
  qc$interaction_library_sum <- unname(interaction_library_sum[qc$sample_id])
  qc$detected_interactions <- unname(detected_interactions[qc$sample_id])
  rnanue_classified_reads <- qc$splits + qc$singletons
  has_explicit_mapped_reads <- !is.na(qc$mapped_reads)
  qc$mapped_read_count <- ifelse(has_explicit_mapped_reads, qc$mapped_reads, rnanue_classified_reads)
  qc$mapped_read_source <- ifelse(
    has_explicit_mapped_reads,
    "mapped_reads_column",
    "rnanue_splits_plus_singletons"
  )
  qc$split_fraction <- qc$splits / rnanue_classified_reads
  qc
}

build_no_ligation_qc <- function(sample_qc, params) {
  no_ligation <- sample_qc[sample_qc$role %in% params$qc_roles, , drop = FALSE]
  if (nrow(no_ligation) == 0) {
    return(no_ligation)
  }
  ligation_reference <- sample_qc[sample_qc$role == params$control_role, , drop = FALSE]
  median_ligation_split_fraction <- stats::median(ligation_reference$split_fraction, na.rm = TRUE)
  median_ligation_interactions <- stats::median(ligation_reference$interaction_library_sum, na.rm = TRUE)
  no_ligation$split_fraction_vs_ligation_control_median <-
    no_ligation$split_fraction / median_ligation_split_fraction
  no_ligation$interaction_sum_vs_ligation_control_median <-
    no_ligation$interaction_library_sum / median_ligation_interactions
  no_ligation
}

sum_feature_counts <- function(feature_ids, sample_counts) {
  if (length(feature_ids) == 0 || is.null(sample_counts)) {
    return(NA_real_)
  }
  present <- intersect(feature_ids, names(sample_counts))
  if (length(present) == 0) {
    return(0)
  }
  sum(sample_counts[present], na.rm = TRUE)
}

build_offset_matrices <- function(counts_model, bedpe, read_summaries,
                                  interaction_transcript_counts, params) {
  sample_ids <- colnames(counts_model)
  pseudocount <- params$offset_pseudocount
  read_summaries <- read_summaries[match(sample_ids, read_summaries$sample_id), , drop = FALSE]
  exposure <- read_summaries$splits
  fallback_exposure <- colSums(counts_model)
  exposure[!is.finite(exposure) | exposure <= 0] <-
    fallback_exposure[!is.finite(exposure) | exposure <= 0]
  exposure[!is.finite(exposure) | exposure <= 0] <- 1
  names(exposure) <- sample_ids

  log_offset <- matrix(
    log(exposure),
    nrow = nrow(counts_model),
    ncol = ncol(counts_model),
    byrow = TRUE,
    dimnames = dimnames(counts_model)
  )
  pair_background_available <- matrix(FALSE, nrow(counts_model), ncol(counts_model), dimnames = dimnames(counts_model))
  unresolved_arm1 <- logical(nrow(counts_model))
  unresolved_arm2 <- logical(nrow(counts_model))
  names(unresolved_arm1) <- rownames(counts_model)
  names(unresolved_arm2) <- rownames(counts_model)

  stable_features <- if (isTRUE(params$pair_background)) {
    read_rnanue_annotation_ids(params$annotations_path)
  } else {
    list(ids = character(), available = FALSE, reason = "pair background disabled by configuration")
  }
  pair_background_applied <- isTRUE(params$pair_background) && stable_features$available

  if (isTRUE(params$pair_background) && !stable_features$available) {
    warning(
      "Pair-background normalization was requested but disabled: ",
      stable_features$reason,
      ". Exposure-only offsets will be used.",
      call. = FALSE
    )
  }

  if (pair_background_applied) {
    for (i in seq_len(nrow(counts_model))) {
      interaction_id <- rownames(counts_model)[i]
      bedpe_row <- bedpe[interaction_id, , drop = FALSE]
      arm1_feature_map <- parse_sample_feature_map(bedpe_row$arm1_features[1])
      arm2_feature_map <- parse_sample_feature_map(bedpe_row$arm2_features[1])
      for (sample_id in sample_ids) {
        arm1 <- feature_ids_for_sample(arm1_feature_map, sample_id, stable_features$ids)
        arm2 <- feature_ids_for_sample(arm2_feature_map, sample_id, stable_features$ids)
        if (length(arm1) == 0) {
          unresolved_arm1[i] <- TRUE
        }
        if (length(arm2) == 0) {
          unresolved_arm2[i] <- TRUE
        }
        if (length(arm1) == 0 || length(arm2) == 0) {
          next
        }

        sample_counts <- interaction_transcript_counts[[sample_id]]
        if (is.null(sample_counts) || length(sample_counts) == 0) {
          next
        }
        total_bg <- sum(sample_counts, na.rm = TRUE)
        if (!is.finite(total_bg) || total_bg <= 0) {
          next
        }
        arm1_sum <- sum_feature_counts(arm1, sample_counts)
        arm2_sum <- sum_feature_counts(arm2, sample_counts)
        multiplier <- if (length(intersect(arm1, arm2)) == 0) 2 else 1
        pair_weight <- multiplier *
          ((arm1_sum + pseudocount) / (total_bg + pseudocount)) *
          ((arm2_sum + pseudocount) / (total_bg + pseudocount))
        pair_weight <- max(pair_weight, .Machine$double.xmin)
        log_offset[interaction_id, sample_id] <-
          log_offset[interaction_id, sample_id] + log(pair_weight)
        pair_background_available[interaction_id, sample_id] <- TRUE
      }
    }
  }

  normalization_diagnostics <- data.frame(
    interaction_id = rownames(counts_model),
    pair_background_available_fraction = rowMeans(pair_background_available),
    unresolved_arm1 = unresolved_arm1[rownames(counts_model)],
    unresolved_arm2 = unresolved_arm2[rownames(counts_model)],
    min_log_offset = apply(log_offset, 1, min),
    max_log_offset = apply(log_offset, 1, max),
    stringsAsFactors = FALSE
  )

  normalization_setup_diagnostics <- data.frame(
    annotations_file = as.character(params$annotations_path),
    stable_feature_id_count = length(stable_features$ids),
    pair_background_requested = isTRUE(params$pair_background),
    pair_background_applied = pair_background_applied,
    fallback_reason = if (pair_background_applied) "" else stable_features$reason,
    stringsAsFactors = FALSE
  )

  list(
    log_offset = log_offset,
    pair_background_available = pair_background_available,
    sample_exposure = exposure,
    normalization_diagnostics = normalization_diagnostics,
    normalization_setup_diagnostics = normalization_setup_diagnostics
  )
}
