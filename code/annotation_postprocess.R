empty_annotation_diagnostics <- function(status = "disabled", message = "GFF annotation disabled.") {
  data.frame(
    metric = c("status", "message"),
    value = c(status, message),
    stringsAsFactors = FALSE
  )
}

split_param_ids <- function(x) {
  ids <- as_character_vector(x, character())
  ids <- unlist(strsplit(ids, ",", fixed = TRUE), use.names = FALSE)
  ids <- trimws(ids)
  unique(ids[nzchar(ids)])
}

strip_seq_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

empty_arm_annotation_table <- function() {
  data.frame(
    interaction_id = character(),
    arm = character(),
    chrom = character(),
    start = integer(),
    end = integer(),
    strand = character(),
    rnanue_arm_features = character(),
    selected_feature_id = character(),
    selected_feature_name = character(),
    selected_feature_type = character(),
    selected_partner_class = character(),
    feature_annotation_source = character(),
    partner_class_source = character(),
    gene = character(),
    locus_tag = character(),
    Alias = character(),
    Parent = character(),
    product = character(),
    gbkey = character(),
    gene_biotype = character(),
    overlap_bp = integer(),
    arm_overlap_fraction = numeric(),
    feature_overlap_fraction = numeric(),
    reciprocal_overlap_fraction = numeric(),
    annotation_status = character(),
    seqname_normalized = logical(),
    result_class = character(),
    high_confidence = logical(),
    stringsAsFactors = FALSE
  )
}

empty_partner_type_pairs <- function() {
  data.frame(
    interaction_id = character(),
    arm1_partner_class = character(),
    arm2_partner_class = character(),
    partner_class_1 = character(),
    partner_class_2 = character(),
    partner_class_pair = character(),
    result_class = character(),
    high_confidence = logical(),
    stringsAsFactors = FALSE
  )
}

empty_partner_type_confusion <- function() {
  data.frame(
    result_set = character(),
    partner_class_1 = character(),
    partner_class_2 = character(),
    interactions = integer(),
    stringsAsFactors = FALSE
  )
}

empty_target_gene_interactions <- function() {
  data.frame(
    target_id = character(),
    interaction_id = character(),
    target_arm = character(),
    target_partner_class = character(),
    other_arm_partner_class = character(),
    selected_feature_id = character(),
    selected_feature_name = character(),
    rnanue_arm_features = character(),
    result_class = character(),
    high_confidence = logical(),
    log2FoldChange = numeric(),
    padj = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_target_feature_matches <- function() {
  data.frame(
    target_id = character(),
    interaction_id = character(),
    target_arm = character(),
    target_partner_class = character(),
    other_arm_partner_class = character(),
    selected_feature_id = character(),
    selected_feature_name = character(),
    selected_feature_type = character(),
    gene = character(),
    locus_tag = character(),
    Alias = character(),
    Parent = character(),
    product = character(),
    gbkey = character(),
    gene_biotype = character(),
    rnanue_arm_features = character(),
    result_class = character(),
    high_confidence = logical(),
    overlap_bp = integer(),
    reciprocal_overlap_fraction = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_target_gene_bedpe <- function() {
  data.frame(
    chr1 = character(),
    start1 = integer(),
    end1 = integer(),
    chr2 = character(),
    start2 = integer(),
    end2 = integer(),
    interaction_id = character(),
    target_ids = character(),
    target_arms = character(),
    result_class = character(),
    high_confidence = logical(),
    stringsAsFactors = FALSE
  )
}

empty_gff_annotation_assessment <- function(status = "disabled", message = "GFF annotation disabled.") {
  list(
    arm_annotations = empty_arm_annotation_table(),
    partner_type_pairs = empty_partner_type_pairs(),
    partner_type_confusion = empty_partner_type_confusion(),
    target_gene_interactions = empty_target_gene_interactions(),
    target_gene_bedpe = empty_target_gene_bedpe(),
    diagnostics = empty_annotation_diagnostics(status, message)
  )
}

parse_gff_attributes <- function(attribute_string) {
  if (length(attribute_string) == 0 || is.na(attribute_string) || !nzchar(attribute_string)) {
    return(character())
  }
  parts <- unlist(strsplit(attribute_string, ";", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  values <- lapply(parts, function(part) {
    key_value <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(key_value) < 2) {
      return(NULL)
    }
    key <- trimws(key_value[1])
    value <- paste(key_value[-1], collapse = "=")
    value <- tryCatch(utils::URLdecode(value), error = function(e) value)
    stats::setNames(value, key)
  })
  values <- unlist(values, use.names = TRUE)
  if (is.null(values)) {
    character()
  } else {
    values
  }
}

metadata_value <- function(metadata_df, column, n) {
  if (column %in% names(metadata_df)) {
    value <- as.character(metadata_df[[column]])
    value[is.na(value)] <- ""
    return(value)
  }
  rep("", n)
}

read_gff_fallback <- function(gff_file) {
  lines <- readLines(gff_file, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  if (length(lines) == 0) {
    return(data.frame())
  }
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 9]
  if (length(fields) == 0) {
    return(data.frame())
  }
  mat <- do.call(rbind, lapply(fields, function(x) x[seq_len(9)]))
  gff <- data.frame(
    seqid = mat[, 1],
    source = mat[, 2],
    type = mat[, 3],
    start = as.integer(mat[, 4]),
    end = as.integer(mat[, 5]),
    score = mat[, 6],
    strand = mat[, 7],
    phase = mat[, 8],
    attributes = mat[, 9],
    stringsAsFactors = FALSE
  )
  bad_ranges <- !is.na(gff$start) & !is.na(gff$end) & gff$end < gff$start
  if (any(bad_ranges)) {
    old_start <- gff$start[bad_ranges]
    gff$start[bad_ranges] <- gff$end[bad_ranges]
    gff$end[bad_ranges] <- old_start
  }
  attrs <- lapply(gff$attributes, parse_gff_attributes)
  attr_names <- unique(unlist(lapply(attrs, names), use.names = FALSE))
  for (name in attr_names) {
    gff[[name]] <- vapply(attrs, function(x) {
      if (name %in% names(x)) {
        x[[name]]
      } else {
        ""
      }
    }, character(1))
  }
  gff
}

read_gff_table <- function(gff_file) {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gr <- tryCatch(
      rtracklayer::import(gff_file, format = "gff3"),
      error = function(e) NULL
    )
    if (!is.null(gr)) {
      metadata <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
      return(data.frame(
        seqid = as.character(GenomicRanges::seqnames(gr)),
        type = as.character(metadata$type %||% metadata$source %||% ""),
        start = GenomicRanges::start(gr),
        end = GenomicRanges::end(gr),
        strand = as.character(GenomicRanges::strand(gr)),
        metadata,
        stringsAsFactors = FALSE
      ))
    }
  }
  read_gff_fallback(gff_file)
}

prepare_gff_features <- function(params, arm_chroms) {
  if (!params$annotation_enabled) {
    return(list(features = NULL, diagnostics = empty_annotation_diagnostics()))
  }
  if (!is_present(params$target_annotations_path) || !file.exists(params$target_annotations_path)) {
    return(list(
      features = NULL,
      diagnostics = empty_annotation_diagnostics(
        "skipped",
        paste("GFF file is missing or not configured:", params$target_annotations_path)
      )
    ))
  }

  gff <- read_gff_table(params$target_annotations_path)
  if (nrow(gff) == 0) {
    return(list(features = NULL, diagnostics = empty_annotation_diagnostics("skipped", "GFF contained no usable feature rows.")))
  }

  feature_type <- as.character(gff$type)
  keep <- !feature_type %in% params$annotation_ignore_feature_types
  if (length(params$annotation_feature_types) > 0) {
    keep <- keep & feature_type %in% params$annotation_feature_types
  }
  gff <- gff[keep, , drop = FALSE]
  feature_type <- feature_type[keep]
  if (nrow(gff) == 0) {
    return(list(features = NULL, diagnostics = empty_annotation_diagnostics("skipped", "No GFF rows remained after feature-type filtering.")))
  }

  seqid <- as.character(gff$seqid)
  arm_chroms <- unique(as.character(arm_chroms))
  stripped_arms <- strip_seq_version(arm_chroms)
  overlap_seqid <- ifelse(
    seqid %in% arm_chroms,
    seqid,
    ifelse(strip_seq_version(seqid) %in% stripped_arms, strip_seq_version(seqid), seqid)
  )

  feature_id <- metadata_value(gff, "ID", nrow(gff))
  feature_name <- metadata_value(gff, "Name", nrow(gff))
  feature_id[!nzchar(feature_id)] <- feature_name[!nzchar(feature_id)]
  feature_id[!nzchar(feature_id)] <- paste0("gff_feature_", seq_len(nrow(gff)))[!nzchar(feature_id)]

  features <- data.frame(
    feature_row = seq_len(nrow(gff)),
    seqid = seqid,
    overlap_seqid = overlap_seqid,
    start = as.integer(gff$start),
    end = as.integer(gff$end),
    strand = ifelse(as.character(gff$strand) %in% c("+", "-"), as.character(gff$strand), "*"),
    feature_id = feature_id,
    feature_name = feature_name,
    feature_type = feature_type,
    gene = metadata_value(gff, "gene", nrow(gff)),
    locus_tag = metadata_value(gff, "locus_tag", nrow(gff)),
    Alias = metadata_value(gff, "Alias", nrow(gff)),
    Parent = metadata_value(gff, "Parent", nrow(gff)),
    product = metadata_value(gff, "product", nrow(gff)),
    gbkey = metadata_value(gff, "gbkey", nrow(gff)),
    gene_biotype = metadata_value(gff, "gene_biotype", nrow(gff)),
    biotype = metadata_value(gff, "biotype", nrow(gff)),
    type_attribute = metadata_value(gff, "type", nrow(gff)),
    stringsAsFactors = FALSE
  )
  features$seqname_normalized <- features$seqid != features$overlap_seqid
  list(features = features, diagnostics = NULL)
}

interaction_arm_table <- function(analysis, result_classes) {
  bedpe <- analysis$bedpe
  common <- data.frame(
    interaction_id = bedpe$interaction_id,
    result_class = result_classes$result_class[match(bedpe$interaction_id, result_classes$interaction_id)],
    high_confidence = result_classes$high_confidence[match(bedpe$interaction_id, result_classes$interaction_id)],
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(
    data.frame(
      interaction_id = bedpe$interaction_id,
      arm = "arm1",
      chrom = bedpe$chr1,
      start = as.integer(bedpe$start1) + 1L,
      end = as.integer(bedpe$end1),
      strand = bedpe$strand1,
      rnanue_arm_features = bedpe$arm1_features,
      stringsAsFactors = FALSE
    ),
    data.frame(
      interaction_id = bedpe$interaction_id,
      arm = "arm2",
      chrom = bedpe$chr2,
      start = as.integer(bedpe$start2) + 1L,
      end = as.integer(bedpe$end2),
      strand = bedpe$strand2,
      rnanue_arm_features = bedpe$arm2_features,
      stringsAsFactors = FALSE
    )
  ) |>
    dplyr::left_join(common, by = "interaction_id")
}

feature_priority <- function(feature_type) {
  priorities <- c(
    sRNA = 1, tRNA = 2, rRNA = 3, ncRNA = 4,
    gene = 5, transcript = 6, CDS = 7, TU = 8,
    mobile_genetic_element = 9
  )
  unname(priorities[feature_type] %||% 99)
}

normalize_annotation_strand <- function(x) {
  x <- as.character(x)
  ifelse(x %in% c("+", "-"), x, "*")
}

strand_policy_matches <- function(arm_strand, feature_strand, policy) {
  arm_strand <- normalize_annotation_strand(arm_strand)
  feature_strand <- normalize_annotation_strand(feature_strand)
  policy <- tolower(policy %||% "both")

  if (identical(policy, "both")) {
    return(rep(TRUE, length(arm_strand)))
  }
  if (identical(policy, "same")) {
    return(arm_strand == feature_strand | arm_strand == "*" | feature_strand == "*")
  }
  if (identical(policy, "opposite")) {
    return(
      (arm_strand == "+" & feature_strand == "-") |
        (arm_strand == "-" & feature_strand == "+")
    )
  }
  stop("annotation.strand_policy must be one of same, opposite, or both.", call. = FALSE)
}

classify_partner_type <- function(feature_type, gene_biotype, biotype, gbkey, type_attribute, status) {
  if (!identical(status, "annotated")) {
    return(status)
  }
  values <- tolower(paste(feature_type, gene_biotype, biotype, gbkey, type_attribute, sep = " "))
  if (grepl("protein_coding|\\bcds\\b", values)) return("protein_coding")
  if (grepl("\\bsrna\\b|small regulatory", values)) return("sRNA")
  if (grepl("\\btrna\\b", values)) return("tRNA")
  if (grepl("\\brrna\\b", values)) return("rRNA")
  if (grepl("\\bncrna\\b|non.?coding", values)) return("ncRNA")
  if (tolower(feature_type) == "tu") return("transcription_unit")
  if (grepl("mobile_genetic_element", values)) return("mobile_genetic_element")
  "ambiguous"
}

make_arm_annotations <- function(arms, features, params) {
  if (is.null(features) || nrow(features) == 0 || nrow(arms) == 0) {
    out <- empty_arm_annotation_table()
    if (nrow(arms) > 0) {
      out <- data.frame(
        interaction_id = arms$interaction_id,
        arm = arms$arm,
        chrom = arms$chrom,
        start = arms$start,
        end = arms$end,
        strand = arms$strand,
        rnanue_arm_features = arms$rnanue_arm_features,
        selected_feature_id = "",
        selected_feature_name = "",
        selected_feature_type = "",
        selected_partner_class = "unannotated",
        feature_annotation_source = "unannotated_no_gff_overlap",
        partner_class_source = "unannotated_no_gff_overlap",
        gene = "",
        locus_tag = "",
        Alias = "",
        Parent = "",
        product = "",
        gbkey = "",
        gene_biotype = "",
        overlap_bp = 0L,
        arm_overlap_fraction = 0,
        feature_overlap_fraction = 0,
        reciprocal_overlap_fraction = 0,
        annotation_status = "unannotated",
        seqname_normalized = FALSE,
        result_class = arms$result_class,
        high_confidence = arms$high_confidence,
        stringsAsFactors = FALSE
      )
    }
    return(out)
  }

  arm_seqid <- if (params$annotation_seqname_normalization == "auto_strip_version") {
    feature_seqids <- unique(features$overlap_seqid)
    ifelse(arms$chrom %in% feature_seqids, arms$chrom, strip_seq_version(arms$chrom))
  } else {
    arms$chrom
  }
  arm_strand <- normalize_annotation_strand(arms$strand)
  feature_strand <- normalize_annotation_strand(features$strand)
  arm_gr <- GenomicRanges::GRanges(
    seqnames = arm_seqid,
    ranges = IRanges::IRanges(start = arms$start, end = arms$end),
    strand = arm_strand
  )
  feature_gr <- GenomicRanges::GRanges(
    seqnames = features$overlap_seqid,
    ranges = IRanges::IRanges(start = features$start, end = features$end),
    strand = feature_strand
  )
  hits <- GenomicRanges::findOverlaps(
    arm_gr,
    feature_gr,
    ignore.strand = TRUE
  )
  if (length(hits) == 0) {
    return(make_arm_annotations(arms, NULL, params))
  }

  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)
  strand_keep <- strand_policy_matches(
    arm_strand[q],
    feature_strand[s],
    params$annotation_strand_policy %||% if (isTRUE(params$annotation_same_strand)) "same" else "both"
  )
  q <- q[strand_keep]
  s <- s[strand_keep]
  if (length(q) == 0) {
    return(make_arm_annotations(arms, NULL, params))
  }
  overlap_start <- pmax(arms$start[q], features$start[s])
  overlap_end <- pmin(arms$end[q], features$end[s])
  overlap_bp <- pmax(0L, overlap_end - overlap_start + 1L)
  candidates <- data.frame(
    arm_row = q,
    feature_row = s,
    overlap_bp = overlap_bp,
    arm_overlap_fraction = overlap_bp / pmax(1L, arms$end[q] - arms$start[q] + 1L),
    feature_overlap_fraction = overlap_bp / pmax(1L, features$end[s] - features$start[s] + 1L),
    feature_priority = vapply(features$feature_type[s], feature_priority, numeric(1)),
    stringsAsFactors = FALSE
  )
  candidates$reciprocal_overlap_fraction <- pmin(
    candidates$arm_overlap_fraction,
    candidates$feature_overlap_fraction
  )
  candidates <- candidates[candidates$overlap_bp >= params$annotation_min_overlap_bp, , drop = FALSE]
  if (nrow(candidates) == 0) {
    return(make_arm_annotations(arms, NULL, params))
  }
  candidates <- candidates[order(
    candidates$arm_row,
    -candidates$reciprocal_overlap_fraction,
    -candidates$arm_overlap_fraction,
    -candidates$overlap_bp,
    candidates$feature_priority
  ), , drop = FALSE]
  best <- candidates[!duplicated(candidates$arm_row), , drop = FALSE]

  rows <- lapply(seq_len(nrow(arms)), function(i) {
    hit <- best[best$arm_row == i, , drop = FALSE]
    if (nrow(hit) == 0) {
      selected <- NULL
      status <- "unannotated"
    } else {
      selected <- features[hit$feature_row[1], , drop = FALSE]
      status <- "annotated"
    }
    partner_class <- if (is.null(selected)) {
      "unannotated"
    } else {
      classify_partner_type(
        selected$feature_type,
        selected$gene_biotype,
        selected$biotype,
        selected$gbkey,
        selected$type_attribute,
        status
      )
    }
    source <- if (identical(status, "annotated")) {
      "gff_overlap"
    } else {
      "unannotated_no_gff_overlap"
    }
    data.frame(
      interaction_id = arms$interaction_id[i],
      arm = arms$arm[i],
      chrom = arms$chrom[i],
      start = arms$start[i] - 1L,
      end = arms$end[i],
      strand = arms$strand[i],
      rnanue_arm_features = arms$rnanue_arm_features[i],
      selected_feature_id = selected$feature_id %||% "",
      selected_feature_name = selected$feature_name %||% "",
      selected_feature_type = selected$feature_type %||% "",
      selected_partner_class = partner_class,
      feature_annotation_source = source,
      partner_class_source = source,
      gene = selected$gene %||% "",
      locus_tag = selected$locus_tag %||% "",
      Alias = selected$Alias %||% "",
      Parent = selected$Parent %||% "",
      product = selected$product %||% "",
      gbkey = selected$gbkey %||% "",
      gene_biotype = selected$gene_biotype %||% "",
      overlap_bp = hit$overlap_bp[1] %||% 0L,
      arm_overlap_fraction = hit$arm_overlap_fraction[1] %||% 0,
      feature_overlap_fraction = hit$feature_overlap_fraction[1] %||% 0,
      reciprocal_overlap_fraction = hit$reciprocal_overlap_fraction[1] %||% 0,
      annotation_status = status,
      seqname_normalized = !identical(arms$chrom[i], as.character(GenomicRanges::seqnames(arm_gr))[i]) ||
        (selected$seqname_normalized %||% FALSE),
      result_class = arms$result_class[i],
      high_confidence = arms$high_confidence[i],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

partner_class_order <- function() {
  c(
    "sRNA", "protein_coding", "tRNA", "rRNA", "ncRNA",
    "transcription_unit", "mobile_genetic_element", "ambiguous", "unannotated"
  )
}

ordered_partner_pair <- function(class1, class2) {
  levels <- partner_class_order()
  rank1 <- match(class1, levels)
  rank2 <- match(class2, levels)
  rank1[is.na(rank1)] <- length(levels) + 1
  rank2[is.na(rank2)] <- length(levels) + 1
  swap <- rank2 < rank1
  data.frame(
    partner_class_1 = ifelse(swap, class2, class1),
    partner_class_2 = ifelse(swap, class1, class2),
    stringsAsFactors = FALSE
  )
}

make_partner_type_pairs <- function(arm_annotations) {
  if (nrow(arm_annotations) == 0) {
    return(empty_partner_type_pairs())
  }
  arm1 <- arm_annotations[arm_annotations$arm == "arm1", , drop = FALSE]
  arm2 <- arm_annotations[arm_annotations$arm == "arm2", , drop = FALSE]
  pairs <- dplyr::inner_join(
    arm1[, c("interaction_id", "selected_partner_class", "result_class", "high_confidence"), drop = FALSE],
    arm2[, c("interaction_id", "selected_partner_class"), drop = FALSE],
    by = "interaction_id",
    suffix = c("_arm1", "_arm2")
  )
  ordered <- ordered_partner_pair(pairs$selected_partner_class_arm1, pairs$selected_partner_class_arm2)
  data.frame(
    interaction_id = pairs$interaction_id,
    arm1_partner_class = pairs$selected_partner_class_arm1,
    arm2_partner_class = pairs$selected_partner_class_arm2,
    partner_class_1 = ordered$partner_class_1,
    partner_class_2 = ordered$partner_class_2,
    partner_class_pair = paste(ordered$partner_class_1, ordered$partner_class_2, sep = "__"),
    result_class = pairs$result_class,
    high_confidence = pairs$high_confidence,
    stringsAsFactors = FALSE
  )
}

ids_for_result_set <- function(analysis, set_name) {
  switch(
    set_name,
    high_confidence = if (!is.null(analysis$high_confidence_results)) analysis$high_confidence_results$interaction_id else character(),
    concordant_same_direction = if (!is.null(analysis$concordance)) {
      analysis$concordance$interaction_id[analysis$concordance$concordance_class == "both_significant_same_direction"]
    } else character(),
    edgeR_significant = if ("edgeR" %in% names(analysis$method_results)) {
      res <- analysis$method_results$edgeR$annotated_results
      res$interaction_id[res$significant]
    } else character(),
    DESeq2_significant = if ("DESeq2" %in% names(analysis$method_results)) {
      res <- analysis$method_results$DESeq2$annotated_results
      res$interaction_id[res$significant]
    } else character(),
    called = called_interaction_ids(analysis),
    character()
  )
}

make_partner_type_confusion <- function(analysis, partner_pairs) {
  if (nrow(partner_pairs) == 0) {
    return(empty_partner_type_confusion())
  }
  sets <- unique(c("high_confidence", "concordant_same_direction", "edgeR_significant", "DESeq2_significant", "called"))
  rows <- lapply(sets, function(set_name) {
    ids <- ids_for_result_set(analysis, set_name)
    set_pairs <- partner_pairs[partner_pairs$interaction_id %in% ids, , drop = FALSE]
    if (nrow(set_pairs) == 0) {
      return(NULL)
    }
    counts <- as.data.frame(table(set_pairs$partner_class_1, set_pairs$partner_class_2), stringsAsFactors = FALSE)
    names(counts) <- c("partner_class_1", "partner_class_2", "interactions")
    counts <- counts[counts$interactions > 0, , drop = FALSE]
    counts$result_set <- set_name
    counts[, c("result_set", "partner_class_1", "partner_class_2", "interactions"), drop = FALSE]
  })
  dplyr::bind_rows(rows)
}

matchable_annotation_values <- function(row) {
  values <- c(
    row$selected_feature_id,
    row$selected_feature_name,
    row$gene,
    row$locus_tag,
    row$Alias,
    row$Parent,
    parse_feature_field(row$rnanue_arm_features)
  )
  values <- unlist(strsplit(values, ",", fixed = TRUE), use.names = FALSE)
  values <- trimws(values)
  unique(values[nzchar(values)])
}

filter_significant_annotation_arms <- function(analysis, arms) {
  called_ids <- called_interaction_ids(analysis)
  called_ids <- called_ids[!is.na(called_ids) & nzchar(called_ids)]
  if (length(called_ids) == 0 || nrow(arms) == 0) {
    return(arms[0, , drop = FALSE])
  }
  arms[arms$interaction_id %in% called_ids, , drop = FALSE]
}

make_target_feature_matches <- function(arm_annotations, target_ids) {
  if (length(target_ids) == 0 || nrow(arm_annotations) == 0) {
    return(empty_target_feature_matches())
  }
  rows <- list()
  for (i in seq_len(nrow(arm_annotations))) {
    values <- matchable_annotation_values(arm_annotations[i, , drop = FALSE])
    matched <- intersect(target_ids, values)
    if (length(matched) == 0) next
    partner <- arm_annotations[
      arm_annotations$interaction_id == arm_annotations$interaction_id[i] &
        arm_annotations$arm != arm_annotations$arm[i],
      ,
      drop = FALSE
    ]
    rows[[length(rows) + 1]] <- data.frame(
      target_id = matched,
      interaction_id = arm_annotations$interaction_id[i],
      target_arm = arm_annotations$arm[i],
      target_partner_class = arm_annotations$selected_partner_class[i],
      other_arm_partner_class = partner$selected_partner_class[1] %||% "",
      selected_feature_id = arm_annotations$selected_feature_id[i],
      selected_feature_name = arm_annotations$selected_feature_name[i],
      selected_feature_type = arm_annotations$selected_feature_type[i],
      gene = arm_annotations$gene[i],
      locus_tag = arm_annotations$locus_tag[i],
      Alias = arm_annotations$Alias[i],
      Parent = arm_annotations$Parent[i],
      product = arm_annotations$product[i],
      gbkey = arm_annotations$gbkey[i],
      gene_biotype = arm_annotations$gene_biotype[i],
      rnanue_arm_features = arm_annotations$rnanue_arm_features[i],
      result_class = arm_annotations$result_class[i],
      high_confidence = arm_annotations$high_confidence[i],
      overlap_bp = arm_annotations$overlap_bp[i],
      reciprocal_overlap_fraction = arm_annotations$reciprocal_overlap_fraction[i],
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(empty_target_feature_matches())
  }
  out[order(out$target_id, out$interaction_id, out$target_arm), , drop = FALSE]
}

make_target_gene_interactions <- function(analysis, target_feature_matches) {
  if (nrow(target_feature_matches) == 0) {
    return(empty_target_gene_interactions())
  }
  method <- primary_result_method(analysis)
  primary_results <- analysis$method_results[[method]]$annotated_results
  stats <- primary_results[, intersect(c("interaction_id", "log2FoldChange", "padj"), names(primary_results)), drop = FALSE]
  out <- dplyr::left_join(
    target_feature_matches[, c(
      "target_id", "interaction_id", "target_arm", "target_partner_class",
      "other_arm_partner_class", "selected_feature_id", "selected_feature_name",
      "rnanue_arm_features", "result_class", "high_confidence"
    ), drop = FALSE],
    stats,
    by = "interaction_id"
  )
  out[order(out$target_id, out$padj, out$interaction_id), , drop = FALSE]
}

make_target_gene_bedpe <- function(analysis, target_interactions) {
  if (nrow(target_interactions) == 0) {
    return(empty_target_gene_bedpe())
  }
  grouped <- split(target_interactions, target_interactions$interaction_id)
  rows <- lapply(names(grouped), function(interaction_id) {
    group <- grouped[[interaction_id]]
    bedpe <- analysis$bedpe[interaction_id, , drop = FALSE]
    data.frame(
      chr1 = bedpe$chr1,
      start1 = bedpe$start1,
      end1 = bedpe$end1,
      chr2 = bedpe$chr2,
      start2 = bedpe$start2,
      end2 = bedpe$end2,
      interaction_id = interaction_id,
      target_ids = paste(unique(group$target_id), collapse = ";"),
      target_arms = paste(unique(group$target_arm), collapse = ";"),
      result_class = group$result_class[1],
      high_confidence = group$high_confidence[1],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

make_annotation_diagnostics <- function(params, features, arm_annotations, target_ids, target_interactions, gff_diag) {
  unmatched_targets <- setdiff(target_ids, target_interactions$target_id)
  data.frame(
    metric = c(
      "status",
      "message",
      "interaction_scope",
      "strand_policy",
      "gff_features_used",
      "annotated_arms",
      "unannotated_arms",
      "ambiguous_arms",
      "seqname_normalized_arms",
      "target_ids_requested",
      "target_ids_matched",
      "unmatched_target_ids"
    ),
    value = c(
      "computed",
      gff_diag$value[match("message", gff_diag$metric)] %||% "GFF annotation computed.",
      "significant_interactions",
      params$annotation_strand_policy %||% if (isTRUE(params$annotation_same_strand)) "same" else "both",
      if (is.null(features)) 0 else nrow(features),
      sum(arm_annotations$annotation_status == "annotated", na.rm = TRUE),
      sum(arm_annotations$annotation_status == "unannotated", na.rm = TRUE),
      sum(arm_annotations$selected_partner_class == "ambiguous", na.rm = TRUE),
      sum(arm_annotations$seqname_normalized, na.rm = TRUE),
      length(target_ids),
      length(intersect(target_ids, target_interactions$target_id)),
      paste(unmatched_targets, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
}

make_gff_annotation_assessments <- function(analysis, result_classes) {
  if (!analysis$params$annotation_enabled) {
    return(empty_gff_annotation_assessment())
  }
  arms <- filter_significant_annotation_arms(
    analysis,
    interaction_arm_table(analysis, result_classes)
  )
  target_ids <- split_param_ids(analysis$params$annotation_target_ids)
  if (nrow(arms) == 0) {
    return(list(
      arm_annotations = empty_arm_annotation_table(),
      partner_type_pairs = empty_partner_type_pairs(),
      partner_type_confusion = empty_partner_type_confusion(),
      target_gene_interactions = empty_target_gene_interactions(),
      target_gene_bedpe = empty_target_gene_bedpe(),
      diagnostics = make_annotation_diagnostics(
        analysis$params,
        NULL,
        empty_arm_annotation_table(),
        target_ids,
        empty_target_gene_interactions(),
        empty_annotation_diagnostics("computed", "No significant interactions available for GFF annotation.")
      )
    ))
  }
  prepared <- prepare_gff_features(analysis$params, arms$chrom)
  if (is.null(prepared$features)) {
    assessment <- empty_gff_annotation_assessment(
      prepared$diagnostics$value[prepared$diagnostics$metric == "status"][1],
      prepared$diagnostics$value[prepared$diagnostics$metric == "message"][1]
    )
    assessment$arm_annotations <- make_arm_annotations(arms, NULL, analysis$params)
    return(assessment)
  }
  arm_annotations <- make_arm_annotations(arms, prepared$features, analysis$params)
  partner_pairs <- make_partner_type_pairs(arm_annotations)
  partner_confusion <- make_partner_type_confusion(analysis, partner_pairs)
  target_feature_matches <- make_target_feature_matches(arm_annotations, target_ids)
  target_interactions <- make_target_gene_interactions(analysis, target_feature_matches)
  target_bedpe <- make_target_gene_bedpe(analysis, target_interactions)
  list(
    arm_annotations = arm_annotations,
    partner_type_pairs = partner_pairs,
    partner_type_confusion = partner_confusion,
    target_gene_interactions = target_interactions,
    target_gene_bedpe = target_bedpe,
    diagnostics = make_annotation_diagnostics(
      analysis$params,
      prepared$features,
      arm_annotations,
      target_ids,
      target_interactions,
      prepared$diagnostics %||% empty_annotation_diagnostics("computed", "GFF annotation computed.")
    )
  )
}
