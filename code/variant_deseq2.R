run_method_DESeq2 <- function(analysis) {
  raw_counts <- as.matrix(analysis$counts_model)
  rounded_counts <- round(raw_counts)
  storage.mode(rounded_counts) <- "integer"
  keep <- custom_count_filter(rounded_counts, analysis$params)
  if (!any(keep)) {
    stop("DESeq2 filtering removed all interactions.", call. = FALSE)
  }

  col_data <- analysis$model_metadata
  rownames(col_data) <- col_data$sample_id
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = rounded_counts[keep, , drop = FALSE],
    colData = col_data,
    design = analysis$design_formula
  )

  nf <- exp(analysis$log_offset[rownames(dds), colnames(dds), drop = FALSE])
  nf <- nf / exp(rowMeans(log(nf)))
  nf[!is.finite(nf) | nf <= 0] <- 1
  DESeq2::normalizationFactors(dds) <- nf

  dds <- tryCatch(
    DESeq2::DESeq(dds, fitType = "local", quiet = TRUE),
    error = function(local_error) {
      warning("DESeq2 fitType='local' failed; retrying with fitType='mean': ", local_error$message, call. = FALSE)
      tryCatch(
        DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE),
        error = function(mean_error) {
          warning(
            "DESeq2 fitType='mean' failed; using gene-wise dispersion estimates and Wald test: ",
            mean_error$message,
            call. = FALSE
          )
          dds_fallback <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
          DESeq2::dispersions(dds_fallback) <- S4Vectors::mcols(dds_fallback)$dispGeneEst
          DESeq2::nbinomWaldTest(dds_fallback, quiet = TRUE)
        }
      )
    }
  )
  raw_res <- DESeq2::results(
    dds,
    contrast = c("contrast_group", analysis$case_condition, analysis$control_condition),
    alpha = analysis$params$padj_threshold
  )

  coef_candidates <- grep("^contrast_group_.*_vs_.*$", DESeq2::resultsNames(dds), value = TRUE)
  shrink_res <- NULL
  if (length(coef_candidates) == 1 && requireNamespace("apeglm", quietly = TRUE)) {
    shrink_res <- tryCatch(
      DESeq2::lfcShrink(dds, coef = coef_candidates, type = "apeglm", quiet = TRUE),
      error = function(e) {
        warning("DESeq2 lfcShrink(apeglm) failed; keeping unshrunken LFC: ", e$message, call. = FALSE)
        NULL
      }
    )
  }

  result_df <- as.data.frame(raw_res, stringsAsFactors = FALSE)
  result_df$interaction_id <- rownames(result_df)
  result_df$log2FoldChange_raw <- result_df$log2FoldChange
  result_df$log2FoldChange_shrunken <- NA_real_
  result_df$lfcSE_shrunken <- NA_real_
  result_df$lfc_shrinkage_ratio <- NA_real_
  if (!is.null(shrink_res)) {
    shrink_df <- as.data.frame(shrink_res, stringsAsFactors = FALSE)
    result_df$log2FoldChange_shrunken <- shrink_df[result_df$interaction_id, "log2FoldChange"]
    result_df$lfcSE_shrunken <- shrink_df[result_df$interaction_id, "lfcSE"]
    result_df$lfc_shrinkage_ratio <- abs(result_df$log2FoldChange_shrunken) /
      pmax(abs(result_df$log2FoldChange_raw), .Machine$double.xmin)
  }
  result_df$method <- "DESeq2"
  result_df <- result_df[order(result_df$padj, result_df$pvalue), , drop = FALSE]

  list(
    method = "DESeq2",
    counts = rounded_counts[keep, , drop = FALSE],
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    filter_keep = keep,
    dataset = dds,
    raw_results = raw_res,
    shrink_results = shrink_res,
    results = result_df,
    rounding = rounding_diagnostics(raw_counts, rounded_counts)
  )
}
