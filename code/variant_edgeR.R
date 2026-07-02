run_method_edgeR <- function(analysis) {
  counts <- as.matrix(analysis$counts_model)
  storage.mode(counts) <- "numeric"

  dge <- edgeR::DGEList(counts = counts, group = analysis$model_metadata$contrast_group)
  keep_custom <- custom_count_filter(dge$counts, analysis$params)
  keep_expr <- tryCatch(
    edgeR::filterByExpr(dge, design = analysis$design),
    error = function(e) {
      warning("edgeR::filterByExpr failed; using only custom count filters: ", e$message, call. = FALSE)
      rep(TRUE, nrow(dge))
    }
  )
  keep <- keep_custom & keep_expr
  filtered_counts <- dge$counts[keep, , drop = FALSE]
  filtered_offset <- analysis$log_offset[rownames(filtered_counts), colnames(filtered_counts), drop = FALSE]
  skip_reason <- method_fit_skip_reason(filtered_counts, filtered_offset)
  if (!is.null(skip_reason)) {
    warning("edgeR skipped: ", skip_reason, call. = FALSE)
    return(empty_method_result("edgeR", dge$counts, keep, reason = skip_reason))
  }

  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR_offset <- analysis$log_offset[rownames(dge$counts), colnames(dge$counts), drop = FALSE]
  edgeR_offset <- edgeR_offset +
    matrix(log(dge$samples$norm.factors), nrow = nrow(edgeR_offset), ncol = ncol(edgeR_offset), byrow = TRUE)
  dge$offset <- edgeR_offset

  dge <- edgeR::estimateDisp(dge, analysis$design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, analysis$design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = analysis$contrast_coef)

  results <- as.data.frame(qlf$table, stringsAsFactors = FALSE)
  results$interaction_id <- rownames(results)
  results$FDR <- stats::p.adjust(results$PValue, method = "BH")
  results$method <- "edgeR"
  results$log2FoldChange <- results$logFC
  results$pvalue <- results$PValue
  results$padj <- results$FDR
  results <- results[order(results$padj, results$pvalue), , drop = FALSE]

  list(
    method = "edgeR",
    counts = dge$counts,
    normalized_counts = offset_normalized_counts(dge$counts, edgeR_offset),
    filter_keep = keep,
    design = analysis$design,
    dataset = dge,
    fit = fit,
    test = qlf,
    results = results,
    rounding = NULL
  )
}
