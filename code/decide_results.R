decide_results = function(obj, ...) {
  UseMethod("decide_results")
}

decide_results.EdgeRAnalysis = function(obj) {
  condition_levels <- levels(obj$metadata_df$condition)
  contrast_str <- paste0(condition_levels[2], " - ", condition_levels[1])
  
  # pass the string directly here:
  contrasts <- limma::makeContrasts(contrasts = contrast_str, levels = obj$edgeR_dataset$design)
  
  obj$edgeR_results = glmQLFTest(obj$edgeR_fit, contrast = contrasts)
  
  obj$results = as.data.frame(topTags(obj$edgeR_results, n = Inf, p.value = 1))
  
  return(obj)
}

decide_results.DESeq2Analysis = function(obj) {
  condition_levels = levels(obj$metadata_df$condition)
  coefficient_string = str_glue("condition_{condition_levels[2]}_vs_{condition_levels[1]}")
  res_shrink = lfcShrink(obj$deseq_dataset, coef=coefficient_string, type="apeglm")
  
  obj$deseq_results = res_shrink
  obj$results = as.data.frame(obj$deseq_results)
  
  return(obj)
}