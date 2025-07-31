run_analysis <- function(obj, ...) UseMethod("run_analysis")

run_analysis.DGEList <- function(obj, metadata_df, ...) {
  condition_levels <- levels(metadata_df$condition)
  contrast_str <- paste0(condition_levels[2], " - ", condition_levels[1])
  
  # pass the string directly here:
  contrasts <- limma::makeContrasts(contrasts = contrast_str, levels = obj$design)
  
  res <- glmQLFTest(fit, contrast = contrasts)
  
  
  return(res)
}

run_analysis.DESeqDataSet <- function(data_set) {
  results = DESeq(data_set, fitType = "local", sfType = "iterate")
  
  return(results)
}

run_analysis.EdgeRAnalysis = function(obj) {
  obj = preprocess_dataset(obj)
  obj = fit_dataset(obj)
  obj = decide_results(obj)
  
  obj$norm_counts = cpm(obj$edgeR_dataset, normalized.lib.sizes = TRUE)
  
  return(obj)
}

run_analysis.DESeq2Analysis = function(obj) {
  obj$deseq_dataset = run_analysis.DESeqDataSet(obj$deseq_dataset)
  
  obj = decide_results(obj)
  
  obj$norm_counts = counts(obj$deseq_dataset, normalized = TRUE)
  
  return(obj)
}

### Helper methods edgeR
preprocess_dataset <- function(edgeR_analysis) {
  design = model.matrix(~0+edgeR_analysis$metadata_df$condition)
  colnames(design) = levels(edgeR_analysis$metadata_df$condition)
  
  edgeR_analysis$edgeR_dataset = calcNormFactors(edgeR_analysis$edgeR_dataset)
  edgeR_analysis$edgeR_dataset = estimateDisp(edgeR_analysis$edgeR_dataset, design, robust=TRUE)
  
  return(edgeR_analysis)
}

fit_dataset = function(edgeR_analysis) {
  edgeR_analysis$edgeR_fit = glmQLFit(edgeR_analysis$edgeR_dataset, edgeR_analysis$edgeR_dataset$design, robust=TRUE)
  
  return(edgeR_analysis)
}

### Generic helpers
