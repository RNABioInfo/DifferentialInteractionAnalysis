run_analysis <- function(obj, ...) UseMethod("run_analysis")

run_analysis.InteractionAnalysis <- function(obj, ...) {
  methods <- requested_methods(obj$params$analysis_method)
  obj$analysis_warnings <- empty_analysis_warnings()

  for (method in methods) {
    captured <- capture_analysis_warnings(
      switch(
        method,
        edgeR = run_method_edgeR(obj),
        DESeq2 = run_method_DESeq2(obj)
      ),
      stage = paste0("run_method_", method),
      method = method
    )
    result <- captured$value
    obj$analysis_warnings <- append_analysis_warnings(obj$analysis_warnings, captured$warnings)
    result$annotated_results <- augment_method_results(obj, result)
    obj$method_results[[method]] <- result
  }

  captured_concordance <- capture_analysis_warnings(
    make_concordance(
      obj$method_results,
      obj$params$padj_threshold,
      obj$native_padj_by_role
    ),
    stage = "make_concordance"
  )
  obj$concordance <- captured_concordance$value
  obj$analysis_warnings <- append_analysis_warnings(obj$analysis_warnings, captured_concordance$warnings)

  captured_high_confidence <- capture_analysis_warnings(
    make_high_confidence_results(obj),
    stage = "make_high_confidence_results"
  )
  obj$high_confidence_results <- captured_high_confidence$value
  obj$analysis_warnings <- append_analysis_warnings(
    obj$analysis_warnings,
    captured_high_confidence$warnings
  )

  captured_postprocess <- capture_analysis_warnings(
    make_postprocess_assessments(obj),
    stage = "make_postprocess_assessments"
  )
  obj$postprocess <- captured_postprocess$value
  obj$analysis_warnings <- append_analysis_warnings(obj$analysis_warnings, captured_postprocess$warnings)

  write_interaction_outputs(obj)
  obj
}
