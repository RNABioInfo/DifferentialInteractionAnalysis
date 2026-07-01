summary.InteractionAnalysis <- function(object, ...) {
  method_counts <- data.frame(
    method = names(object$method_results),
    tested_interactions = vapply(object$method_results, function(x) nrow(x$annotated_results), integer(1)),
    significant_interactions = vapply(
      object$method_results,
      function(x) sum(x$annotated_results$significant, na.rm = TRUE),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  print(method_counts)
  invisible(method_counts)
}
