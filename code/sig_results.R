sig_results <- function(obj, ...) UseMethod("sig_results")

sig_results.InteractionAnalysis <- function(obj, method = NULL, ...) {
  method <- method %||% names(obj$method_results)[1]
  if (is.null(method) || !method %in% names(obj$method_results)) {
    stop("Requested method has not been run: ", method, call. = FALSE)
  }
  res <- obj$method_results[[method]]$annotated_results
  res[res$significant, , drop = FALSE]
}
