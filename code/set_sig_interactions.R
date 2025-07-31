set_sig_interactions = function(obj, interactions_bed) {
  UseMethod("set_sig_interactions")
}

set_sig_interactions.default = function(obj, interactions_bed) {
  
  sig_results = sig_results(obj)
  sig_interactions = interactions_bed[rownames(sig_results_df), ]
  
  obj$sig_interaction_bed = sig_interactions
  return(obj)
}