sample_replicate_str = function(obj) {
  UseMethod("sample_replicate_str")
}

sample_replicate_str.default = function(obj) {
  sample_count = summary(obj$metadata_df$condition)
  sample_names = names(sample_count)
  
  sample_string = str_glue("{sample_names[1]} N: {as.character(sample_count[1])}; {sample_names[2]} N: {as.character(sample_count[2])}")

  return(sample_string)
}