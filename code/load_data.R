read_counts.default <- function(obj) {
  counts = read_delim(obj$params$interaction_counts_path, 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE, skip = 2, name_repair = "check_unique")
  
  sample_subset = levels(obj$metadata_df$name)
  
  if(missing(sample_subset)) {
    counts$Description = NULL
    counts_df = data.frame(counts, row.names = "Name")
  } else {
    filter = c("Name", sample_subset)
    counts_df = data.frame(subset(counts, select = filter), row.names = "Name")
  }
  
  colnames(counts_df) = levels(obj$metadata_df$name)
  return(counts_df)
}

read_metadata <- function(metadata_file_path) {
  metadata_df = as.data.frame(read_delim(metadata_file_path, 
                                         delim = ",", escape_double = FALSE, 
                                         trim_ws = TRUE, col_types = cols(.default = col_factor())), stringsAsFactors = TRUE)
  rownames(metadata_df) = metadata_df$name
  return(metadata_df)
}