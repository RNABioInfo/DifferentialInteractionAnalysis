make_analysis <- function(params) {
  metadata_df = read_metadata(params$metadata_file)
  
  print(metadata_df)
  
  if (params$analysis_method == "edgeR") {
    analysis = EdgeRAnalysis(params, metadata_df)
  } else if (params$analysis_method == "DESeq2") {
    analysis = DESeq2Analysis(params, metadata_df)
  } else {
    stop(str_glue("Analysis methods did not match anly of the options: {params$de_variant}"))
  }
  
  return(analysis)
}


get_info_with_ids_from_gff <- function(info, ids, gff_df, alt_info = "") {
  # sanity‐check
  if (!("X9" %in% names(gff_df))) {
    stop("gff_df must have a column named X9")
  }
  
  # for each gene_id, pull out its hits, extract `info=…`, then collapse
  results <- vapply(ids, FUN.VALUE = character(1), USE.NAMES = FALSE, function(gene_id) {
    # build ID=…; pattern
    pat_id <- str_glue("ID={gene_id};")
    hits   <- grep(pat_id, gff_df$X9)
    
    if (length(hits) == 0) {
      return("*")   # no hit → empty string
    }
    if (length(hits) > 1) {
      message(str_glue("Found more than 1 hit for {gene_id}"))
    }
    
    pattern = str_glue(".*{info}=([^;]+)(?:;|$).*")
    alt_pattern = str_glue(".*{alt_info}=([^;]+)(?:;|$).*")
    
    # extract the info field from each matched row
    info_vals <- ifelse(
      grepl(pattern, gff_df$X9[hits]),
      sub(pattern,
          "\\1",
          gff_df$X9[hits]
      ),
      sub(alt_pattern,
          "\\1",
          gff_df$X9[hits]
      )
    )
      
      
    # collapse into one comma‐separated string
    return(paste(info_vals, collapse = ","))
  })
  
  # returns a plain character vector, same length/order as `ids`
  return(paste(results, collapse = ","))
}

get_gene_ids <- function(col) {
  gene_ids = str_split(col, ";")
  gene_ids = lapply(gene_ids, function(x) x[-length(x)])
  gene_ids = lapply(gene_ids, function(vec) {
    sapply(strsplit(vec, ","), `[`, 1)
  })
  gene_ids = lapply(gene_ids, function(vec) {
    sub("^[^:]*:(.*)$", "\\1", vec)
  })
  return(lapply(gene_ids, unique))
}

lists_to_df <- function(...) {
  # capture the call and the raw args
  mc   <- match.call(expand.dots = FALSE)$...
  args <- list(...)
  
  # determine a name for each arg: 
  #  • if the user wrote foo=…, use “foo” 
  #  • otherwise, deparse the symbol they passed
  arg_names <- names(args)
  for (i in seq_along(arg_names)) {
    if (is.null(arg_names) || arg_names[i] == "") {
      arg_names[i] <- deparse(mc[[i]])
    }
  }
  names(args) <- arg_names
  
  # now extract each list-of-length-1 into a bare vector,
  # checking that it really is a list of length 1 of character:
  cols <- lapply(args, function(lst) {
    if (!is.list(lst) || length(lst) != 1 || !is.character(lst[[1]])) {
      stop("Each argument must be a list of length 1 containing a character vector")
    }
    lst[[1]]
  })
  
  # finally bind into a data.frame, no factors, no name mangling
  as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}

# 1. Prepare and filter the data
prepare_unambiguous_results = function(df) {
  # keep only rows without commas or asterisks in biotypes
  good = subset(
    df,
    !grepl(",|[*]", gene_biotype_1) &
      !grepl(",|[*]", gene_biotype_2)
  )
  
  # add a factor column for intra vs inter
  good$interaction_type = factor(
    ifelse(good$gene_name_1 == good$gene_name_2, "intra", "inter"),
    levels = c("intra", "inter")
  )
  
  return(good)
}

# 2. Build a biotype count matrix for a given interaction type
get_biotype_matrix = function(df, interaction) {
  subset_df = subset(df, interaction_type == interaction)
  count_mat = table(
    subset_df$gene_biotype_1,
    subset_df$gene_biotype_2
  )
  return(as.matrix(count_mat))
}

plot_biotype_heatmap = function(mat, main_title) {
  # ensure enough margin for axis labels
  setHook(
    "grid.newpage",
    function() pushViewport(
      viewport(
        x = 1, y = 1,
        width = 0.9, height = 0.9,
        name = "vp",
        just = c("right","top")
      )
    ),
    action = "prepend"
  )
  
  pheatmap(
    mat,
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    display_numbers  = TRUE,
    fontsize_number  = 10,
    angle_col        = 45,
    color            = colorRampPalette(c("white", "steelblue"))(50),
    main             = main_title
  )
  
  # reset hook and add axis labels
  setHook("grid.newpage", NULL, "replace")
  grid.text("Biotype 2", x = 0.4, y = -0.05, gp = gpar(fontsize = 12))
  grid.text("Biotype 1", x = -0.02, rot = 90,       gp = gpar(fontsize = 12))
}