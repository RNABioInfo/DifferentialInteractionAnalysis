annotate_interactions = function(interactions_df, annotation_gff) {
  sig_gene_ids_1 = get_gene_ids(interactions_df$X12)
  
  gene_names_1 = lapply(sig_gene_ids_1, function(elem) {
    return(get_info_with_ids_from_gff("Name", elem, annotation_gff))
  })
  gene_biotype_1 = lapply(sig_gene_ids_1, function(elem) {
    return(get_info_with_ids_from_gff("gene_biotype", elem, annotation_gff, "gbkey"))
  })
  
  sig_gene_ids_2 = get_gene_ids(interactions_df$X13)
  
  gene_names_2 = lapply(sig_gene_ids_2, function(elem) {
    return(get_info_with_ids_from_gff("Name", elem, annotation_gff))
  })
  gene_biotype_2 = lapply(sig_gene_ids_2, function(elem) {
    return(get_info_with_ids_from_gff("gene_biotype", elem, annotation_gff, "gbkey"))
  })
  
  sig_annotation = data.frame(unlist(gene_names_1), unlist(gene_biotype_1), unlist(gene_names_2), unlist(gene_biotype_2))
  
  colnames(sig_annotation) = c("gene_name_1", "gene_biotype_1", "gene_name_2", "gene_biotype_2")
  
  results_df = interactions_df[,1:10] 
  
  results_df = cbind(results_df, sig_annotation)
  
  return(results_df)
}