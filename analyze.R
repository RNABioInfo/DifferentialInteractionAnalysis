source("load_libraries.R")

read_counts <- function(counts_file_path, sample_subset) {
  counts = read_delim(counts_file_path, 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE, skip = 2)
  
  if(missing(sample_subset)) {
    counts$Description = NULL
    return(data.frame(counts, row.names = "Name"))
  } else {
    filter = c("Name", sample_subset)
    
    return(data.frame(subset(counts, select = filter), row.names = "Name"))
  }
}

read_metadata <- function(metadata_file_path) {
  metadata_df = as.data.frame(read_delim(metadata_file_path, 
                        delim = ",", escape_double = FALSE, 
                        trim_ws = TRUE, col_types = cols(.default = col_factor())), stringsAsFactors = TRUE)
  rownames(metadata_df) = metadata_df$name
  return(metadata_df)
}

metadata_df = read_metadata("metadata.csv")

interaction_counts_df = read_counts("/Volumes/chris_ssd/projects/external/sRNA_mRNA_interactome_hyderabad/sequencing_run_2/light_to_dark/05_postprocess/complete_super_interaction_transcript_counts.gct",
                          levels(metadata_df$name))

interaction_counts_df_t <-  as.data.frame(t(interaction_counts_df), stringsAsFactors = FALSE)


df_int <- counts_only %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns  = ~ as.integer(round(.))
  ))

pca_res = prcomp(counts_only_t)

fviz_eig(pca_res)
fviz_pca_ind(pca_res,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)



dds = DESeqDataSetFromMatrix(df_int, colData = coldata, design = ~condition)
keep <- rowSums(counts(dds)) >= 5

keep2 <- rowSums(counts(dds) > 0) >= 2

combined = keep & keep2

summary(combined)

dds <- dds[combined,]

dds <- DESeq(dds)

res <- results(dds)

resLFC <- lfcShrink(dds, coef="condition_HighLight_vs_Control")

summary(resLFC, alpha=0.05)

sum(resLFC$padj < 0.05, na.rm=TRUE)

plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05)

EnhancedVolcano(resLFC,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05)

resOrdered <- res[order(res$pvalue),]


select <- order(res$padj)
select_padj = res$padj[select] <= 0.05

select = subset(select, select_padj)

df <- as.data.frame(colData(dds)[,c("condition", "name")])

vsd <- vst(dds)

mat <- assay(vsd)[ select, ]
rownames(df)  <- colnames(mat)

pheatmap(mat, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)[ select, ]))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$name, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

write.csv(as.data.frame(resOrdered), file = "condition_HighLight_vs_Control.tsv", )

sig_res = subset(resOrdered, padj < 0.05)

interactions_bed = read_delim("/Volumes/chris_ssd/projects/external/sRNA_mRNA_interactome_hyderabad/sequencing_run_2/light_to_dark/05_postprocess/complete_super_interaction_regions.bedpe", 
                                                          delim = "\t", escape_double = FALSE, 
                                                          trim_ws = TRUE, skip = 2, col_names = FALSE)

filtered_interactions_bed <- interactions_bed[interactions_bed[[7]] %in% sig_res@rownames, ]
write_tsv(filtered_interactions_bed,  file = "condition_HighLight_vs_Control_sig_padj05.bedpe", )

target_interaction = interactions_bed[(grepl("ncr0700", interactions_bed$X12) | grepl("ncr0700", interactions_bed$X13)), ]$X7

norm_counts = counts(dds, normalized = TRUE)

subsample_rownames = norm_counts[rownames(norm_counts) %in% target_interaction, ]

subsample_rownames = subsample_rownames[order(rowMeans(subsample_rownames)), ]
subsample_rownames = log(subsample_rownames + 1)

mat <- assay(subsample_rownames)
rownames(df)  <- colnames(mat)

pheatmap(subsample_rownames, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


target = res[rownames(res) %in% target_interaction, ]
target = target[order(target$padj),]
