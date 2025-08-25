# ========================================================================================
# Fig. 5b (Part 3): Correlation between replicated/unreplicated status and gene expression
# ========================================================================================

# For spliced data  -----
# Set output directory and create if needed
output_dir <- "/Volumes/scRR-seq_GEO/codes/Fig5/data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Load RT and spliced RNA counts --------

# Load replication timing (RT) genes data
RT_genes <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_avrRTperGene.txt")
rownames(RT_genes) <- RT_genes$Geneid
RT_genes$Geneid <- NULL

# Load spliced RNA counts
spliced_midS_df <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_spliced.txt")

# Ensure rownames are set for spliced data
rownames(spliced_midS_df) <- spliced_midS_df$Geneid

# Select genes common to RT and spliced data
common_genes <- intersect(rownames(RT_genes), rownames(spliced_midS_df))
spliced_midS_sel <- spliced_midS_df[common_genes,]
RT_genes_sel <- RT_genes[common_genes,]

# Confirm matching order
stopifnot(all(rownames(RT_genes_sel) == rownames(spliced_midS_sel)))

# -------- Filtering and adjust data before comparing --------

# Filter low-abundance genes (keep genes with counts > 0 in more than 28 samples)
gene_expr_mat <- spliced_midS_sel[, 8:ncol(spliced_midS_sel)]
to_keep <- rowSums(gene_expr_mat > 0) > 28
filtered_mat_gene <- gene_expr_mat[to_keep, ]
filtered_mat_rt <- RT_genes_sel[to_keep, ]

# Adjust gene expression matrix (can apply log1p if desired)
norm_mat <- filtered_mat_gene

# Adjust RT matrix: replace 0 with NA and -1 with 0
rt_mat <- filtered_mat_rt
rt_mat[rt_mat == 0] <- NA
rt_mat[rt_mat == -1] <- 0


# -------- Comparing gene expression between groups --------
# Threshold for minimum group size (~20% of samples)
min_group_size <- round(ncol(norm_mat) * 0.2)

# Initialize result containers
med_df <- NULL
listofdf <- list()

# Loop over each gene
for (i in seq_len(nrow(norm_mat))) {
  rt_values <- rt_mat[i, ]
  expr_values <- norm_mat[i, ]
  
  rep_samples <- names(rt_values)[which(rt_values == 1)]
  unrep_samples <- names(rt_values)[which(rt_values == 0)]
  
  expr_rep <- as.numeric(expr_values[rep_samples])
  expr_unrep <- as.numeric(expr_values[unrep_samples])
  
  # Calculate p-value if both groups have enough samples
  if (length(expr_rep) < min_group_size || length(expr_unrep) < min_group_size) {
    pval <- NA
  } else {
    pval <- wilcox.test(expr_rep, expr_unrep)$p.value
  }
  
  # Compile summary statistics for this gene
  gene_summary <- data.frame(
    Geneid = rownames(norm_mat)[i],
    rep_bins = length(expr_rep),
    unrep_bins = length(expr_unrep),
    unrep_ratio = length(expr_unrep) / (length(expr_rep) + length(expr_unrep)),
    rep_med = median(expr_rep),
    unrep_med = median(expr_unrep),
    rep_mean = mean(expr_rep),
    unrep_mean = mean(expr_unrep),
    rep_sd = sd(expr_rep),
    unrep_sd = sd(expr_unrep),
    rep_var = var(expr_rep),
    unrep_var = var(expr_unrep),
    pvalue = pval,
    gene_med = median(as.numeric(expr_values)),
    gene_mean = mean(as.numeric(expr_values))
  )
  
  med_df <- rbind(med_df, gene_summary)
  listofdf[[rownames(norm_mat)[i]]] <- list(rep = expr_rep, unrep = expr_unrep)
}

# Export results
write.table(filtered_mat_gene, paste0(output_dir, "/spliced_informative_genes.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(filtered_mat_rt, paste0(output_dir, "/spliced_informative_RT.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(med_df, paste0(output_dir, "/spliced_informative_summary_genes_RT_Median_Mean_ttest.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(listofdf, paste0(output_dir, "/spliced_informative_rawnumber_RT_Median_Mean_ttest.rds"))

# -------- Summarize results --------
# Remove rows with NA p-values (stable RT bins)
med_df_rmna <- na.omit(med_df)

# Summarize conditions
summary_table <- data.frame(
  Condition = c(
    "Total (NA removed)",
    "unrep < rep",
    "unrep > rep",
    "unrep < rep & p < 0.05",
    "unrep > rep & p < 0.05",
    "unrep < rep & p < 0.01",
    "unrep > rep & p < 0.01"
  ),
  Count = c(
    nrow(med_df_rmna),
    nrow(subset(med_df_rmna, unrep_med < rep_med)),
    nrow(subset(med_df_rmna, unrep_med > rep_med)),
    nrow(subset(med_df_rmna, unrep_med < rep_med & pvalue < 0.05)),
    nrow(subset(med_df_rmna, unrep_med > rep_med & pvalue < 0.05)),
    nrow(subset(med_df_rmna, unrep_med < rep_med & pvalue < 0.01)),
    nrow(subset(med_df_rmna, unrep_med > rep_med & pvalue < 0.01))
  )
)

print(summary_table)
write.table(summary_table, paste0(output_dir, "/spliced_informative_summary_table.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



# For unspliced data -----
# Set output directory and create if needed
output_dir <- "/Volumes/scRR-seq_GEO/codes/Fig5/data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Load RT and spliced RNA counts --------

# Load replication timing (RT) genes data
RT_genes <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_avrRTperGene.txt")
rownames(RT_genes) <- RT_genes$Geneid
RT_genes$Geneid <- NULL

# Load spliced RNA counts
unspliced_midS_df <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_unspliced.txt")

# Ensure rownames are set for spliced data
rownames(unspliced_midS_df) <- unspliced_midS_df$Geneid

# Select genes common to RT and spliced data
common_genes <- intersect(rownames(RT_genes), rownames(unspliced_midS_df))
unspliced_midS_sel <- unspliced_midS_df[common_genes,]
RT_genes_sel <- RT_genes[common_genes,]

# Confirm matching order
stopifnot(all(rownames(RT_genes_sel) == rownames(unspliced_midS_sel)))

# -------- Filtering and adjust data before comparing --------

# Filter low-abundance genes (keep genes with counts > 0 in more than 28 samples)
gene_expr_mat <- unspliced_midS_sel[, 8:ncol(unspliced_midS_sel)]
to_keep <- rowSums(gene_expr_mat > 0) > 28
filtered_mat_gene <- gene_expr_mat[to_keep, ]
filtered_mat_rt <- RT_genes_sel[to_keep, ]

# Adjust gene expression matrix (can apply log1p if desired)
norm_mat <- filtered_mat_gene

# Adjust RT matrix: replace 0 with NA and -1 with 0
rt_mat <- filtered_mat_rt
rt_mat[rt_mat == 0] <- NA
rt_mat[rt_mat == -1] <- 0


# -------- Comparing gene expression between groups --------
# Threshold for minimum group size (~20% of samples)
min_group_size <- round(ncol(norm_mat) * 0.2)

# Initialize result containers
med_df <- NULL
listofdf <- list()

# Loop over each gene
for (i in seq_len(nrow(norm_mat))) {
  rt_values <- rt_mat[i, ]
  expr_values <- norm_mat[i, ]
  
  rep_samples <- names(rt_values)[which(rt_values == 1)]
  unrep_samples <- names(rt_values)[which(rt_values == 0)]
  
  expr_rep <- as.numeric(expr_values[rep_samples])
  expr_unrep <- as.numeric(expr_values[unrep_samples])
  
  # Calculate p-value if both groups have enough samples
  if (length(expr_rep) < min_group_size || length(expr_unrep) < min_group_size) {
    pval <- NA
  } else {
    pval <- wilcox.test(expr_rep, expr_unrep)$p.value
  }
  
  # Compile summary statistics for this gene
  gene_summary <- data.frame(
    Geneid = rownames(norm_mat)[i],
    rep_bins = length(expr_rep),
    unrep_bins = length(expr_unrep),
    unrep_ratio = length(expr_unrep) / (length(expr_rep) + length(expr_unrep)),
    rep_med = median(expr_rep),
    unrep_med = median(expr_unrep),
    rep_mean = mean(expr_rep),
    unrep_mean = mean(expr_unrep),
    rep_sd = sd(expr_rep),
    unrep_sd = sd(expr_unrep),
    rep_var = var(expr_rep),
    unrep_var = var(expr_unrep),
    pvalue = pval,
    gene_med = median(as.numeric(expr_values)),
    gene_mean = mean(as.numeric(expr_values))
  )
  
  med_df <- rbind(med_df, gene_summary)
  listofdf[[rownames(norm_mat)[i]]] <- list(rep = expr_rep, unrep = expr_unrep)
}

# Export results
write.table(filtered_mat_gene, paste0(output_dir, "/unspliced_informative_genes.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(filtered_mat_rt, paste0(output_dir, "/unspliced_informative_RT.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(med_df, paste0(output_dir, "/unspliced_informative_summary_genes_RT_Median_Mean_ttest.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(listofdf, paste0(output_dir, "/unspliced_informative_rawnumber_RT_Median_Mean_ttest.rds"))

# -------- Summarize results --------
# Remove rows with NA p-values (stable RT bins)
med_df_rmna <- na.omit(med_df)

# Summarize conditions
summary_table <- data.frame(
  Condition = c(
    "Total (NA removed)",
    "unrep < rep",
    "unrep > rep",
    "unrep < rep & p < 0.05",
    "unrep > rep & p < 0.05",
    "unrep < rep & p < 0.01",
    "unrep > rep & p < 0.01"
  ),
  Count = c(
    nrow(med_df_rmna),
    nrow(subset(med_df_rmna, unrep_med < rep_med)),
    nrow(subset(med_df_rmna, unrep_med > rep_med)),
    nrow(subset(med_df_rmna, unrep_med < rep_med & pvalue < 0.05)),
    nrow(subset(med_df_rmna, unrep_med > rep_med & pvalue < 0.05)),
    nrow(subset(med_df_rmna, unrep_med < rep_med & pvalue < 0.01)),
    nrow(subset(med_df_rmna, unrep_med > rep_med & pvalue < 0.01))
  )
)

print(summary_table)
write.table(summary_table, paste0(output_dir, "/unspliced_informative_summary_table.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

