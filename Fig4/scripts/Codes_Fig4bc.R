# ===================================================================================
# Fig. 4b and c: Expression Haplotype Ratios of Autosomal and X-linked Genes in RPE1
# ===================================================================================

# Load required libraries
options(scipen = 100)
library(dplyr)
library(reshape2)
library(ggplot2)

# Load gene annotation
asset_hg19 <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/asset.Homo_sapiens.GRCh37.75.Coord.ID.name.txt")

# Load haplotype-specific read counts
RPE1a <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/04_scRR-seq-RNA_RPE1_wholeS_merged_emase_haplotype_a.gene.alignment.counts_hg19.txt")
RPE1b <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/04_scRR-seq-RNA_RPE1_wholeS_merged_emase_haplotype_b.gene.alignment.counts_hg19.txt")

rownames(RPE1a) <- RPE1a$X.target_id
RPE1a$X.target_id <- NULL

rownames(RPE1b) <- RPE1b$X.target_id
RPE1b$X.target_id <- NULL

# Sum total allele counts per gene across both haplotypes
RPE1_sum <- RPE1a+RPE1b

# Identified genes that are informative (SNP-containing reads ≥ 6, in more that 50% of cells)
threshold_detection = 6
nCellRatio = 0.5

#  Get informative genes after filtering genes with low reads
filter_reads <- function(df, asset){
  nCell = round(dim(df)[2] * nCellRatio)
  print(sprintf("nCell: %d", nCell))
  selRowPre_df = df[rowSums(df > threshold_detection) > nCell,]
  informative_genes = rownames(selRowPre_df)
  
  informative_genes_df <- selRowPre_df
  informative_genes_df$gene_id <- rownames(selRowPre_df)
  informative_genes_final <- merge(asset, informative_genes_df)
}

RPE_scRR_set <- filter_reads(RPE1_sum, asset_hg19)
dim(RPE_scRR_set) #[1] 1979   54

write.table(RPE_scRR_set, paste0("/Volumes/scRR-seq_GEO/codes/Fig4/data/informative_genes_RPE_all_set.txt"), sep = "\t", row.names = F, col.names = T, quote = F)


# Calculate allele expression ratios (a / (a + b))
a_tpm <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/04_scRR-seq-RNA_RPE1_wholeS_merged_emase_haplotype_a.gene.tpm_hg19.txt")
b_tpm <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/04_scRR-seq-RNA_RPE1_wholeS_merged_emase_haplotype_b.gene.tpm_hg19.txt")

rownames(a_tpm) <- a_tpm$X.target_id
a_tpm$X.target_id <- NULL
rownames(b_tpm) <- b_tpm$X.target_id
b_tpm$X.target_id <- NULL

ratio_df <- a_tpm/(a_tpm+b_tpm)

# Annotate with gene coordinates
ratio_df$gene_id <- rownames(ratio_df)
ratio_df_gene_list <- merge(asset_hg19, ratio_df, by = "gene_id")
rownames(ratio_df_gene_list) <- ratio_df_gene_list$gene_id

# Select informative genes
ratio_df_gene_list_informative <- ratio_df_gene_list[RPE_scRR_set$gene_id,]
dim(ratio_df_gene_list_informative) #[1] 1979   54

# Find median and mean ratio for each gene
selected <- ratio_df_gene_list_informative
selected$median_ratio <- apply(ratio_df_gene_list_informative[6:ncol(ratio_df_gene_list_informative)], 1, function(x) median(x, na.rm = TRUE))
selected$mean_ratio <- apply(ratio_df_gene_list_informative[6:ncol(ratio_df_gene_list_informative)], 1, function(x) mean(x, na.rm = TRUE))

# Classify allelic expression patterns
selected <- selected %>%
  mutate(category_median = case_when(
    median_ratio >= 0.85 | median_ratio <= 0.15 ~ "monoallelic",
    median_ratio >= 0.7 & median_ratio < 0.85 ~ "rpe1_a",
    median_ratio >= 0.15 & median_ratio <= 0.3 ~ "rpe1_b",
    median_ratio > 0.3 & median_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

selected <- selected %>%
  mutate(category_mean = case_when(
    mean_ratio >= 0.85 | mean_ratio <= 0.15 ~ "monoallelic",
    mean_ratio >= 0.7 & mean_ratio < 0.85 ~ "rpe1_a",
    mean_ratio >= 0.15 & mean_ratio <= 0.3 ~ "rpe1_b",
    mean_ratio > 0.3 & mean_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

write.table(selected, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Summarize gene categories [Median]
counts <- selected %>%
  count(category_median)
print(counts)

write.table(counts, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression_category_median.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_median    n
#1       biallelic 1890
#2     monoallelic   39
#3          rpe1_a   33
#4          rpe1_b   17

# # Summarize gene categories [Mean]
counts <- selected %>%
  count(category_mean)
print(counts)

write.table(counts, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression_category_mean.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_mean    n
#1     biallelic 1905
#2   monoallelic   39
#3        rpe1_a   22
#4        rpe1_b   13


# Remove abnormal chromosome (chr10)
new_df_selected_rmna_filterd <- selected[!grepl("chr10", selected$chr),]
write.table(new_df_selected_rmna_filterd, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression_rmchr10.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Count categories after filtering [Median]
counts_filtered <- new_df_selected_rmna_filterd %>%
  count(category_median)
print(counts_filtered)
write.table(counts_filtered, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression_category_median_rmchr10.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_median    n
#1       biallelic 1808
#2     monoallelic   38
#3          rpe1_a   15
#4          rpe1_b   17

# Count categories after filtering [Mean]
counts_filtered <- new_df_selected_rmna_filterd %>%
  count(category_mean)
print(counts_filtered)

write.table(counts_filtered, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RPE1_mono_bias_bi_expression_category_mean_rmchr10.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_mean    n
#1     biallelic 1817
#2   monoallelic   38
#3        rpe1_a   10
#4        rpe1_b   13


# Plot haplotype ratio for selected imprinted genes
list <- c("MEG3","MEST", "PEG10", "PRIM2","IGF2R",  "DNMT1", "MCM5" )

df_plot <- new_df_selected_rmna_filterd[ which(new_df_selected_rmna_filterd$gene_name %in% list),]
df_plot$chr <- NULL
df_plot$start <- NULL
df_plot$end <- NULL
df_plot_melt <- melt(df_plot, id.vars = c("gene_id", "category_median", "gene_name", "median_ratio",
                                          "category_mean","mean_ratio"))

# Order according to repliscores
repliscores <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-DNA_RPE1_wholeS_pctreplicationscore_selectedsamples_all_hg38.txt")

sort_samples <- repliscores$Sample[order(repliscores$pct_repliscore, decreasing = F)]
df_plot_melt$variable <- factor(df_plot_melt$variable, levels = sort_samples)

# Plot 
ggplot(df_plot_melt, aes(x = value, y = variable)) +
  geom_point() +  # Plot points
  facet_grid(~gene_name, scales = "free_y") +
  labs(x = "Gene", y = "ratio") + # Labels for axes
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + # Add red dashed line at y = 0.5
  theme_bw()

# Plots for chrX only
chrX_only <- new_df_selected_rmna_filterd[new_df_selected_rmna_filterd$chr == "X",]
df_plot_chrX <- chrX_only
df_plot_chrX$chr <- NULL
df_plot_chrX$start <- NULL
df_plot_chrX$end <- NULL
df_plot_melt_chrX <- melt(df_plot_chrX, id.vars = c("gene_id", "category_median", "gene_name", "median_ratio",
                                                    "category_mean","mean_ratio"))

# Order according to repliscores
repliscores <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-DNA_RPE1_wholeS_pctreplicationscore_selectedsamples_all_hg38.txt")
sort_samples <- repliscores$Sample[order(repliscores$pct_repliscore, decreasing = F)]
df_plot_melt_chrX$variable <- factor(df_plot_melt_chrX$variable, levels = sort_samples)

# Plot
ggplot(df_plot_melt_chrX, aes(x = value, y = variable)) +
  geom_point() +  # Plot points
  facet_grid(~gene_name, scales = "free_y") +
  labs(x = "Gene", y = "ratio") + # Labels for axes
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + # Add red dashed line at y = 0.5
  theme_bw()

#
