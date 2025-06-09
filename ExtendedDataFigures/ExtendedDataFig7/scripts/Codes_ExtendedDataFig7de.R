# ==================================================================================
# Extended Data Fig. 7d and e: Expression Allelic Ratios of Autosomal Genes in CBMS1
# ==================================================================================

# Load required libraries
options(scipen = 100)
library(dplyr)
library(reshape2)
library(ggplot2)

# Load gene annotation
asset_mm10 <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/asset.Mus_musculus.GRCm38.75.Coord.ID.name.txt")

# Load haplotype-specific read counts
cba <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/03_scRR-seq-RNA_CBMS1_mESC_allelic_merged_emase_cba.gene.alignment.counts_mm10.txt")
msm <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/03_scRR-seq-RNA_CBMS1_mESC_allelic_merged_emase_msm.gene.alignment.counts_mm10.txt")

rownames(cba) <- cba$X.target_id
cba$X.target_id <- NULL

rownames(msm) <- msm$X.target_id
msm$X.target_id <- NULL

# Sum total allele counts per gene across both haplotypes
cbms1_sum <- cba+msm

# Identified genes that are informative (SNP-containing reads ≥ 10, in more that 50% of cells)
threshold_detection = 10
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

CBMS1_scRR_set <- filter_reads(cbms1_sum, asset_mm10)
dim(CBMS1_scRR_set) #[1] 6530   89

write.table(CBMS1_scRR_set, paste0("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/informative_genes_CBMS1_all_set.txt"), sep = "\t", row.names = F, col.names = T, quote = F)


# Calculate allele expression ratios (a / (a + b))
cba_tpm <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/03_scRR-seq-RNA_CBMS1_mESC_allelic_merged_emase_cba.gene.tpm_mm10.txt")
msm_tpm <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/03_scRR-seq-RNA_CBMS1_mESC_allelic_merged_emase_msm.gene.tpm_mm10.txt")

rownames(cba_tpm) <- cba_tpm$X.target_id
cba_tpm$X.target_id <- NULL
rownames(msm_tpm) <- msm_tpm$X.target_id
msm_tpm$X.target_id <- NULL

ratio_df <- cba_tpm/(cba_tpm+msm_tpm)

# Annotate with gene coordinates
ratio_df$gene_id <- rownames(ratio_df)
ratio_df_gene_list <- merge(asset_mm10, ratio_df, by = "gene_id")
rownames(ratio_df_gene_list) <- ratio_df_gene_list$gene_id

# Select informative genes
ratio_df_gene_list_informative <- ratio_df_gene_list[CBMS1_scRR_set$gene_id,]
dim(ratio_df_gene_list_informative) #[1] 6530   89

# Find median and mean ratio for each gene
selected <- ratio_df_gene_list_informative
selected$median_ratio <- apply(ratio_df_gene_list_informative[6:ncol(ratio_df_gene_list_informative)], 1, function(x) median(x, na.rm = TRUE))
selected$mean_ratio <- apply(ratio_df_gene_list_informative[6:ncol(ratio_df_gene_list_informative)], 1, function(x) mean(x, na.rm = TRUE))

# Classify allelic expression patterns
selected <- selected %>%
  mutate(category_median = case_when(
    median_ratio >= 0.85 | median_ratio <= 0.15 ~ "monoallelic",
    median_ratio >= 0.7 & median_ratio < 0.85 ~ "biased_cba",
    median_ratio >= 0.15 & median_ratio <= 0.3 ~ "biased_msm",
    median_ratio > 0.3 & median_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

selected <- selected %>%
  mutate(category_mean = case_when(
    mean_ratio >= 0.85 | mean_ratio <= 0.15 ~ "monoallelic",
    mean_ratio >= 0.7 & mean_ratio < 0.85 ~ "biased_cba",
    mean_ratio >= 0.15 & mean_ratio <= 0.3 ~ "biased_msm",
    mean_ratio > 0.3 & mean_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

write.table(selected, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression.txt", sep = "\t", row.names = F, col.names = T, quote = F)
#

# Summarize gene categories [Median]
counts <- selected %>%
  count(category_median)
print(counts)

write.table(counts, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression_category_median.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_median    n
#1       biallelic 5980
#2      biased_cba  190
#3      biased_msm  150
#4     monoallelic  210

# # Summarize gene categories [Mean]
counts <- selected %>%
  count(category_mean)
print(counts)

write.table(counts, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression_category_mean.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_mean    n
#1     biallelic 6153
#2    biased_cba  148
#3    biased_msm  117
#4   monoallelic  112


# Remove abnormal chromosome (chr8 and chrX)
new_df_selected_rmna_filterd <- selected[!grepl("chr8|X", selected$chr),]
write.table(new_df_selected_rmna_filterd, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression_rmchr8chrX.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Count categories after filtering [Median]
counts_filtered <- new_df_selected_rmna_filterd %>%
  count(category_median)
print(counts_filtered)
write.table(counts_filtered, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression_category_median_rmchr8chrX.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_median    n
#1       biallelic 5512
#2      biased_cba  171
#3      biased_msm  134
#4     monoallelic  155

# Count categories after filtering [Mean]
counts_filtered <- new_df_selected_rmna_filterd %>%
  count(category_mean)
print(counts_filtered)

write.table(counts_filtered, "/Volumes/scRR-seq_GEO/codes/ExtendedDataFig7/data/CBMS1_mono_bias_bi_expression_category_mean_rmchr8chrX.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#category_mean    n
#1     biallelic 5622
#2    biased_cba  142
#3    biased_msm  105
#4   monoallelic  103


# Plot for imprinted genes
list <- c("Peg3","Peg10", "Snrpn", "Rian",  "Igf2r", "Slc38a4")

df_plot <- new_df_selected_rmna_filterd[ which(new_df_selected_rmna_filterd$gene_name %in% list),]
df_plot$chr <- NULL
df_plot$start <- NULL
df_plot$end <- NULL
df_plot_melt <- melt(df_plot, id.vars = c("gene_id", "category_median", "gene_name", "median_ratio",
                                          "category_mean","mean_ratio"))

# Order according to repliscores
repliscores <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/03_scRR-seq-DNA_CBMS1_mESC_pctreplicationscore_selectedsamples_all_mm10.txt")

sort_samples <- repliscores$Sample[order(repliscores$pct_repliscore, decreasing = T)]
df_plot_melt$variable <- factor(df_plot_melt$variable, levels = sort_samples)

# Plot
ggplot(df_plot_melt, aes(x = value, y = variable)) +
  geom_point() +  # Plot points
  facet_grid(~gene_name, scales = "free_y") +
  labs(x = "Gene", y = "ratio") + # Labels for axes
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") + # Add red dashed line at y = 0.5
  theme_bw()

#


