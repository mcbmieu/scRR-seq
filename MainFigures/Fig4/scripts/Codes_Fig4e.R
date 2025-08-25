# ==============================================================================
# Fig. 4e: Allelic Gene Expression Comparison in CBMS1 Bulk RNA-seq vs scRR-seq
# ==============================================================================

# Load required libraries
library(dplyr)
library(eulerr)

# Load data
cbms1_bulk <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/GSE108556_CBMS1_allele_specific_genetable.txt", skip = 10)

# Compute mean allelic ratio and mean number of informative SNPs
cbms1_bulk$mean_ratio <- rowMeans(cbms1_bulk[6:8])
cbms1_bulk$mean_snps <- rowMeans(cbms1_bulk[9:11])

# Filter out genes on chr8 and chrX (due to known abnormalities)
cbms1_bulk_rm <- cbms1_bulk[cbms1_bulk$X.chrom != "chr8" & cbms1_bulk$X.chrom != "chrX",]

# Classify gene allelic patterns (based on mean_ratio)
cbms1_bulk_rm <- cbms1_bulk_rm %>%
  mutate(category = case_when(
    mean_ratio >= 0.85 | mean_ratio <= 0.15 ~ "monoallelic",
    mean_ratio >= 0.7 & mean_ratio < 0.85 ~ "biased_cba",
    mean_ratio >= 0.15 & mean_ratio <= 0.3 ~ "biased_msm",
    mean_ratio > 0.3 & mean_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

# Filter for genes with strong SNP data (mean_snps > 10)
cbms1_bulk_snp10 <- cbms1_bulk_rm[cbms1_bulk_rm$mean_snps > 10,] #[1] 8231   17

# Export
write.table(cbms1_bulk_snp10, "/Volumes/scRR-seq_GEO/codes/Fig4/data/informative_cbms1_bulk_snp10.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Find overlap of informative genes
scrr_informative_genes <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig13/data/CBMS1_mono_bias_bi_expression_rmchr8chrX.txt")
all_gene_sets <- list(cbms1_bulk_snp10 = unique(cbms1_bulk_snp10$X.name),
                      scrr_informative_genes = unique(scrr_informative_genes$gene_name))

# Create a Venn diagram
euler_plot <- euler(all_gene_sets)
plot(euler_plot, legend = T, quantities = TRUE)

## 

