# ==============================================
# Fig. 2: Bulk vs scRNA-seq (scRR-seq & G&T-seq)
# ==============================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(UpSetR)

# --- Fig. 2a: PCA - Bulk vs scRNA-seq ---
# Load expression data
df <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig2/data/scRNAsvsBulk_rsem_tpm_mouse_8cell_embryo.txt")
head(df)
# log10 transformation
df_log <- log10(t(df+1))

# PCA
df_log_pc <- prcomp(df_log)
df_log_pca <- data.frame(df_log_pc$x)

# Calculate explained variance
pca_var <- df_log_pc$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 2)

# Annotate data
library(dplyr)
df_log_pca <- df_log_pca %>%
  mutate(
    stage = case_when(
      grepl("scRR", rownames(df_log_pca)) ~ "8cell",
      grepl("GTseq", rownames(df_log_pca)) ~ "8cell",
      grepl("1cell", rownames(df_log_pca)) ~ "1cell",
      grepl("2cell", rownames(df_log_pca)) ~ "2cell",
      grepl("4cell", rownames(df_log_pca)) ~ "4cell",
      grepl("8cell", rownames(df_log_pca)) ~ "8cell",
      grepl("blastocyst", rownames(df_log_pca)) ~ "blastocyst"),
    method = case_when(
      grepl("scRR", rownames(df_log_pca)) ~ "scRR",
      grepl("GTseq", rownames(df_log_pca)) ~ "GTSeq",
      grepl("SRR", rownames(df_log_pca)) ~ "Bulk"),
    exp = case_when(
      grepl("scRR3", rownames(df_log_pca)) ~ "scRR3",
      grepl("scRR1", rownames(df_log_pca)) ~ "scRR1",
      grepl("GTseq", rownames(df_log_pca)) ~ "GTSeq",
      grepl("SRR", rownames(df_log_pca)) ~ "Bulk")
  )

# PCA plot
ggplot(df_log_pca, aes(PC1, PC2, col = stage, shape = exp))+
  geom_point(alpha = 0.7) +
  xlab(paste0("PC1=",pca_var_percent[1],"%"))+
  ylab(paste0("PC2=",pca_var_percent[2],"%"))+
  theme_classic()+
  scale_shape_manual(values = c("Bulk" = 19, "GTSeq" = 8, "scRR3" = 12, "scRR1" = 11))  # Specify shapes

#
write.table(df_log_pca[1:2], "/Volumes/scRR-seq_GEO/source_data/source_data_fig2/fig2a_pca.txt",
            sep = "\t", quote = F, col.names = T, row.names = T)


# --- Fig. 2b: Detected Transcripts and Genes ---

# Load RSEM TPM data
rsem_iso <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig2/data/scRR-seq-RNA_8cell_mouse_transcripts_genes_rsem_selectedsamples_TPM_mm10.txt")
rsem_gene <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig2/data/scRR-seq-RNA_8cell_mouse_genes_rsem_selectedsamples_TPM_mm10.txt")

# Set rownames and remove ID columns
rownames(rsem_iso) <- rsem_iso$transcript.id
rownames(rsem_gene) <- rsem_gene$gene.id
rsem_iso$transcript.id <- NULL
rsem_gene$gene.id <- NULL

# Define TPM thresholds
multithresholds <- seq(0, 5, by = 0.5)

# Count detected isoforms and genes above TPM thresholds
result_count <- NULL
for (thresh in multithresholds) {
  # Apply the function to calculate the count for each column
  iso_count_values <- apply(rsem_iso[, 1:ncol(rsem_iso)], 2, function(x) sum(x > thresh))
  gene_count_values <- apply(rsem_gene[, 1:ncol(rsem_gene)], 2, function(x) sum(x > thresh))
  # Append the results to the data frame
  result_count <- rbind(result_count, data.frame(threshold = thresh, iso_count = iso_count_values, gene_count = gene_count_values, Sample = names(iso_count_values)))
}


# Annotate samples with experiment type
anno <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig2/data/Annotation_8cell_embryo.txt")
result_count2 <- merge(result_count, anno, by = "Sample", sort = F)

# Compute mean and standard deviation by group
result_iso_merge <- as.data.table(aggregate(iso_count ~ threshold + exp_RNA, result_count2, function(x) c(Mean = mean(x), SD = sd(x))))
colnames(result_iso_merge) <- c("threshold","exp_RNA","Mean","SD")
result_gene_merge <- as.data.table(aggregate(gene_count ~ threshold + exp_RNA, result_count2, function(x) c(Mean = mean(x), SD = sd(x))))
colnames(result_gene_merge) <- c("threshold","exp_RNA","Mean","SD")

print(result_iso_merge)
print(result_gene_merge)


# Plot
ggplot(result_iso_merge, aes(x = threshold, y = Mean, group = exp_RNA))+
  #geom_point(aes(col = Method)) +
  geom_line(aes(col = exp_RNA)) +
  geom_vline(xintercept = 1, col = "gray", linetype = "dashed" ) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = exp_RNA), alpha = 0.3 ) +
  ylab("Detected transcripts") + xlab("TMP >") + ylim(0,40000) +
  theme_classic()

ggplot(result_gene_merge, aes(x = threshold, y = Mean, group = exp_RNA))+
  #geom_point(aes(col = Method)) +
  geom_line(aes(col = exp_RNA)) +
  geom_vline(xintercept = 1, col = "gray", linetype = "dashed" ) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = exp_RNA), alpha = 0.3 ) +
  ylab("Detected genes") + xlab("TMP >") +ylim(0,20000) +
  theme_classic()

#
write.table(result_iso_merge, "/Volumes/scRR-seq_GEO/source_data/source_data_fig2/fig2b_result_iso_merge.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(result_gene_merge, "/Volumes/scRR-seq_GEO/source_data/source_data_fig2/fig2b_esult_gene_merge.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)



# --- Fig. 2c: Find shared expressed genes (TPM > 1, 100% samples) ***among 8-cell stage*** ---

# Load data
df <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig2/data/scRNAsvsBulk_rsem_tpm_mouse_8cell_embryo.txt")
colnames(df)
# Subset by experiment type
scRR3 <- df[ , grep("scRR3", colnames(df))]
scRR1 <- df[ , grep("scRR1", colnames(df))]
GTseq <- df[ , grep("GTseq", colnames(df))]
bulk <- df[ , grep("BL6J_8cell", colnames(df))]

# Subset genes with TPM > 1 in 100% of samples
data <- df
filter_genes <- function(data, threshold_percentage) {
  # Calculate the threshold count based on the percentage
  threshold_count <- ncol(data) * (threshold_percentage / 100)
  # Apply the filtering condition
  filtered_genes <- rownames(data)[rowSums(data > 1) >= threshold_count]
  
  return(filtered_genes)
}

# Set the threshold percentage
threshold_percentage <- 100

# Apply filter
selected_genes_scRR3 <- filter_genes(scRR3, threshold_percentage)
selected_genes_scRR1 <- filter_genes(scRR1, threshold_percentage)
selected_genes_GTseq <- filter_genes(GTseq, threshold_percentage)
selected_genes_bulk <- filter_genes(bulk, threshold_percentage)

# Create list of gene sets
gene_sets <- list(scRR3 = selected_genes_scRR3,scRR1 = selected_genes_scRR1, GTseq = selected_genes_GTseq, Bulk = selected_genes_bulk)

# Plot shared genes using UpSet
shared <- gene_sets
upset_df <- upset(fromList(shared), order.by = "freq")

#





