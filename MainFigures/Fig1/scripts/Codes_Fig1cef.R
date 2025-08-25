# ============================================
# Fig. 1: Detect genes and transcripts
# ============================================

# --- Fig. 1c: Coverage along transcripts ---

# Load additional libraries
library(reshape2)

# Load coverage data
# If analysis done by ramdaq, coverage file is stored in the following location
# "Multiqc/multiqc_data/mqc_rseqc_gene_body_coverage_plot_Percentages.txt"
coverage <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data//mqc_rseqc_gene_body_coverage_plot_Percentages_RPE1midS.txt")

# Reshape and clean coverage data
coverage_re <- melt(coverage)
coverage_re$quantile <- as.numeric(gsub("X","",coverage_re$variable))

# Plot: Transcript coverage
ggplot(coverage_re, aes(quantile, value, col = exp, group = Sample))+
  #geom_point(alpha = 0.3) + 
  geom_line(alpha = 0.3) +
  ylab("Mean read coverage (%)") + xlab(" Transcript (%)") +
  theme_classic()

#
write.table(coverage_re, "/Volumes/scRR-seq_GEO/source_data/source_data_fig1/fig1c_coverage.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)


# --- Fig. 1e: Detected isoforms and genes ---

# Load required libraries
library(dplyr)
library(data.table)
library(ggplot2)

# Load expression data
rsem_iso <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data/01_scRR-seq-RNA_RPE1_midS_transcripts_rsem_selectedsamples_TPM_hg38.txt", header = T)
rsem_gene <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data/01_scRR-seq-RNA_RPE1_midS_genes_rsem_selectedsamples_TPM_hg38.txt", header = T)

# Set rownames and remove ID columns
rownames(rsem_iso) <- rsem_iso$transcript.id
rownames(rsem_gene) <- rsem_gene$Gene.id
rsem_iso$transcript.id <- NULL
rsem_gene$Gene.id <- NULL

# Define TPM thresholds
multithresholds <- seq(0, 5, by = 0.5)

# Count detected isoforms and genes above TPM thresholds
result_count <- NULL
for (thresh in multithresholds) {
  iso_count_values <- apply(rsem_iso[, 1:ncol(rsem_iso)], 2, function(x) sum(x > thresh))
  gene_count_values <- apply(rsem_gene[, 1:ncol(rsem_gene)], 2, function(x) sum(x > thresh))
  result_count <- rbind(result_count, data.frame(threshold = thresh, iso_count = iso_count_values, gene_count = gene_count_values, Sample = names(iso_count_values)))
}

# Annotate samples with experiment type
library(dplyr)
result_count2 <- result_count %>%
  mutate(exp = case_when(
    grepl("scRR3", Sample) ~ "scRR3",
    grepl("scRR1", Sample) ~ "scRR1",
    grepl("scRamDA", Sample) ~ "scRamDA")
    )

# Compute mean and standard deviation by group
result_iso_merge <- as.data.table(aggregate(iso_count ~ threshold + exp, result_count2, function(x) c(Mean = mean(x), SD = sd(x))))
colnames(result_iso_merge) <- c("threshold","exp","Mean","SD")
result_gene_merge <- as.data.table(aggregate(gene_count ~ threshold + exp, result_count2, function(x) c(Mean = mean(x), SD = sd(x))))
colnames(result_gene_merge) <- c("threshold","exp","Mean","SD")

print(result_iso_merge)
print(result_gene_merge)


# Plot
ggplot(result_iso_merge, aes(x = threshold, y = Mean, group = exp))+
  #geom_point(aes(col = Method)) +
  geom_line(aes(col = exp)) +
  geom_vline(xintercept = 1, col = "gray", linetype = "dashed" ) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = exp), alpha = 0.3 ) +
  ylab("Detected transcripts") + xlab("TMP >") + ylim(0,65000) +
  theme_classic()

ggplot(result_gene_merge, aes(x = threshold, y = Mean, group = exp))+
  #geom_point(aes(col = Method)) +
  geom_line(aes(col = exp)) +
  geom_vline(xintercept = 1, col = "gray", linetype = "dashed" ) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = exp), alpha = 0.3 ) +
  ylab("Detected genes") + xlab("TMP >") +ylim(0,25000) +
  theme_classic()
#

write.table(result_iso_merge, "/Volumes/scRR-seq_GEO/source_data/source_data_fig1/fig1e_result_iso_merge.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(result_gene_merge, "/Volumes/scRR-seq_GEO/source_data/source_data_fig1/fig1e_result_gene_merge.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)










# --- Fig. 1f: Non-poly(A) ---
# Load expression data (using those from featureCounts)
nuclear_subset <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data/fig1f_nonpolyA_table.txt")
head(nuclear_subset)

nuclear_subset_df_log10 <- log10(nuclear_subset+1)
nuclear_subset_df_log10$median <- apply(nuclear_subset_df_log10, 1, median)
nuclear_subset_df_log10_sort <- nuclear_subset_df_log10[order(nuclear_subset_df_log10$median, decreasing = F),]
nuclear_subset_df_log10_sort_df <- nuclear_subset_df_log10_sort
nuclear_subset_df_log10_sort_df$median <- NULL


heatmap(as.matrix(nuclear_subset_df_log10_sort_df), 
        col = colorRampPalette(c('black','red','yellow','white'))(100),  # Define color palette
        scale = "none",  # Do not scale data
        Rowv = NA,  # Do not reorder rows
        Colv = NA,  # Do not reorder columns
        main = "Sensitivity for detecting non-poly(A) genes",
        breaks = seq(0, 5, length.out = 101), # Title of the plot
        cexRow = 0.5, cexCol = 0.3)

#
write.table(nuclear_subset_df_log10_sort_df, "/Volumes/scRR-seq_GEO/source_data/source_data_fig1/fig1f_nonpolyA_log10plus1.txt",
            sep = "\t", quote = F, col.names = T, row.names = T)

