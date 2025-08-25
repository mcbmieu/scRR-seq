# =========================================================
# Supplementary Fig. 5: Detected genes in RPE1 midS samples
# =========================================================

# --- Supplementary Fig. 5a (top panel): Mean expression level of detected genes ---

# Load libraries
options(scipen = 999)
library(dplyr)
library(rtracklayer)
library(ggplot2)


# Load data
rsem_gene <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data/01_scRR-seq-RNA_RPE1_midS_genes_rsem_selectedsamples_TPM_hg38.txt", header = T)
rownames(rsem_gene) <- rsem_gene$Gene.id

# Subset by experimental conditions
scRR1 <- rsem_gene[ , grep("scRR1", colnames(rsem_gene))]
scRR3 <- rsem_gene[ , grep("scRR3", colnames(rsem_gene))]
ramda <- rsem_gene[ , grep("scRamDA", colnames(rsem_gene))]

# Function to define formative genes: genes that have TPM > 0.1 in more than 50% of sample
informative_category <- function(df){
  thres = 0.1
  all = round(dim(df)[2])
  nCell = round(dim(df)[2]*0.5)
  df <- df %>%
    mutate(status = case_when(
      rowSums(df > thres) > nCell ~ "informative",
      rowSums(df == 0) == all ~ "nocounts",
      TRUE ~ "not_informative")
    )
}

# Apply to datasets
scRR1_status <- informative_category(scRR1)
scRR3_status <- informative_category(scRR3)
ramda_status <- informative_category(ramda)

# Calculate mean expression
scRR1_status$mean <- apply(scRR1, 1, mean)
scRR3_status$mean <- apply(scRR3, 1, mean)
ramda_status$mean <- apply(ramda, 1, mean)

# Summarize
scRR1_count <- scRR1_status %>%
  count(status)
scRR3_count <- scRR3_status %>%
  count(status)
ramda_count <- ramda_status %>%
  count(status)

# Merged status and mean data
merge_df <- data.frame(scRR1_status$status, scRR1_status$mean,
                       scRR3_status$status, scRR3_status$mean,
                       ramda_status$status, ramda_status$mean, Geneid = rownames(ramda_status))

## Remove genes with no counts in any condition
merge_df_nocount <- merge_df[merge_df$scRR1_status.status == "nocounts" & 
                               merge_df$scRR3_status.status == "nocounts" &
                               merge_df$ramda_status.status == "nocounts",]

merge_df_rmnocount <- merge_df[-which(merge_df$Geneid %in% merge_df_nocount$Geneid),]


# Calculate percentiles (based on scRamDA mean)
data = merge_df_rmnocount
percentiles <- sapply(1:nrow(data), function(i) {
  value <- data[i, 6]
  column_values <- data[, 6]
  percentile <- sum(column_values < value) / length(column_values) * 100
  return(percentile)
})

merge_df_rmnocount$percentile <- percentiles


# Genes informative in RamDA but NOT in others
merge_df_rmnocount_informative <- merge_df_rmnocount[merge_df_rmnocount$ramda_status.status == "informative" &
                                                     merge_df_rmnocount$scRR1_status.status != "informative" & 
                                                     merge_df_rmnocount$scRR3_status.status != "informative",]
boxplot(log10(merge_df_rmnocount_informative$ramda_status.mean+1),
        log10(merge_df_rmnocount_informative$scRR3_status.mean+1),
        log10(merge_df_rmnocount_informative$scRR1_status.mean+1),  ylim = c(0,5), names = c("Ramda", "scRR3", "scRR1"),
        main = paste0("informative genes found in RamDA,\n but not in scRR \n", dim(merge_df_rmnocount_informative)[1], " genes"))


# Get formative genes from scRamda AND also others
merge_df_rmnocount_informative_both <- merge_df_rmnocount[merge_df_rmnocount$ramda_status.status == "informative" &
                                                      merge_df_rmnocount$scRR1_status.status == "informative" & 
                                                      merge_df_rmnocount$scRR3_status.status == "informative",]
boxplot(log10(merge_df_rmnocount_informative_both$ramda_status.mean+1),
        log10(merge_df_rmnocount_informative_both$scRR3_status.mean+1),
        log10(merge_df_rmnocount_informative_both$scRR1_status.mean+1),  ylim = c(0,5), names = c("Ramda", "scRR3", "scRR1"),
        main = paste0("informative genes found in,\n both RamDA and scRR \n", dim(merge_df_rmnocount_informative_both)[1], " genes"))


# --- Supplementary Fig. 5b: Histogram of percentiles ---

hist(merge_df_rmnocount_informative_both$percentile, col=rgb(0,0,0,1/4), main = "hist of percentile", xlim = c(0,100))
hist(merge_df_rmnocount_informative$percentile, col=rgb(1,0,0,1/4), add = T)

# --- Supplementary Fig. 5c: Gene type enrichment ----

onlyramda <- merge_df_rmnocount_informative$Geneid
detectboth <- merge_df_rmnocount_informative_both$Geneid

# Gene annotation
asset <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.human.gencode.v38.primary_assembly.annotation.txt")

# Select only subset of genes
onlyramda_name <- asset[asset$Geneid %in% onlyramda,]
detectboth_name <- asset[asset$Geneid %in% detectboth,]

# Check gene_type for these genes
gtf_file <- "/Volumes/scRR-seq_GEO/codes/asset/gencode.v38.primary_assembly.annotation.gtf"
gtf_data <- import.gff(gtf_file)

# Extract gene_id and gene_name
gene_info <- subset(gtf_data, type == "gene", select = c("gene_id", "gene_name", "gene_type"))

# Remove duplicates (if any)
gene_info <- unique(gene_info)

# Extract gene list
gene_list_ramda <- onlyramda_name$Geneid
gene_list_both <- detectboth_name$Geneid

# Select only subset of genes
subset_gencode_ramda <- subset(gene_info, gene_id %in% gene_list_ramda)
subset_gencode_both <- subset(gene_info, gene_id %in% gene_list_both)

# Plot
gene_types_ramda <- data.frame(subset_gencode_ramda$gene_type)
gene_type_ramda_counts <- table(gene_types_ramda)

gene_types_both <- data.frame(subset_gencode_both$gene_type)
gene_type_both_counts <- table(gene_types_both)

# Convert gene_type_counts to data frame for easier plotting
gene_type_df_ramda <- as.data.frame(gene_type_ramda_counts, stringsAsFactors = FALSE)
gene_type_df_both <- as.data.frame(gene_type_both_counts, stringsAsFactors = FALSE)

# Plot
ggplot(gene_type_df_ramda, aes(x = reorder(subset_gencode_ramda.gene_type, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = Freq), hjust = -0.5 ,vjust = 0.5, size = 3.5, color = "black") +  # Add count annotations
  ylim(0,1800) +
  labs(x = "Gene Type", y = "Count", title = "Number of Genes by Gene Type: RamDA-seq") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

ggplot(gene_type_df_both, aes(x = reorder(subset_gencode_both.gene_type, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = Freq), hjust = -0.5 ,vjust = 0.5, size = 3.5, color = "black") +  # Add count annotations
  ylim(0,15000) +
  labs(x = "Gene Type", y = "Count", title = "Number of Genes by Gene Type:Both") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

#


# --- Supplementary Fig. 5a (middle and bottom panels): Mean expression level of detected lncRNAs and protein coding genes ---

# Subset only lncRNAs and protein coding genes
subset_gencode_ramda_lncRNA <- subset(subset_gencode_ramda, grepl("lncRNA", gene_type))
subset_gencode_both_lncRNA <- subset(subset_gencode_both, grepl("lncRNA", gene_type))

subset_gencode_ramda_protein <- subset(subset_gencode_ramda, grepl("protein_coding", gene_type))
subset_gencode_both_protein <- subset(subset_gencode_scRR, grepl("protein_coding", gene_type))

# Get expression data
data_subset_gencode_ramda_lncRNA <- data[data$Geneid%in%subset_gencode_ramda_lncRNA$gene_id,]
data_subset_gencode_both_lncRNA <- data[data$Geneid%in%subset_gencode_both_lncRNA$gene_id,]

data_subset_gencode_ramda_protein <- data[data$Geneid%in%subset_gencode_ramda_protein$gene_id,]
data_subset_gencode_both_protein <- data[data$Geneid%in%subset_gencode_both_protein$gene_id,]


boxplot(log10(data_subset_gencode_ramda_lncRNA$ramda_status.mean+1),
        log10(data_subset_gencode_ramda_lncRNA$scRR3_status.mean+1),
        log10(data_subset_gencode_ramda_lncRNA$scRR1_status.mean+1),  ylim = c(0,5), names = c("RamDA", "scRR3", "scRR1"),
        col=rgb(0,0,0,1/4),
        main = paste0("informative lncRNAs found in RamDA,\n but not in scRR \n", dim(data_subset_gencode_ramda_lncRNA)[1], "genes"))

boxplot(log10(data_subset_gencode_both_lncRNA$ramda_status.mean+1),
        log10(data_subset_gencode_both_lncRNA$scRR3_status.mean+1),
        log10(data_subset_gencode_both_lncRNA$scRR1_status.mean+1),  ylim = c(0,5), names = c("RamDA", "scRR3", "scRR1"),
        col=rgb(0,0,0,1/4),
        main = paste0("informative lncRNAs found in,\n both RamDA and scRR \n", dim(data_subset_gencode_both_lncRNA)[1], "genes"))

boxplot(log10(data_subset_gencode_ramda_protein$ramda_status.mean+1),
        log10(data_subset_gencode_ramda_protein$scRR3_status.mean+1),
        log10(data_subset_gencode_ramda_protein$scRR1_status.mean+1),  ylim = c(0,5), names = c("RamDA", "scRR3", "scRR1"),
        col=rgb(0,0,0,1/4),
        main = paste0("informative protein-coding found in RamDA,\n but not in scRR \n", dim(data_subset_gencode_ramda_protein)[1], "genes"))

boxplot(log10(data_subset_gencode_both_protein$ramda_status.mean+1),
        log10(data_subset_gencode_both_protein$scRR3_status.mean+1),
        log10(data_subset_gencode_both_protein$scRR1_status.mean+1),  ylim = c(0,5), names = c("RamDA", "scRR3", "scRR1"),
        col=rgb(0,0,0,1/4),
        main = paste0("informative protein-coding found in,\n both RamDA and scRR \n", dim(data_subset_gencode_both_protein)[1], "genes"))
#
