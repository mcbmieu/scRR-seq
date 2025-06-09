# =====================================================================
# Extended Data Fig, 2: Correlation, PCA, and tSNE of RPE1 midS samples
# =====================================================================

# --- Extended Data Figure 2d: Sample correlation by edgeR ---

# Load required libraries
library(limma)
library(edgeR)
library(data.table)
library(gplots)
library(ggplot2)
library(dplyr)
library(Rtsne)

# Load count column from all files into a list of data frames
data <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig2/data/merged_featureCounts_allgene_RPE1midS.txt")

# Convert data frame to edgeR DGE object
dataDGE <- DGEList(counts=data.matrix(data))

# Normalise counts
dataNorm <- calcNormFactors(dataDGE)

# Get the log counts per million values using cpm
logcpm <- cpm(dataNorm, prior.count=2, log=TRUE)

# Calculate the euclidean distances between samples
dists = dist(t(logcpm))
str(dists)

# Plot a heatmap of correlations
heatmap.2(as.matrix(cor(logcpm, method="pearson")),
          key.title="Pearson's Correlation", trace="none", scale="none",
          dendrogram="row", margin=c(9, 9), cexRow=0.5, cexCol=0.5)



# --- Extended Data Figure 1e: PCA ---

# Load data
rsem_gene <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig1/data/scRR-seq-RNA_RPE1_midS_genes_rsem_selectedsamples_TPM_hg38.txt", header = T)

# PCA
inputdf <- rsem_gene[2:ncol(rsem_gene)]
rna.pca = prcomp(t(log10(inputdf+1)))
rna.pca.layout = data.frame(rna.pca$x)
rna.pca.layout$Sample = rownames(rna.pca.layout)

# Add annotation
rna.pca.layout <- rna.pca.layout %>% 
  mutate(exp = case_when(
    grepl("scRR3", Sample) ~ "scRR3",
    grepl("scRR1", Sample) ~ "scRR1",
    grepl("scRamDA", Sample) ~ "scRamDA"),
    batch = case_when(
      grepl("Lot1", Sample) ~ "batch1",
      grepl("Lot2", Sample) ~ "batch2",
      TRUE ~ "batch1")
  )

# Plot
ggplot(rna.pca.layout, aes(x = PC1, y = PC2, color = exp, shape = batch))+
  geom_point(aes(fill = exp), alpha = 0.75) + theme_classic()


# --- Extended Data Figure 1f: tSNE ---
set.seed(42)
tsne_out <- Rtsne(t(log10(inputdf+1)), dims = 2, perplexity = 5)
tsne_plot <- data.frame(tSNE1 = tsne_out$Y[,1], 
                        tSNE2 = tsne_out$Y[,2],
                        Sample = colnames(inputdf),
                        exp = rna.pca.layout$exp,
                        batch = rna.pca.layout$batch)

# Plot
ggplot(tsne_plot, aes(x = tSNE1, y = tSNE2, color = exp, shape = batch))+
  geom_point(aes(fill = exp), alpha = 0.75) + theme_classic()

# 
