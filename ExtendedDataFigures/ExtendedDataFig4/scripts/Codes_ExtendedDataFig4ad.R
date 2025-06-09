# ==========================================================================
# Extended Data Fig. 4a: Sample correlation of mouse 8-cell embryo samples
# ==========================================================================

# Load libraries
library("limma")
library("edgeR")
library("data.table")
library("gplots")

# --- Extended Data Fig. 4a: Sample correlation of mouse 8-cell embryo samples ---

# Load count data
data <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig4/data/merged_featureCounts_allgene_8cellembryos_scRNAsvsBulk.txt")

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
          dendrogram="col", margin=c(9, 9), cexRow=0.5, cexCol=0.5)



# --- Extended Data Fig. 4d: Detection of histone and non-poly(A) genes in mouse 8-cell embryo samples ---

# Load histone gene data (TPM)
histone_df <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig4/data/merged_featureCounts_allgene_TPM_8cellembryos_scRNAsvsBulk_histone.txt")
rownames(histone_df) <- histone_df$gene_name

# Log10 transformation
histone_df[1:3] <- NULL
histone_df_log10 <- log10(histone_df+1)

# Find median for each gene and sort genes according to median
histone_df_log10$median <- apply(histone_df_log10, 1, median)
histone_df_log10_sort <- histone_df_log10[order(histone_df_log10$median, decreasing = F),]
histone_df_log10_sort_df <- histone_df_log10_sort
histone_df_log10_sort_df$median <- NULL

heatmap(as.matrix(histone_df_log10_sort_df), 
        col = colorRampPalette(c('black','red','yellow','white'))(100),  # Define color palette
        scale = "none",  # Do not scale data
        Rowv = NA,  # Do not reorder rows
        Colv = NA,  # Do not reorder columns
        main = "Sensitivity for detecting histone genes", # Title of the plot
        breaks = seq(0, 4, length.out = 101), 
        cexRow = 0.5, cexCol = 0.3)

# Load non-poly(A) gene data (TPM)
data <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig4/data/merged_featureCounts_allgene_TPM_8cellembryos_scRNAsvsBulk.txt")

# Subset data for known non-poly(A) genes in mice
#ENSMUSG00000092274.3 # Neat1
#ENSMUSG00000115420.1 # Rmrp
#ENSMUSG00000092837.1 # Rpph1
#ENSMUSG00000092341.3 # Malat1
#ENSMUSG00000065037.1 # Rn7sk
#ENSMUSG00000099021.1 # Rn7s1


nonpolA <- c("ENSMUSG00000092274.3","ENSMUSG00000115420.1","ENSMUSG00000092837.1","ENSMUSG00000092341.3", "ENSMUSG00000065037.1", "ENSMUSG00000099021.1")
data_nonpolA <- data[nonpolA,]
rownames(data_nonpolA) <- c("Neat1","Rmrp","Rpph1","Malat1","Rn7sk","Rn7s1")

# Log10 transformation
data_nonpolA_log10 <- log10(data_nonpolA+1)

heatmap(as.matrix(data_nonpolA_log10), 
        #col = color_palette,
        #col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color palette
        col = colorRampPalette(c('black','red','yellow','white'))(100),  # Define color palette
        #col = heat.colors(20),
        scale = "none",  # Do not scale data
        Rowv = NA,  # Do not reorder rows
        Colv = NA,  # Do not reorder columns
        main = "Sensitivity for detecting non-poly(A) genes",
        breaks = seq(0, 4, length.out = 101), # Title of the plot
        cexRow = 0.5, cexCol = 0.3)

#


