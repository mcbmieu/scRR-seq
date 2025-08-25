# ====================================================================
# Supplementary Fig. 3: Correlation and PCA of RT of RPE1 midS samples
# ====================================================================


# --- Supplementary Fig. 3e: Plot a heatmap of correlations ---

# Load libraries
library(gplots)
library(ggplot2)

# Load data
inputdf <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig3/data/scRR-seq-DNA_RPE1_midS_binarized_selectedsamples_all_sorted_hg38.txt")
inputinfo <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig3/data/RPE1midS_info_repliscores.txt")

# Prepare and plot correlation matrix
cordat <- inputdf[4:ncol(inputdf)]
heatmap.2(as.matrix(cor(cordat, method="pearson")),
          key.title="Pearson's Correlation", trace="none", scale="none",
          dendrogram="row", margin=c(9, 9))


# --- Supplementary Fig. 3f: PCA ---

# Load data
inputdf <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig3/data/scRR-seq-DNA_RPE1_midS_binarized_selectedsamples_all_sorted_hg38.txt")
inputinfo <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig3/data/RPE1midS_info_repliscores.txt")

# PCA computation
cordat <- inputdf[4:ncol(inputdf)]
repli.pca = prcomp(t(cordat))
repli.pca.layout = data.frame(repli.pca$x)
repli.pca.layout$Sample = rownames(repli.pca.layout)

# Add annotation
repli.pca.layout.anno <- merge(repli.pca.layout, inputinfo, by = "Sample")

# Plot
ggplot(repli.pca.layout.anno, aes(x = PC1, y = PC2, color = exp))+
  geom_point(aes(fill = exp), alpha = 0.75) + theme_classic()
ggplot(repli.pca.layout.anno, aes(x = PC1, y = PC2))+
  geom_point(aes(fill = repliscore, col = repliscore), alpha = 0.75) + theme_classic()+
  scale_color_gradientn(values = c(0,1), colours = colorRampPalette(c("gold","blue"))(10))+
  scale_fill_gradientn(values = c(0,1), colours = colorRampPalette(c("gold","blue"))(10))

#