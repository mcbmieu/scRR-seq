# ===============================================================
# Extended Data Fig. 3b: New S-phase progression markers for RPE1
# ===============================================================

# Load Required Libraries
library(data.table)
library(dplyr)
library(magrittr)
library(dtplyr)
library(mgcv)
library(flashClust)
library(ggplot2)

# Load scRR-RNA data
rpe1 = read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-RNA_RPE1_wholeS_genes_rsem_selectedsamples_TPM_hg38.txt")
rownames(rpe1) <- rpe1$gene.id
rpe1$gene.id <- NULL

# Remove G1 cells
rpe1_rmG1 <- rpe1[,!grepl("G0",colnames(rpe1))]

# Load scRR-RNA data's annotation
cellanno <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-DNA_RPE1_wholeS_pctreplicationscore_selectedsamples_all_hg38.txt")

# Remove G1 cells to avoid complication
cellanno_all <- cellanno[!grepl("G0",cellanno$Sample),]

# Log10 transformation
rpe1_log10 <- log10(rpe1_rmG1+1)

# Sort data according to repliscore
sorted_repli <- cellanno_all$Sample[order(cellanno_all$pct_repliscore)]
rpe1_log10_sort <- rpe1_log10[,sorted_repli]

matEx <- as.matrix(rpe1_log10_sort)
repliscore <- cellanno_all$pct_repliscore[order(cellanno_all$pct_repliscore)]

####################################
# Modified from https://github.com/yuifu/tutorial-RamDA-paper-fugures
# Fiting with GAM
## 00_GAM_fitting ----
# Load libraries
library(data.table); library(dplyr); library(magrittr); library(dtplyr); library(mgcv)

####################################
# Fit generalized additive model (GAM) for each Tx
idx = 1:nrow(matEx)

listFit0 = lapply(idx, function(i){
  y = unlist(matEx[i,])
  gam(y~1, family=gaussian(link=identity))
})

listFit = lapply(idx, function(i){
  y = unlist(matEx[i,])
  gam(y~s(repliscore), family=gaussian(link=identity))
})

head(listFit0)
head(listFit)

#####################################
# Define functions
extractPvalue = function (fit) {
  (summary(fit)$s.table)[4]
}
extractFittedValue = function (fit) {
  fit$fitted.values
}


####################################
# Calculate false discovery rate (FDR)
resP = sapply(listFit, extractPvalue)
fdr = p.adjust(resP, method = "BH")

####################################
# Calculate AIC and BIC
# https://www.r-bloggers.com/2021/10/model-selection-in-r-aic-vs-bic/
# (the Akaike Information Criterion and the Bayesian Information Criterion described in Fabozzi et al. (2014)
aic = sapply(listFit, AIC)
bic = sapply(listFit, BIC)

aic = sapply(listFit, AIC)
bic = sapply(listFit, BIC)

aic0 = sapply(listFit0, AIC)
bic0 = sapply(listFit0, BIC)

####################################
# Extract fitted expression levels
fittedValues = sapply(listFit, extractFittedValue) %>% t
dim(fittedValues) #[1] 55487    84
colnames(fittedValues) = colnames(matEx)
head(fittedValues)

####################################
# Save fitting results (p-value, FDR, AIC, BIC)
outdir="/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data"
dir.create(outdir, recursive = T)

dtRes = data.frame(gene_id = rownames(rpe1_log10_sort))
dtRes$pvalue = resP
dtRes$fdr = fdr
dtRes$aic = aic
dtRes$bic = bic
dtRes$aic0 = aic0
dtRes$bic0 = bic0
write.table(dtRes, paste0(outdir,"/rpe1_fit_gam.txt"), row.names = F, quote = F, sep = "\t")

# Save fitted values
dtRes2 = data.frame(gene_id = rownames(rpe1_log10_sort))
dtRes2 = cbind(dtRes2, fittedValues)
write.table(dtRes2, gzfile(paste0(outdir,"/rpe1_fittedValues.txt.gz")), row.names = F, quote = F, sep = "\t")

###
### 01_Clustering_expression ---
# Define File paths
dtRes <- read.table("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/rpe1_fit_gam.txt", header = T)
dtRes2 <- read.table("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/rpe1_fittedValues.txt.gz", header = T)

# Load libraries
#install.packages("flashClust")
library(data.table); library(dplyr); library(magrittr); library(dtplyr); library(flashClust); library(ggplot2)

# Define samples based on repliscores
selCol = sorted_repli

# Load GAM fitting data
dtFdr = dtRes

# Load fitted expression value data
dtEx = dtRes2

# Select genes used for clustering
## Select genes whose AIC value was lower than that of intercept model (aic0)
### Also, filter low-expressed genes.
#### Note that dtEx value is already log-transformed [log10(TPM+1)]

threshold_fdr = 0.05
threshold_detection = log10(1+1)
nCellRatio = 0.2
nCluster = 3

#  Filtered low reads
nCell = round(dim(cellanno_all)[1] * nCellRatio)
print(sprintf("nCell: %d", nCell))
selRowPre_df = rpe1_log10_sort[rowSums(rpe1_log10_sort > threshold_detection) > nCell,]
selRowPre = rownames(selRowPre_df)

# Select significant genes
selRow = dtFdr$gene_id[dtFdr$fdr < threshold_fdr & dtFdr$gene_id %in% selRowPre & dtFdr$aic < dtFdr$aic0]
print(length(selRow)) #[1] 52

# Get fitted values for significant genes
dtEx_df <- data.frame(dtEx)
rownames(dtEx_df) <- dtEx_df$gene_id
dtEx_df$gene_id <- NULL

matEx_fit = as.matrix(dtEx_df[selRow, selCol])
print(dim(matEx_fit)) #[1] 52  47

# Calculate Pearson distance between genes
d = as.dist((1 - cor(t(matEx_fit)))/2)

# Perform clustering
clusters = flashClust(d, method = "ward", members=NULL)
clusterCut <- cutree(clusters, nCluster)

# Save clustering results
dtClust = data.table(gene_id = selRow, cluster = clusterCut)
dtClust %>% setkey(gene_id)
ofile = sprintf("%s/%s", outdir, "rpe1_dtClust.txt")
write.table(dtClust, ofile, row.names = F, quote = F, sep = "\t")

# Add name
asset_df <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.human.gencode.v38.primary_assembly.annotation.txt")
colnames(asset_df)[1] <- "gene_id"
dtClust_name <- merge(dtClust, asset_df, by = "gene_id")
write.table(dtClust_name, paste0(outdir, "/rpe1_dtClust_name.txt"), row.names = F, quote = F, sep = "\t")

# Save each cluster
lapply(dtClust[, unique(cluster)], function(cl){
  ofile = sprintf("%s/rpe1_gene_id.cluster_%d.txt", outdir, cl)
  write.table(dtClust[cluster == cl, .(gene_id)], ofile, row.names = F, quote = F, sep = "\t")
})

# Calculate mean and sd (standard deviation) of expression values for each cluster
colMeansWithScaling = sapply(unique(clusterCut), function(cl){
  colMeans(t(scale(t(matEx_fit[clusterCut==cl,]))))
})
colSdWithScaling = sapply(unique(clusterCut), function(cl){
  apply(t(scale(t(matEx_fit[clusterCut==cl,]))), 2, sd) %>% t
}) 

colnames(colMeansWithScaling) = sprintf("cluster_%d", unique(clusterCut))
colnames(colSdWithScaling) = sprintf("cluster_%d", unique(clusterCut))
rownames(colSdWithScaling) = rownames(colMeansWithScaling)

dtCluterMeanWithScaling = as.data.table(colMeansWithScaling, keep.rownames = TRUE)
dtCluterMeanWithScaling %>% setnames("rn", "Sample")

dtTmp = as.data.table(colSdWithScaling, keep.rownames = TRUE)
dtTmp %>% setnames("rn", "Sample")
dtTmp = dtTmp  %>% melt(id.vars = c("Sample"), value.name = "sd")


# Create annotation
annot <- cellanno_all
dtClustRowMeans = merge(annot, dtCluterMeanWithScaling, by = "Sample")
dtClustRowMeans = dtClustRowMeans %>% melt(id.vars = c("Sample", "pct_repliscore"))
dtClustRowMeans = merge(dtClustRowMeans, dtTmp, by=c("Sample", "variable"))

ofile = sprintf("%s/%s", outdir, "rpe1_cluster_mean_withScaling_withSd.txt")
write.table(dtClustRowMeans, ofile, row.names = F, quote = F, sep = "\t")

# Plot mean and sd (standard deviation) of expression values for each cluster along pseudotime
str(dtClustRowMeans)

dtClustRowMeans$value <- as.numeric(dtClustRowMeans$value)

g = ggplot(dtClustRowMeans, aes(pct_repliscore, value))
g = g + facet_wrap(~variable, ncol=1, scale = "free_y")
g = g + geom_ribbon(aes(ymax = value + sd, ymin = value - sd), alpha = 0.8, fill = "skyblue")
g = g + geom_line()
g = g + theme_bw()
g = g + theme(panel.grid.minor = element_blank())
print(g)

# Plot fitted value and gene expression in each cell for specific gene
plot(repliscore, matEx[rownames(matEx) == "ENSG00000065328.17",],
     ylab = "log10(TPM+1)", xlab="repliscore")
lines(repliscore, dtEx_df[rownames(dtEx_df) == "ENSG00000065328.17",], col = "red")

#
