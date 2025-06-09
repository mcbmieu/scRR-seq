# ============================================================
# Fig. 3: Cell Cycle Assignment using Seurat on scRR-seq data
# ============================================================

# --- Fig. 3e: Seurat Cell Cycle Scoring on RPE1 scRR-seq data ---

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)

# Load data
rpe1_df <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-RNA_RPE1_wholeS_genes_rsem_selectedsamples_TPM_hg38.txt")
rpe1_repli <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/04_scRR-seq-DNA_RPE1_wholeS_pctreplicationscore_selectedsamples_all_hg38.txt")

# Remove G1 cells
rpe1_df <- rpe1_df[,!grepl("G0",colnames(rpe1_df))]
rownames(rpe1_df) <- rpe1_df$gene.id
rpe1_df$gene.id <- NULL
rpe1_repli <- rpe1_repli[!grepl("G0",rpe1_repli$Sample),]

# Prepare Seurat object
info_data <- data.frame(Sample = rpe1_repli$Sample, repliscore = rpe1_repli$pct_repliscore)
rownames(info_data) <- info_data$Sample

scrna_tpm <- rpe1_df
all.genes <- rownames(scrna_tpm)
share_df <- CreateSeuratObject(counts = scrna_tpm, meta.data = info_data)
share_df <- NormalizeData(share_df, normalization.method = "LogNormalize",
                          scale.factor = 10000)
share_df <- FindVariableFeatures(share_df, 
                                 selection.method = "vst",
                                 nfeatures = 2000, 
                                 verbose = FALSE)
share_df <- ScaleData(share_df, features = all.genes)

# Perform PCA
share_df <- RunPCA(share_df, npcs = 20)
DimPlot(share_df, reduction = "pca")

## Cell-cycle assignment
# https://satijalab.org/seurat/articles/cell_cycle_vignette
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Prepare asset
human_genename <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.human.gencode.v38.primary_assembly.annotation.txt")

# Convert genename ot geneid
s.genes.id <- human_genename[human_genename$gene_name %in% s.genes, ,drop = F]
g2m.genes.id <- human_genename[human_genename$gene_name %in% g2m.genes, ,drop = F]

# CellCycleScoring
share_df <- CellCycleScoring(share_df, s.features = s.genes.id$Geneid, g2m.features = g2m.genes.id$Geneid, set.ident = TRUE)
cc.seurat <- share_df@meta.data

# Plot: replication score and Seurat cell-cycle phase
cc.seurat_sort <- cc.seurat[order(cc.seurat$repliscore),]
cc.seurat_sort$num <- 1:nrow(cc.seurat)

p1 <- ggplot(cc.seurat_sort, aes(num, y = 1, fill = repliscore))+
  geom_tile(color="gray") + theme_void()+
  scale_fill_gradientn(limits = c(0,100), colours = rainbow(6),
                       breaks = seq(0,100,by = 20))

p2 <-ggplot(cc.seurat_sort, aes(num, y = 1, fill = Phase))+
  geom_tile(color="gray") + theme_void()

combined_plot <- p1 / p2
print(combined_plot)
# # Or use p1 | p2 for side-by-side





# --- Fig. 3f: Seurat Cell Cycle Scoring on CBMS1 scRR-seq data ---
# Load data
cbms1_df <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/03_scRR-seq-RNA_CBMS1_mESC_genes_rsem_selectedsamples_TPM_mm10.txt")
cbms1_repli <- read.delim("/Volumes/scRR-seq_GEO/codes/ExtendedDataFig6/data/03_scRR-seq-DNA_CBMS1_mESC_pctreplicationscore_selectedsamples_all_mm10.txt")

# Prepare data
rownames(cbms1_df) <- cbms1_df$Gene.id
cbms1_df$Gene.id <- NULL
cbms1_repli <- cbms1_repli[!grepl("SRR",cbms1_repli$Sample),]


# Prepare Seurat object
info_data <- data.frame(Sample = cbms1_repli$Sample, repliscore = cbms1_repli$pct_repliscore)
rownames(info_data) <- info_data$Sample

scrna_tpm <- cbms1_df
all.genes <- rownames(scrna_tpm)

share_df <- CreateSeuratObject(counts = scrna_tpm, meta.data = info_data)
share_df <- NormalizeData(share_df, normalization.method = "LogNormalize",
                          scale.factor = 10000)
share_df <- FindVariableFeatures(share_df, 
                                 selection.method = "vst",
                                 nfeatures = 2000, 
                                 verbose = FALSE)
share_df <- ScaleData(share_df, features = all.genes)

# Perform PCA
share_df <- RunPCA(share_df, npcs = 20)
DimPlot(share_df, reduction = "pca")

## Cell-cycle assignment
# https://satijalab.org/seurat/articles/cell_cycle_vignette
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Prepare asset
mouse_genename <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.mouse.gencode.vM25.primary_assembly.annotation.txt")

# Convert genename ot geneid
library(stringr)
s.genes.sentence <- str_to_sentence(s.genes)
g2m.genes.sentence <- str_to_sentence(g2m.genes)

s.genes.id <- mouse_genename[mouse_genename$gene_name %in% s.genes.sentence, ,drop = F]
g2m.genes.id <- mouse_genename[mouse_genename$gene_name %in% g2m.genes.sentence, ,drop = F]

# CellCycleScoring
share_df <- CellCycleScoring(share_df, s.features = s.genes.id$Geneid, g2m.features = g2m.genes.id$Geneid, set.ident = TRUE)
head(share_df@meta.data)
cc.seurat <- share_df@meta.data

# Plot: replication score and Seurat cell-cycle phase
cc.seurat_sort <- cc.seurat[order(cc.seurat$repliscore),]
cc.seurat_sort$num <- 1:nrow(cc.seurat)

p1 <- ggplot(cc.seurat_sort, aes(num, y = 1, fill = repliscore))+
  geom_tile(color="gray") + theme_void()+
  scale_fill_gradientn(limits = c(0,100), colours = rainbow(6),
                       breaks = seq(0,100,by = 20))

p2 <-ggplot(cc.seurat_sort, aes(num, y = 1, fill = Phase))+
  geom_tile(color="gray") + theme_void()

combined_plot <- p1 / p2
print(combined_plot)
# # Or use p1 | p2 for side-by-side



#



