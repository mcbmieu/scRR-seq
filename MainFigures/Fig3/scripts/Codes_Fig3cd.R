# =============================================
# Fig. 3: DMA using S-phase progression markers
# =============================================

# Load required libraries
library(destiny)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Fig. 3c: DMA - RPE1 & CBMS1 using respective markers ---
# Load RNA and repliscore data
rpe1_df <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/04_scRR-seq-RNA_RPE1_wholeS_genes_rsem_selectedsamples_TPM_hg38.txt")
rpe1_repli <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/04_scRR-seq-DNA_RPE1_wholeS_pctreplicationscore_selectedsamples_all_hg38.txt")
cbms1_df <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/03_scRR-seq-RNA_CBMS1_mESC_genes_rsem_selectedsamples_TPM_mm10.txt")
cbms1_repli <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/03_scRR-seq-DNA_CBMS1_mESC_pctreplicationscore_selectedsamples_all_mm10.txt")

# Load S-phase progression marker genes
rpe1_marker <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/rpe1_dtClust_name.txt")
cbms1_marker <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig11/data/CBMS1_dtClust_name.txt")

## DMA for RPE1 with its own markers
rownames(rpe1_df) <- rpe1_df$gene.id
rpe1_matEx <- rpe1_df[2:ncol(rpe1_df)]
rpe1_genelist <- rpe1_marker$gene_id
rpe1_scrna_tpm_cc <- rpe1_matEx[rpe1_genelist,]
rpe1_scrna_tpm_cc_log = log10(rpe1_scrna_tpm_cc+1)

# Make DMA
set.seed(6)
rpe1_scrna_tpm_dm <- DiffusionMap(t(rpe1_scrna_tpm_cc_log))
rpe1_dm_plot <- plot(rpe1_scrna_tpm_dm, 1:2)
rpe1_scrna_tpm_dm_tmp <- data.frame(Sample = rownames(rpe1_dm_plot$data),
                               DC1 = rpe1_dm_plot$data$DC1,
                               DC2 = rpe1_dm_plot$data$DC2)
rpe1_scrna_tpm_dm_tmp_info <- merge(rpe1_scrna_tpm_dm_tmp, rpe1_repli, by = "Sample")

# Remove G1 cells to avoid complication
rpe1_scrna_tpm_dm_tmp_info <- rpe1_scrna_tpm_dm_tmp_info[!grepl("G0",rpe1_scrna_tpm_dm_tmp_info$Sample),]

# Plot
ggplot(rpe1_scrna_tpm_dm_tmp_info, aes(x = DC1, y = DC2, col = pct_repliscore)) +
  geom_point() + #scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() +
  scale_color_gradientn(limits = c(0,100), colours = rainbow(6),
                        breaks = seq(0,100,by = 20))

#
write.table(rpe1_scrna_tpm_dm_tmp_info, "/Volumes/scRR-seq_GEO/source_data/source_data_fig3/fig3c_DMA_RPE1.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)


## DMA for CBMS1 with its own markers
rownames(cbms1_df) <- cbms1_df$Gene.id
cbms1_matEx <- cbms1_df[2:ncol(cbms1_df)]
cbms1_genelist <- cbms1_marker$gene_id
cbms1_scrna_tpm_cc <- cbms1_matEx[cbms1_genelist,]
cbms1_scrna_tpm_cc_log = log10(cbms1_scrna_tpm_cc+1)

# Make DMA
set.seed(6)
cbms1_scrna_tpm_dm <- DiffusionMap(t(cbms1_scrna_tpm_cc_log))
cbms1_dm_plot <- plot(cbms1_scrna_tpm_dm, 1:2)
cbms1_scrna_tpm_dm_tmp <- data.frame(Sample = rownames(cbms1_dm_plot$data),
                               DC1 = cbms1_dm_plot$data$DC1,
                               DC2 = cbms1_dm_plot$data$DC2)
cbms1_scrna_tpm_dm_tmp_info <- merge(cbms1_scrna_tpm_dm_tmp, cbms1_repli, by = "Sample")

# Plot
ggplot(cbms1_scrna_tpm_dm_tmp_info, aes(x = DC1, y = DC2, col = pct_repliscore)) +
  geom_point() + #scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() +
  scale_color_gradientn(limits = c(0,100), colours = rainbow(6),
                        breaks = seq(0,100,by = 20))

#
#
write.table(cbms1_scrna_tpm_dm_tmp_info, "/Volumes/scRR-seq_GEO/source_data/source_data_fig3/fig3c_DMA_CBMS1.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)


# --- Reciprocal Test: Cross-species marker application ---

# Load gene ID mappings
human_asset <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.human.gencode.v38.primary_assembly.annotation.txt")
mouse_asset <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.mouse.gencode.vM25.primary_assembly.annotation.txt")

# Convert marker names
rpe1_marker_converted <- str_to_sentence(rpe1_marker$gene_name)
cbms1_marker_converted <- str_to_upper(cbms1_marker$gene_name)

# Get gene_id where gene_name matches the list
rpe1_marker_converted_mouseid <- mouse_asset[mouse_asset$gene_name %in% rpe1_marker_converted, ]
cbms1_marker_converted_humanid <- human_asset[human_asset$gene_name %in% cbms1_marker_converted, ]

## RPE1 using CBMS1 markers
rownames(rpe1_df) <- rpe1_df$gene.id
rpe1_matEx <- rpe1_df[2:ncol(rpe1_df)]
rpe1_scrna_tpm_cc <- rpe1_matEx[cbms1_marker_converted_humanid$Geneid,]
rpe1_scrna_tpm_cc_log = log10(rpe1_scrna_tpm_cc+1)

# Make DMA
library(destiny)
set.seed(3)
rpe1_scrna_tpm_dm <- DiffusionMap(t(rpe1_scrna_tpm_cc_log))
rpe1_dm_plot <- plot(rpe1_scrna_tpm_dm, 1:2)
rpe1_scrna_tpm_dm_tmp <- data.frame(Sample = rownames(rpe1_dm_plot$data),
                                    DC1 = rpe1_dm_plot$data$DC1,
                                    DC2 = rpe1_dm_plot$data$DC2)
rpe1_scrna_tpm_dm_tmp_info <- merge(rpe1_scrna_tpm_dm_tmp, rpe1_repli, by = "Sample")

# Remove G1 cells to avoid complication
rpe1_scrna_tpm_dm_tmp_info <- rpe1_scrna_tpm_dm_tmp_info[!grepl("G0",rpe1_scrna_tpm_dm_tmp_info$Sample),]
dim(rpe1_scrna_tpm_dm_tmp_info) #[1] 47 4

# Plot
ggplot(rpe1_scrna_tpm_dm_tmp_info, aes(x = DC1, y = DC2, col = pct_repliscore)) +
  geom_point() + #scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() +
  scale_color_gradientn(limits = c(0,100), colours = rainbow(6),
                        breaks = seq(0,100,by = 20))


#
write.table(rpe1_scrna_tpm_dm_tmp_info, "/Volumes/scRR-seq_GEO/source_data/source_data_fig3/fig3c_DMA_RPE1_reciprocal.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)


## CBMS1 using RPE1 markers
rownames(cbms1_df) <- cbms1_df$Gene.id
cbms1_matEx <- cbms1_df[2:ncol(cbms1_df)]
cbms1_scrna_tpm_cc <- cbms1_matEx[rpe1_marker_converted_mouseid$Geneid,]
cbms1_scrna_tpm_cc_log = log10(cbms1_scrna_tpm_cc+1)

# Make DMA
set.seed(10)
cbms1_scrna_tpm_dm <- DiffusionMap(t(cbms1_scrna_tpm_cc_log))
cbms1_dm_plot <- plot(cbms1_scrna_tpm_dm, 1:2)
cbms1_scrna_tpm_dm_tmp <- data.frame(Sample = rownames(cbms1_dm_plot$data),
                                     DC1 = cbms1_dm_plot$data$DC1,
                                     DC2 = cbms1_dm_plot$data$DC2)
cbms1_scrna_tpm_dm_tmp_info <- merge(cbms1_scrna_tpm_dm_tmp, cbms1_repli, by = "Sample")

# Plot
ggplot(cbms1_scrna_tpm_dm_tmp_info, aes(x = DC1, y = DC2, col = pct_repliscore)) +
  geom_point() + #scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() +
  scale_color_gradientn(limits = c(0,100), colours = rainbow(6),
                        breaks = seq(0,100,by = 20))

#
write.table(cbms1_scrna_tpm_dm_tmp_info, "/Volumes/scRR-seq_GEO/source_data/source_data_fig3/fig3c_DMA_CBMS1_reciprocal.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)




# --- Fig. 3d: DMA - CBMS1 markers applied to Hayashi et al. 2018 scRNA-seq ---
# Load scRNA-seq data (Hayashi et al., 2018)
mESC_ramda = read.delim("/Volumes/scRR-seq_GEO/codes/Fig3/data/merged_rsem_genes_TPM_Hayashi2018.txt")
rownames(mESC_ramda) <- mESC_ramda$gene_id
mESC_ramda[c(1,2)] <- NULL # remove gene_id and gene_name

# Apply CBMS1 markers
matEx <- as.matrix(mESC_ramda)
scrna_tpm_cc <- matEx[cbms1_genelist,]
scrna_tpm_cc_log = log10(scrna_tpm_cc+1)

# Make DMA
set.seed(2)
scrna_tpm_dm <- DiffusionMap(t(scrna_tpm_cc_log))
dm_plot <- plot(scrna_tpm_dm, 1:2)
scrna_tpm_dm_tmp <- data.frame(Sample = rownames(dm_plot$data),
                                     DC1 = dm_plot$data$DC1,
                                     DC2 = dm_plot$data$DC2)
# Use Plk1 as a marker
# ENSMUSG00000030867.7
scrna_tpm_dm_tmp$Plk1 <- scrna_tpm_cc_log["ENSMUSG00000030867.7",]

# Plot
ggplot(scrna_tpm_dm_tmp, aes(x = DC1, y = DC2, col = Plk1)) +
  geom_point() + #scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic() +
  scale_color_gradientn(
    colours = colorRampPalette(c("yellow", "red"))(10))

#
write.table(cbms1_scrna_tpm_dm_tmp_info, "/Volumes/scRR-seq_GEO/source_data/source_data_fig3/fig3d_DMA_mESC_Hayashi2018.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)


