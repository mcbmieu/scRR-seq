# ====================================================================
# Supplementary Fig. 5: Detected genes in RPE1 midS samples (continue)
# ====================================================================

# Load libraries
options(scipen = 999)
library(ggplot2)

# --- Supplementary Fig. 5de ---

# Load data
# Load lncATLAS_all_data_RCI (downloaded from https://lncatlas.crg.eu/ [Get Raw Data > RCI]; then converted to tabdelim format)
localization <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig5/data/2025-04-22_lncATLAS_all_data_RCI_tab.txt")

# Get only CNRCI data type
localization_CN <- localization[localization$Data.Type == "CNRCI",]

# List of genes found only in RamDA-seq or both
onlyramda <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig5/data/RPE1_06_Genes_detected_onlyinRamDA.txt")
both <- read.delim("/Volumes/scRR-seq_GEO/codes/SupplementaryFig5/data/RPE1_06_Genes_detected_Both.txt")

# Get main gene_id
onlyramda$ENSEMBL.ID <- gsub("\\..*", "", onlyramda$Geneid)
both$ENSEMBL.ID <- gsub("\\..*", "", both$Geneid)

# merge with localization_CN
onlyramda_CN <- merge(localization_CN, onlyramda, by = "ENSEMBL.ID")
both_CN <- merge(localization_CN, both, by = "ENSEMBL.ID")

# Combine two datasets
onlyramda_CN$Group <- "only_ramda"
both_CN$Group <- "both"

df_combine <- rbind(onlyramda_CN, both_CN)

# Plot
ggplot(df_combine, aes(Data.Source, Value, fill = Group))+
  #facet_wrap(~gene_type)+
  facet_wrap(Biotype~Group, nrow = 1)+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  ylim(-7, 7) +
  labs(
    x = "Cell lines/Cell types",
    y = "CN RCI",
    #fill = "Gene Type",
    title = "Cytoplasmic/Nuclear Localisation: RCI (all cell types)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(df_combine, aes(Group, Value, fill = Group))+
  #facet_wrap(~gene_type)+
  facet_wrap(~Biotype, nrow = 1)+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  ylim(-7, 7) +
  labs(
    x = "Groups",
    y = "CN RCI",
    #fill = "Gene Type",
    title = "Cytoplasmic/Nuclear Localisation: RCI (all cell types)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#


