# =======================================
# Fig. 5 (Part 2): Get Unspliced RNA Data
# =======================================

# -------- Load gene and exon-level counts --------
gene_count <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_merged_featureCounts_allgene_gene.txt")
exon_count <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_merged_featureCounts_allgene_exon.txt")

# Keep only common genes
common_gene_id <- intersect(rownames(gene_count), rownames(exon_count))
gene_count_df <- gene_count[common_gene_id, ]
exon_count_df <- exon_count[common_gene_id, ]

# Subtract exon counts from gene counts (Stuart Lee et al., NAR Genomics 2020)
subtract_df <- gene_count_df - exon_count_df
subtract_df[subtract_df < 0] <- 0  # Set negative values to zero

# -------- Add coordinates --------
asset <- read.delim("/Volumes/scRR-seq_GEO/codes/asset/asset.human.gencode.v38.primary_assembly.annotation.featurecounts.txt")
subtract_df$Geneid <- rownames(subtract_df)
exon_count_df$Geneid <- rownames(exon_count_df)

# Merge coordinates with unspliced/spliced counts
pos_unspliced_count <- merge(asset, subtract_df, by = "Geneid", sort = FALSE)
pos_spliced_count   <- merge(asset, exon_count_df, by = "Geneid", sort = FALSE)

# Export
write.table(pos_unspliced_count, "/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_unspliced.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(pos_spliced_count, "/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_spliced.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# ===============================
# RT Data Integration by Gene
# ===============================

# -------- Load RT binarized data and prepare gene coordinates --------
rt_binary <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_binarized_RT.txt")

geneid_coord <- pos_unspliced_count[c(2, 3, 4, 1)]  # chr, start, end, Geneid
colnames(geneid_coord) <- c("chr", "start", "end", "Geneid")

library(GenomicRanges)
rt_binary_gr     <- GRanges(rt_binary)
geneid_coord_gr  <- GRanges(geneid_coord)

# -------- Overlap and merge RT info to genes --------
hits <- findOverlaps(rt_binary_gr, geneid_coord_gr)

result <- cbind(
  as.data.frame(rt_binary_gr[queryHits(hits)]),
  Geneid = as.data.frame(geneid_coord_gr[subjectHits(hits)])$Geneid
)

# -------- Compute average RT signal per gene --------
library(dplyr)

inputfile <- result
output <- inputfile %>%
  group_by(Geneid) %>%
  summarize() %>%
  as.data.frame()

# Loop over columns with sample data (exclude metadata columns)
for (i in 6:(ncol(inputfile) - 1)) {
  df <- data.frame(sample = inputfile[[i]], Geneid = inputfile$Geneid)
  
  df_avr <- df %>%
    group_by(Geneid) %>%
    summarize(avr = mean(sample, na.rm = TRUE)) %>%
    as.data.frame()
  
  # Binarize
  df_avr$avr <- ifelse(df_avr$avr < 0, -1,
                       ifelse(df_avr$avr > 0, 1, 0))
  
  # Rename to match original column
  colnames(df_avr)[2] <- colnames(inputfile)[i]
  
  # Combine if Geneid order matches
  if (all(output$Geneid == df_avr$Geneid)) {
    output <- cbind(output, df_avr[2])
  } else {
    stop("Geneid list does not match; cannot cbind.")
  }
}

# -------- Export averaged RT per gene --------
write.table(output, "/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_avrRTperGene2.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
