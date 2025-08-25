
# ======================================================
# Fig. 5b (Part 2): MidS_RT and Median gene expression
# ======================================================

# Load data
# MidS-RT value of all genes
gene_RT <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/MidSRT_genes.txt")
head(gene_RT)
dim(gene_RT) #[1] 42996 8

# Sort accordint ot midS-RT
gene_RT_sort <- gene_RT[order(gene_RT$mids, decreasing = T),]

# Classify genes according to their corresponding RT
range <- seq(min(gene_RT$mids), max(gene_RT$mids), length = 9)

genes_early <- gene_RT_sort[gene_RT_sort$mids > range[7],]
genes_late <- gene_RT_sort[gene_RT_sort$mids < range[3],]
genes_mid <- gene_RT_sort[gene_RT_sort$mids > range[3] & gene_RT_sort$mids < range[7],]

gene_RT_sort <- gene_RT_sort %>% mutate(
  class = case_when(
    mids > range[7] ~ "early",
    mids < range[3] ~ "late",
    mids > range[3] & mids < range[7] ~ "mid"
  )
)

print(gene_RT_sort)

# Examine expression of genes in RT groups: Spliced
exp <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_spliced.txt")
rownames(exp) <- exp$Geneid
exp_df <- exp
exp_df[1:7] <- NULL

# log10
exp_log <- log10(exp_df +1 )
head(exp_log)

# Find median expression for genes
exp_log$median <- apply(exp_log[1:ncol(exp_df)], 1, median)

# If 50% of data has value, defined as formative
exp_log$formative <- rowSums(exp_log[1:ncol(exp_df)] > 0) > 28 
exp_log$Geneid <- rownames(exp_log)

# Include midS RT data
exp_log_df <- merge(gene_RT_sort, exp_log, by = "Geneid")
dim(exp_log_df) #[1] 42996    67
head(exp_log_df)

#
write.table(data.frame(exp_log_df$Geneid, exp_log_df$class, exp_log_df$median),
            "/Volumes/scRR-seq_GEO/source_data/source_data_fig5/fig5b_midS_RTclass_median_expression_spliced.txt",
            sep="\t", col.names = T, row.names = F, quote = F)

# Check stat
summary_df <- as.data.frame(exp_log_df[exp_log_df$formative == T,] %>%
                              group_by(class = class) %>%
                              summarize(sum =n(), med_median = median(median)))
print(summary_df)
#class  sum med_median
#1 early 7765  1.1786338
#2  late  857  0.9209018
#3   mid 2070  1.1461280


# Plot
exp_log_df$class <- factor(exp_log_df$class, levels = c("early", "mid", "late"))
ggplot(exp_log_df[exp_log_df$formative == T,], aes(x= class, y = median))+
  geom_violin() + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + ylim(0,4) +
  theme_classic() + ggtitle("Median of formative genes: spliced")





# Examine expression of genes in RT groups: Unpliced
exp <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig5/data/HAP1midS_1Mscale_featureCounts_midS_unspliced.txt")
rownames(exp) <- exp$Geneid
exp_df <- exp
exp_df[1:7] <- NULL

# log10
exp_log <- log10(exp_df +1 )
head(exp_log)

# Find median expression for genes
exp_log$median <- apply(exp_log[1:ncol(exp_df)], 1, median)

# If 50% of data has value, defined as formative
exp_log$formative <- rowSums(exp_log[1:ncol(exp_df)] > 0) > 28 
exp_log$Geneid <- rownames(exp_log)

# Include midS RT data
exp_log_df <- merge(gene_RT_sort, exp_log, by = "Geneid")
dim(exp_log_df) #[1] 42996    67
head(exp_log_df)

#
write.table(data.frame(exp_log_df$Geneid, exp_log_df$class, exp_log_df$median),
            "/Volumes/scRR-seq_GEO/source_data/source_data_fig5/fig5b_midS_RTclass_median_expression_unspliced.txt",
            sep="\t", col.names = T, row.names = F, quote = F)

# Check stat
summary_df <- as.data.frame(exp_log_df[exp_log_df$formative == T,] %>%
                              group_by(class = class) %>%
                              summarize(sum =n(), med_median = median(median)))
print(summary_df)
#class  sum med_median
#1 early 4654  0.4771213
#2  late 1099  0.2041200
#3   mid 1812  0.3010300


# Plot
exp_log_df$class <- factor(exp_log_df$class, levels = c("early", "mid", "late"))
ggplot(exp_log_df[exp_log_df$formative == T,], aes(x= class, y = median))+
  geom_violin() + 
  geom_boxplot(width = 0.05, outlier.shape = NA) + ylim(0,4) +
  theme_classic() + ggtitle("Median of formative genes: unspliced")
