# =======================================================
# Fig. 4f (part 3): Allelic gene expression ratio vs. RT
# =======================================================

# --- Plot results ----

# Load required libraries
library(dplyr)
library(ggplot2)

# Load intersect RT and gene expression data (see: Codes_Fig4f_part2.txt)
rt_ge <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/Intersect_allelic_RT_GE_CBMS1.txt", header = T)

# Keep only genes with a unique gene_id
unique_samples <- rt_ge %>%
  group_by(gene_id) %>%
  filter(n() == 1) %>%
  ungroup()

unique_samples_df <- as.data.frame(unique_samples)

# Keep only significantly different RT regions (Fisher's test)
rt_ge_sig <- unique_samples_df[is.na(unique_samples_df$p.fishertest) == F & unique_samples_df$p.fishertest < 0.05,]

# Classify gene allelic patterns (based on median_ratio)
rt_ge_sig <- rt_ge_sig %>%
  mutate(category_median_2 = case_when(
    ge_ratio >= 0.85 ~ "monoallelic_cba",
    ge_ratio <= 0.15 ~ "monoallelic_msm",
    ge_ratio >= 0.7 & ge_ratio < 0.85 ~ "biased_cba",
    ge_ratio >= 0.15 & ge_ratio <= 0.3 ~ "biased_msm",
    ge_ratio > 0.3 & ge_ratio < 0.7 ~ "biallelic",
    TRUE ~ "Unknown"
  ))

# Plot
ggplot(rt_ge_sig, aes(x=ge_ratio, y=RTdiff, col = category_median_2))+
  geom_boxplot() + geom_point(alpha = 0.5) +
  theme_classic() + ylab("midS_RT: CBA-MSM") + xlab("TPM: CBA/(CBA+MSM)")

#