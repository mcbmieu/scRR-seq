# ======================================================= 
# Fig. 4f (part 1): Allelic gene expression ratio vs. RT
# =======================================================

# --- RT differences between CBA and MSM ----

# Load required libraries
library(dplyr)
library(ggplot2)


# Load and summarize RT binarized data (CBA)
cba <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/03_scRR-seq-DNA_CBMS1_mESC_allelic_binarized_selectedsamples_CBA_sorted_mm10.txt", header = T)
cba_df <- cba[1:3]
cba_df$cba_none <- rowSums(cba[4:ncol(cba)] == 0)
cba_df$cba_rep <- rowSums(cba[4:ncol(cba)] == 1)
cba_df$cba_unrep <- rowSums(cba[4:ncol(cba)] == -1)

# Label informative bins
rowlength = ncol(cba)-3
cba_df <- cba_df %>%
  mutate(cba_informative = case_when(
    cba_none/rowlength > 0.2  ~ "no",
    TRUE ~ "yes"
  ))



# Load and summarize RT binarized data (MSM)
msm <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/03_scRR-seq-DNA_CBMS1_mESC_allelic_binarized_selectedsamples_MSM_sorted_mm10.txt", header = T)
msm_df <- msm[1:3]
msm_df$msm_none <- rowSums(msm[4:ncol(msm)] == 0)
msm_df$msm_rep <- rowSums(msm[4:ncol(msm)] == 1)
msm_df$msm_unrep <- rowSums(msm[4:ncol(msm)] == -1)
rowlength = ncol(msm)-3

# Label informative bins
msm_df <- msm_df %>%
  mutate(msm_informative = case_when(
    msm_none/rowlength > 0.2  ~ "no",
    TRUE ~ "yes"
  ))



# Merge CBA and MSM RT summaries
merged_df <- merge(cba_df,msm_df, by = c(1:3))

# Perform Fisher's exact test to detect RT differences
fisher_p <- function(x) {
  if (x["cba_informative"] == "yes" & x["msm_informative"] == "yes") {
    # Need to convert to numeric, if not it fails
    cba_rep <- as.numeric(x["cba_rep"])
    cba_unrep <- as.numeric(x["cba_unrep"])
    msm_rep <- as.numeric(x["msm_rep"])
    msm_unrep <- as.numeric(x["msm_unrep"])
    df <- matrix(c(cba_rep, cba_unrep, msm_rep, msm_unrep), nrow = 2, byrow = TRUE)
    p <- fisher.test(df, alternative = "two.sided")
    pvalue.f <- p$p.value
  } else {
    pvalue.f <- NA
  }
  return(pvalue.f)
}

merged_df$p.fishertest <- apply(merged_df, 1, fisher_p)
merged_df$logp.fishertest <- log10(merged_df$p.fishertest)

# Export Fisher's test results
write.table(merged_df, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RT_differences_CBA_MSM_pfishertest.txt", sep = "\t", row.names = F, col.names = T, quote = F)
#



# Merge with average mid-S RT values
midS_cba <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/03_scRR-seq-DNA_CBMS1_mESC_allelic_midS_Avg_RT_CBA_mm10.bedGraph", header = F)
midS_msm <- read.delim("/Volumes/scRR-seq_GEO/codes/Fig4/data/03_scRR-seq-DNA_CBMS1_mESC_allelic_midS_Avg_RT_CBA_mm10.bedGraph", header = F)

colnames(midS_cba) <- c("chr", "start", "end", "midS_cba")
colnames(midS_msm) <- c("chr", "start", "end", "midS_msm")

midS_cba_msm <- merge(midS_cba, midS_msm, by = c(1:3))
midS_cba_msm$RTdiff <- midS_cba_msm$midS_cba-midS_cba_msm$midS_msm

midS_cba_msm_df <- merge(merged_df, midS_cba_msm, by = c(1:3))

write.table(midS_cba_msm_df, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RT_differences_CBA_MSM_midSdata.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(midS_cba_msm_df, "/Volumes/scRR-seq_GEO/codes/Fig4/data/RT_differences_CBA_MSM_midSdata_noheader.txt", sep = "\t", row.names = F, col.names = F, quote = F)
#

# Continue part 2
