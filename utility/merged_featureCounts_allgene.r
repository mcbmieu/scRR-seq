#!/usr/bin/env Rscript

# RP code for RamdaQ 240201
# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

inputdir <- args[1]
outputdir <- args[2]

# Make out dir
dir.create(outputdir, recursive = T)

# Extract data for isoform.results
inputfile <- list.files(path = inputdir, full.names = T, pattern = "allgene.featureCounts")

for (i in 1:length(inputfile)){
  out_df <- read.table(inputfile[1], header = T, skip = 1)
  out_df <- out_df[c(1,6,7)]
}

for (j in inputfile){
  df <- read.table(j, header = T)
  df <- df[8]
  sample <- gsub(".allgene.featureCounts.txt", "", basename(j))
  colnames(df) <- sample
  out_df <- cbind(out_df, df)
}

write.table(out_df, paste0(outputdir,"/merged_featureCounts_allgene.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# Printing sessioninfo to standard out
print("Done")
sessionInfo()







