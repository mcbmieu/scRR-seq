# =================================================
# Fig. 6 (Part 2): Plot CNVs as heatmap genome-wide
# =================================================

# Run in R

# -------- Load gene and exon-level counts --------
library(AneuFinder) #v.1.2.1

Aneufinder(inputfolder='~/scRR_G1_IMR90_hg38/bam', # Use bam files from scRR-seq-DNA analysis
           outputfolder='~/scRR_G1_IMR90_hg38/results',
           numCPU=4,
           method=c('dnacopy'),
           correction.method = 'mappability',
           blacklist = '~/programs/scRepliseq_Pipeline_RPv1/blacklist/hg38-blacklist.v2.mod.bed',
           mappability.reference = '~/scRR_G1_IMR90_hg38/bam_G1/merged_control_1.hg38.clean_srt_markdup.bam',
           binsizes = c(5e4,1e5,2e5,5e5,1e6),
           eps = 0.1,
           max.time = 60,
           max.iter = 3000,
           chromosomes = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10", "chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                           "chr18","chr19","chr20","chr21","chr22","chrX"))

# Plot CNV for a whole genome
results <- "~/scRR_G1_IMR90_hg38/results/MODELS/method-dnacopy"
files <- list.files(results, full.names=TRUE, pattern="binsize_1e\\+05_") #
cl <- clusterByQuality(files, measures=c('spikiness'))
plot(cl$Mclust, what='classification')
selected.files <- unlist(cl$classification[1:2]) # Or list manually

heatmapGenomewide(selected.files, ylabels = NULL, classes = NULL,
                  reorder.by.class = T, # 'F' if not need
                  classes.color = NULL, file = NULL,
                  cluster = T, # 'F' if not need
                  plot.SCE = TRUE, hotspots = NULL)

heatmapGenomewide(manual2, ylabels = NULL, classes = NULL,
                  reorder.by.class = F, classes.color = NULL, file = NULL,
                  cluster = F, plot.SCE = TRUE, hotspots = NULL)


