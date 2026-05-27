# qrsh -l h_data=10G,h_rt=08:00:00
# module load R/3.6.1
# R

## Install DiffBind if necessary, has to use R 3.6.1 to install this version, same as my MAC version.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DiffBind")

library(DiffBind)
## Data input: sample sheet and peaksets
#setwd("/u/project/geschwind/jennybea/ATAC/NEW_ASD_parietal/DiffBind_Consensus_peaks/")
setwd("/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/4_ConsensusPeaks/Diffbind/")
Peaksets <- dba(sampleSheet="SampleSheet.csv") # took 1 min 
Peaksets # To check the compiled peaksets. 
# 38 Samples, 195373 sites in matrix (303131 total)

## Generate a correlation heatmap which gives an initial clustering of the samples using the cross-correlations of each row of the binding matrix, which contains peaks that overlap in at least two of the samples.
if (! dir.exists("plots")) {dir.create("plots")}
Peaksets$config$treatment <- "AgeGroup"
Peaksets$config$condition <- "Diagnosis"
pdf("plots/Peaksets_Correlation_Heatmap_Batch.pdf", height = 10, width = 10)
plot(Peaksets) 
dev.off()
# Observation: 
# 0) Min correlation is over 0.8, Hooray!  
# 1) The 6 pairs of replicates are always clustered together.
# 2) Tissue/Batch has an impact. The Parietal/Batch1 samples tend to form subclusters within large clusters.
# 3) AgeGroups are scattered all around.

## Occupancy analysis and overlaps
# Remove replicates
SS=read.csv("SampleSheet.csv")
idx_rm = which(SS$SampleID %in% c("B4334-1", "B4721A", "B5000B", "B5718D", "B5813B", "CQ56-2"))
SS = SS[-idx_rm,]
write.csv(SS,file="SampleSheet_rmTechRep.csv",quote=F,row.names = F)

Peaksets <- dba(sampleSheet="SampleSheet_rmTechRep.csv") # took 1 min 
Peaksets # To check the compiled peaksets. 
# 32 Samples, 189783 sites in matrix (294760 total)

olap.rate <- dba.overlap(Peaksets,mode=DBA_OLAP_RATE)
olap.rate # list how many peaks are identified in 1-38 samples
pdf("plots/Number_of_peaks_Overlap_in_Number_of_samples.pdf",height = 5,width = 10)
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
dev.off()
# The overlap rate starts to decrease in a arithmetric (as apposed to geometrically) progression at 5/32 (>15%) samples.

## Generating binding affinity matrix
Peaksets$config$treatment <- "AgeGroup"
Peaksets$config$condition <- "Diagnosis"

for (m in 3:7) {
  print(m)
  Peaksets <- dba.count(Peaksets, minOverlap = m) # Default minOverlap=2. This step takes time!
  print("Done peakset overlap.")
  save(Peaksets, file = paste0("Peaksets_minOverlap",m,".rda")) # FRiP showing fraction of reads in peak. Bryois et al requires at least 1% for each sample. Here i have 14%-43%
  pdf(paste0("plots/Peaksets_minOverlap",m,"_Correlation_Heatmap_affinity_score.pdf"), height = 10, width = 10)
  plot(Peaksets) # This generates correlation heatmap which gives an initial clustering of the samples using the cross-correlations of each row of the binding matrix, which contains peaks that overlap in at least two of the samples.
  dev.off()
}
# 2022/08/11 3:07-3:59 PM

# Observation:
# The clustering gets stable after minOverlap=5. So use minOverlap=5

## Get consensus peaks
load("Peaksets_minOverlap5.rda")
dataset = as.data.frame(Peaksets$peaks[[1]])
dim(dataset) # 128260 consensus ATAC peaks x 8 columns. Yuyan obtained a total of 151630 consensus peaks. 

write.table(dataset[,1:3], file = "Peaksets_minOverlap5.bed", quote = F, row.names = F, col.names = F, sep = "\t")
# Copy to mac

## Average peak width
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/4_ConsensusPeaks/1_DiffBind/")

dataset = read.table("Peaksets_minOverlap5.bed")
colnames(dataset) = c("CHR", "START", "END")
dataset$length = dataset$END - dataset$START
mean(dataset$length) # 1183bp, Yuyan's is 1282bp.

pdf("plots/Peaksets_minOverlap5_peaklength.pdf", height = 5, width = 8)
hist(dataset$length, breaks = 30, xlab = "Consensus ATAC peak length")
dev.off()

## Annotation using ChIPseeker
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ChIPseeker")
# Y to all install questions, n to Update all/some/none.
# Tutorial, Refer to http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# BiocManager::install("org.Hs.eg.db")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(clusterProfiler)

peak <- readPeakFile("Peaksets_minOverlap5.bed")
peak

# Promoter TSS enrichment heatmap
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red") # looks normal, like the tutorial

# Peak annotation
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
save(peakAnno, file = "../2_ChIPseeker_annotation/peakAnno.rda")
pdf("../2_ChIPseeker_annotation/peakAnno.pdf", height = 4, width = 6.5)
plotAnnoPie(peakAnno) # 17.43% Promoter +- 3kb (14% TSS +- 1kb)
dev.off()

# ----  used peakAnno.pdf rather than the one generated using the following commands -----
# Since UCSC-style chromosome names are used we have to change the style of the chromosome names from Ensembl to UCSC:
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

peakAnno.edb <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                             TxDb=edb, annoDb="org.Hs.eg.db")

pdf("../2_ChIPseeker_annotation/peakAnnoUCSC.pdf", height = 4, width = 6.5)
plotAnnoPie(peakAnno.edb) # 30.54% Promoter +- 3kb (21.84% TSS +- 1kb)
dev.off()
# ----------------