
##### pipeline #####
### (1) Determine the scATAC clusters that have enough reads for differential analysis - 20-30M reads per sample per cell cluster.
### (2) Load the differential scATAC peaks for each selected cluster
### (3) Assess the overlap between bulk DiffATAC and differential scATAC peaks

rm(list = ls())
DIR = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/scATACseq/results/04_ArchR_minTSS4_minFrag1000/10_Pseudobulk_SubclusterIdent_rm3BadSamples_NR4_minTSS5/BulkDiffATAC_overlap_diffSCatac/"
dir.create(DIR); setwd(DIR)

library(Repitools)
library(GenomicRanges)
library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(liftOver)

# --- On hoffman2: ----
##### (1) Determine the scATAC clusters that have enough reads for differential analysis 
# Yuyan: 20-30M reads per pseudobulk replicate per cell cluster.

library(ArchR)
archr_project = "/u/project/geschwind/jennybea/ASD_scATAC/04_ArchR/PeakCalling_SubclusterIdent_rm3BadSamples_NR4_minTSS5_C9"
project <- loadArchRProject(archr_project)

# project@ # tab
tmp = project@peakSet; str(tmp)
unique(tmp$GroupReplicate) # "ASD._.M9H3" "ASD._.2" "ASD._.1" "CTL._.Rep1" "CTL._.Rep2"  "ASD._.285A"
# but there is no read depth infor
tmp2 = project@cellColData; str(tmp2) # nFrags, Clusters for each cell, ReadsInPeaks
tmp3 = project@cellMetadata; str(tmp3) # nothing
# Hard to find out which cell was used for which pseudobulk, and whether the read depth is enough.

# --- give up after the above explorations ----

##### (2) Load the differential scATAC peaks for each selected cluster

### Load cluster-celltype match table
rm(list = ls())
lnames = load("../BulkDiffATAC_DistributionIn_scATACclusters/10_st3_v4_02_scATACmergedPeaks_membershipOfclusters.rda") # hits_scATAC_collapse, remapClust, remapClust2, remapClust3

# come back to correct the name of na3-ODClike cluster:
remapClust["C24"] = "EXT-ODC"

### Load differential scATAC peaks

fs = list.files(path = "../DiffPeaksUsingPseudobulkReplicates_SubclusterIdent_rm3BadSamples_NR4_minTSS5_allClusters", pattern = "marker_feature_list.*")
fs[1] # marker_feature_list_PeakCalling_SubclusterIdent_rm3BadSamples_NR4_minTSS5_C1_ASD.rds
cls = gsub("marker_feature_list_PeakCalling_SubclusterIdent_rm3BadSamples_NR4_minTSS5_", "", fs)
cls = gsub("_ASD\\.rds", "", cls)

i = 1; cl = cls[i]
diffSC = readRDS(paste0("../DiffPeaksUsingPseudobulkReplicates_SubclusterIdent_rm3BadSamples_NR4_minTSS5_allClusters/", fs[i]))
diffSC = as.data.frame(diffSC)
range(diffSC$FDR) # FDR < 0.1
colnames(diffSC)[1] = "chr"
diffSC = diffSC[,c(1:3, 6:7)]
diffSC$cl = cl

diffSC_all = diffSC

for (i in 2:length(cls)) {
  cl = cls[i]
  diffSC = readRDS(paste0("../DiffPeaksUsingPseudobulkReplicates_SubclusterIdent_rm3BadSamples_NR4_minTSS5_allClusters/", fs[i]))
  diffSC = as.data.frame(diffSC)
  
  if(nrow(diffSC) > 0) {
    colnames(diffSC)[1] = "chr"
    diffSC = diffSC[,c(1:3, 6:7)]
    diffSC$cl = cl
    
    diffSC_all = rbind(diffSC_all, diffSC)
  } else {
    ct = remapClust[cl]
    print(paste0(cl, "-", ct, " has 0 diff. scATAC peaks"))
  }
}
# print outs:
# [1] "C11-INT_ADARB2 has 0 diff. scATAC peaks"
# [1] "C8-BBB_Endo has 0 diff. scATAC peaks"

table(diffSC_all$cl)
#  C1   C10   C14   C19   C20   C22   C24   C25    C9 
# 116    29     3   629  8365   150 19904  5636  1258 

### Remember that scATAC is based on hg38, to overlap with bulk ATAC, I need to liftover to hg19
chain <- import.chain("~/Documents/Documents/Geschwind_lab/LAB/Computational_tools/HiC_related/R/liftOver/chain/hg38ToHg19.over.chain") # hg38ToHg19.over.chain downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/

diffSC_all_gr = annoDF2GR(diffSC_all)
diffSC_all_hg19_gr = unlist(liftOver(diffSC_all_gr, chain))
diffSC_all_hg19_df = as.data.frame(diffSC_all_hg19_gr)
colnames(diffSC_all_hg19_df)[1] = "chr"

save(diffSC_all, diffSC_all_hg19_df, remapClust, remapClust2, remapClust3, file = "11_02_diff_scATACpeaks_by_cluster_and_remapCluster123.rda")

##### (3) Assess the overlap between bulk DiffATAC and differential scATAC peaks

rm(list = ls())
lnames = load("11_02_diff_scATACpeaks_by_cluster_and_remapCluster123.rda") # "diffSC_all" "diffSC_all_hg19_df" "remapClust"  "remapClust2" "remapClust3"

diffSC_all = diffSC_all_hg19_df[,-c(4,5)]; rm(diffSC_all_hg19_df)
diffSC_all_gr = annoDF2GR(diffSC_all)
colnames(diffSC_all)[4:5] = c("scATAC_logFC", "scATAC_FDR")

### Load bulk DiffATAC peaks
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/DiffATAC.rda") # DiffATAC

fdr_t = 0.05; logfc_t = 0
ATACup = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC > 0), 1:4]
ATACdown = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC < 0), 1:4]

ATACup$bulk_category = "Up"; ATACdown$bulk_category = "Down"
ATACdiff = rbind(ATACup, ATACdown)
ATACdiff_gr = annoDF2GR(ATACdiff)
ATACdiff$DiffATAC = rownames(ATACdiff)

### findOverlaps
hits = as.data.frame(findOverlaps(ATACdiff_gr, diffSC_all_gr))
hits_df_DiffATAC = ATACdiff[hits$queryHits,]
hits_df_scATAC = diffSC_all[hits$subjectHits,]
hits_df = unique(cbind(hits_df_DiffATAC[,c("DiffATAC", "bulk_category", "ATAC_logFC")],
                hits_df_scATAC[,c("scATAC_logFC", "scATAC_FDR", "cl")]))

### Add major cell type memberships
hits_df$Clst = remapClust[match(hits_df$cl, names(remapClust))]
hits_df$MajorCt = remapClust2[match(hits_df$cl, names(remapClust2))]
hits_df$BroadCt = remapClust3[match(hits_df$cl, names(remapClust3))]

unique(hits_df$Clst) # "ODC-like (now EXT-ODC)" "ASTRO_3" "ASTRO_1" "ODC" "MG" "ASTRO-like" "OPC" "EXT"
#hits_df$Clst = factor(hits_df$Clst, levels = c("ASTRO_1", "ASTRO_3", "MG", "EXT", "ODC-like", "ODC"))
hits_df$Clst = factor(hits_df$Clst, levels = c("ASTRO_1", "ASTRO_3", "ASTRO-like", "MG", "EXT", "EXT-ODC", "ODC", "OPC"))

### Add colors
lnames = load("../BulkDiffATAC_DistributionIn_scATACclusters/10_st3_v4_02_count_scATACmergedPeaks_membershipOfMajorCts.rda") # count_scATAC, remapCol
rm(count_scATAC)
#hits_df$Ct_col = remapCol[match(hits_df$cl_mem, names(remapCol))]

### Plot logFC bulk vs scATAC, color by cell cluster
lm_res = lm(scATAC_logFC ~ 0 + ATAC_logFC, data = hits_df)
lm_summary = summary(lm_res)
lm_coef = round(lm_summary$coefficients[1,1],1) # 1.7
lm_p = formatC(lm_summary$coefficients[1,4], digits = 0, format = "e") # 5e-64
lm_adjr2 = round(lm_summary$adj.r.squared,2) # 0.16

lm_res2 = lm(scATAC_logFC ~ 0 + ATAC_logFC, data = hits_df[hits_df$scATAC_FDR < 0.05,])
lm_summary2 = summary(lm_res2)
lm_coef2 = round(lm_summary2$coefficients[1,1],1) # 2
lm_p2 = formatC(lm_summary2$coefficients[1,4], digits = 0, format = "e") # 4e-50
lm_adjr2_2 = round(lm_summary2$adj.r.squared,2) # 0.21
# Better use scATAC FDR < 0.05

lm_res3 = lm(scATAC_logFC ~ 0 + ATAC_logFC, data = hits_df[hits_df$scATAC_FDR < 0.1 & hits_df$Clst != "EXT-ODC",])
lm_summary3 = summary(lm_res3)
lm_coef3 = round(lm_summary3$coefficients[1,1],1) # 0.7
lm_p3 = formatC(lm_summary3$coefficients[1,4], digits = 0, format = "e") # 6e-2
lm_adjr2_3 = round(lm_summary3$adj.r.squared,2) # 0.01
# Better use scATAC FDR < 0.05

pdf("11_03_ggplot_logFC_bulk_vs_scATACpeaks_ColorBy_cluster.pdf", width = 6, height = 4)

hits_df[hits_df$scATAC_FDR < 0.05,] %>% # not as good, only ASTRO_3, EXT, ODC-like, ODC are left
  ggplot(aes(x = ATAC_logFC, y = scATAC_logFC)) +
  geom_point(aes(col = Clst), alpha = 0.8) +
  geom_abline(slope = lm_coef2, intercept = 0) +
  theme_bw() +
  scale_color_manual(values = c("yellow", "gold", "lightgoldenrod1", "darkorchid1", "brown1", "blue", "green", "darkolivegreen1")) +
  annotate("text", x = 0, y = -5, label = paste0("adjusted R2 = ", lm_adjr2_2, "\n", "p = ", lm_p2)) +
  xlab("Log2FC in bulk tissue ATAC-seq") +
  ylab("Log2FC in single-cell ATAC") +
  labs(col = "Cell cluster") +
  ggtitle("scATAC FDR < 0.05")

hits_df %>%
# hits_df[hits_df$scATAC_FDR < 0.05,] %>% # not as good, only ASTRO_3, EXT, ODC-like, ODC are left
  ggplot(aes(x = ATAC_logFC, y = scATAC_logFC)) +
  geom_point(aes(col = Clst), alpha = 0.8) +
  geom_abline(slope = lm_coef, intercept = 0) +
  theme_bw() +
  scale_color_manual(values = c("yellow", "gold", "lightgoldenrod1", "darkorchid1", "brown1", "blue", "green", "darkolivegreen1")) +
  annotate("text", x = 0, y = -5, label = paste0("adjusted R2 = ", lm_adjr2, "\n", "p = ", lm_p)) +
  xlab("Log2FC in bulk tissue ATAC-seq") +
  ylab("Log2FC in single-cell ATAC") +
  labs(col = "Cell cluster") +
  ggtitle("scATAC FDR < 0.1")

hits_df[hits_df$scATAC_FDR < 0.1 & hits_df$Clst != "EXT-ODC",] %>% # not as good, only ASTRO_3, EXT, ODC-like, ODC are left
  ggplot(aes(x = ATAC_logFC, y = scATAC_logFC)) +
  geom_point(aes(col = Clst), alpha = 0.8) +
  geom_abline(slope = lm_coef3, intercept = 0) +
  theme_bw() +
  scale_color_manual(values = c("yellow", "gold", "lightgoldenrod1", "darkorchid1", "brown1", "green", "darkolivegreen1")) +
  annotate("text", x = 0, y = -5, label = paste0("adjusted R2 = ", lm_adjr2_3, "\n", "p = ", lm_p3)) +
  xlab("Log2FC in bulk tissue ATAC-seq") +
  ylab("Log2FC in single-cell ATAC") +
  labs(col = "Cell cluster") +
  ggtitle("scATAC FDR < 0.1 (C24-na3 cluster excluded)")
# looks very bad without C24. The positive correlation are mostly driven by the up-reg ATAC peaks in EXT.

dev.off()
# Observation: Many ASD down-regulated bulk and single cell ATAC peaks are "ODC-like" (now named EXT-ODC), and that's driving the significantly positive correlation! Use FDR < 0.05

# No need to make bar plot given only 1668 bulk-sc differential ATAC pairs.

save(hits_df, lm_res, lm_res2, file = "11_03_logFC_bulk_vs_scATACpeaks_by_cluster_and_lmRes.rda")



