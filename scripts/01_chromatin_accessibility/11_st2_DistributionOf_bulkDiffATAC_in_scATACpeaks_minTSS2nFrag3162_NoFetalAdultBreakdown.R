
### Goal: Break down cell-type specificity of my bulk DiffATAC using diff. scATAC-seq peaks. (use FDR < 0.05 as this gives better cor in 11_scATAC_differential_analysis_usingPseudobulkReplicates_minTSS2nFrag3162.R)
### Note that Brie's scATAC-seq was done on GRCh38, while my bulk ATAC was on hg19.

##### pipeline #####
### (1) Load the GRanges objects of scATAC diff peaks
### (2) Table the counts of scATAC diff peaks
### (3) Load and build GRanges object for bulk diffATAC peaks (up-/down-reg) 
### (4) Make a dataframe of piechart distribution of bulk diffATAC peaks in scATAC diff peaks (based on major celltypes)
### (5) Plot the distribution as stacked bars

rm(list = ls())
#DIR = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/scATACseq/results/03_ArchR_minTSS2_minFrag3162/10_Pseudobulk/BulkDiffATAC_overlap_scDiffATAC/"
DIR = "/Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Projects/Project_ASD/scATACseq/results/03_ArchR_minTSS2_minFrag3162/10_Pseudobulk/BulkDiffATAC_overlap_scDiffATAC/"
dir.create(DIR, recursive = T); setwd(DIR)

library(Repitools)
library(GenomicRanges)
library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
#install.packages("bedr")
library(bedr)
check.binary("bedtools") # T
#BiocManager::install("liftOver")
library(liftOver)

##### (1) Load the GRanges objects of scATAC diff peaks 

lnames = load("11_02_diff_scATACpeaks_by_cluster_and_remapCluster.rda") # diffSC_all, diffSC_all_hg19_df, remapClust, remapClust2, remapClust3, remapClust4
diffSC_all_hg19_gr = annoDF2GR(diffSC_all_hg19_df[,-4])

##### (2) Table the counts of scATAC diff peaks

### Get merged scATAC diff peaks 
# Note: just for counting diff peaks and not for intersecting bulk DiffATAC peaks

diffSC_all_hg19_vector = paste0(diffSC_all_hg19_df$chr, ":", diffSC_all_hg19_df$start, "-", diffSC_all_hg19_df$end)

a.sort <- bedr.sort.region(diffSC_all_hg19_vector)
a.merge <- bedr.merge.region(a.sort) 
#  * Collapsing 82596 --> 71130 regions... NOTE, cool!
merged_peaks = as.data.frame(str_split_fixed(a.merge, pattern = ":|-", 3))
colnames(merged_peaks) = c("chr", "start", "end")
merged_peaks$start = as.integer(merged_peaks$start);
merged_peaks$end = as.integer(merged_peaks$end)
scATAC_merged_diffpeaks_gr = annoDF2GR(merged_peaks)

### Table the celltype distribution of scATAC diff peaks
hits_scATAC = as.data.frame(findOverlaps(scATAC_merged_diffpeaks_gr, diffSC_all_hg19_gr))
hits_scATAC$cl = diffSC_all_hg19_df$cl[hits_scATAC$subjectHits]

### Expand cell type information
hits_scATAC$Clst = remapClust[match(hits_scATAC$cl, names(remapClust))]
any(is.na(hits_scATAC$Clst)) # F, great.
hits_scATAC$MajorCt = remapClust2[match(hits_scATAC$cl, names(remapClust2))]
any(is.na(hits_scATAC$MajorCt)) # F, great.
hits_scATAC$stgBroadCt = remapClust3[match(hits_scATAC$cl, names(remapClust3))]
any(is.na(hits_scATAC$stgBroadCt)) # F, great.
hits_scATAC$BroadCt = remapClust4[match(hits_scATAC$cl, names(remapClust4))]
any(is.na(hits_scATAC$BroadCt)) # F, great.

### Make Venn Diagram of the clusters - ignore

colnames(hits_scATAC)
# [1] "queryHits" "subjectHits" "cl" "Clst" "MajorCt" "stgBroadCt" "BroadCt"
hits_scATAC_collapse = hits_scATAC %>%
  group_by(queryHits) %>%
  dplyr::select(queryHits, cl, Clst, MajorCt, stgBroadCt, BroadCt) %>%
  dplyr::mutate(
    Clsts = paste0(unique(Clst)[order(unique(Clst))], collapse = ","),
    MajorCts = paste0(unique(MajorCt)[order(unique(MajorCt))], collapse = ","),
    stgBroadCts = paste0(unique(stgBroadCt)[order(unique(stgBroadCt))], collapse = ","),
    BroadCts = paste0(unique(BroadCt)[order(unique(BroadCt))], collapse = ",")
  ) %>%
  ungroup()

length(unique(hits_scATAC_collapse$Clsts)) # 234 combinations of Cluster memberships (individual scATAC clusters)
length(unique(hits_scATAC_collapse$MajorCts)) # 64 combinations of Major celltype/cluster memberships (EXT/INT/MG/ASTRO/OPC/ODC/ENDO from unique(remapClust2))
length(unique(hits_scATAC_collapse$stgBroadCts)) # 15 combinations of Broad celltype memberships (Neuron/ASTRO/OPC-ODC/Endo). Turns out only BroadCts needs order.
length(unique(hits_scATAC_collapse$BroadCts)) # 3 combinations of Broad celltype memberships (Neuron/Glia/Endo). Turns out only BroadCts needs order.

len_Clsts = sapply(hits_scATAC_collapse$Clsts, function(x) length(unlist(str_split(x, ","))))
len_MajorCts = sapply(hits_scATAC_collapse$MajorCts, function(x) length(unlist(str_split(x, ","))))
len_stgBroadCts = sapply(hits_scATAC_collapse$stgBroadCts, function(x) length(unlist(str_split(x, ","))))
len_BroadCts = sapply(hits_scATAC_collapse$BroadCts, function(x) length(unlist(str_split(x, ","))))
hits_scATAC_collapse$Membership = ifelse(len_MajorCts == 1, hits_scATAC_collapse$MajorCts, ifelse(len_stgBroadCts == 1, hits_scATAC_collapse$stgBroadCts, hits_scATAC_collapse$BroadCts))

save(scATAC_merged_diffpeaks_gr, merged_peaks, hits_scATAC, hits_scATAC_collapse, remapClust, remapClust2, remapClust3, remapClust4, file = "11_st2_02_scATAC_merged_peaks_hg19_clusterMembership.rda")

### Count the celltype membership of the merged scATAC peaks 
hits_scATAC_collapse2 = unique(hits_scATAC_collapse[,c("queryHits", "Membership")])

count_scATAC = as.data.frame(table(hits_scATAC_collapse2$Membership))
colnames(count_scATAC) = c("cl_mem", "nATACpeaks")
count_scATAC = count_scATAC %>%
  arrange(desc(nATACpeaks))

# Add a row with cl_mem = "None"
count_scATAC[1,]
count_scATAC_0 = data.frame(cl_mem = "None", nATACpeaks = 0)
count_scATAC = rbind(count_scATAC, count_scATAC_0)

### Add colors 
# copied from 10_st3_DistributionOf_bulkDiffATAC_in_scATACclusters_minTSS2nFrag3162_NoFetalAdultBreakdown.R
# Color based on UMAP 03_ArchR_minTSS2_minFrag3162/05_harmony_cluster/08_cluster_umap.png

unique(count_scATAC$cl_mem)
# [1] EXT         ASTRO       OPC         Glia,Neuron MG          ODC        
# [7] Glia        INT         Neuron      OPC-ODC     INT_EXT     None       

remapCol <- c(
  #"ASTRO1" = "purple","ASTRO2" = "salmon","ENDO" = "",
  "INT_EXT" = "orchid",
  "INT_LHX6" = "yellow",
  "MG" = "aquamarine4",
  #"ODC1" = "lightblue","ODC2" = "lightgreen","ODC3" = "yellowgreen",
  "OPC" = "pink",
  #darkgoldenrod3,cornflowerblue,firebrick2

  "ASTRO" = "purple", # there is no unique peak in ASTRO1 or ASTRO2 alone
  "EXT" = "red",
  "INT" = "darkblue", # because there is no unique peak in INT_ADARB2 or INT_LHX6 alone
  "ODC" = "lightblue", # because there is no unique peak in ODC1 or ODC2 alone
  "OPC-ODC" = "deepskyblue",
  
  "Neuron" = "gray95",
  "Glia" = "gray75",
  #"Endo,Glia" = "gray75",Endo,Neuron" = "gray65",  "Endo,Glia,Neuron" = "gray45",
  "Glia,Neuron" = "gray55",
  
  "None" = "black"  
)
count_scATAC$Ct_col = remapCol[match(count_scATAC$cl_mem, names(remapCol))]
any(is.na(count_scATAC$Ct_col)) # F, great
#unique(count_scATAC$cl_mem[is.na(count_scATAC$Ct_col)])

count_scATAC$cl_mem = factor(count_scATAC$cl_mem, levels = rev(unique(count_scATAC$cl_mem)))

save(hits_scATAC_collapse2, count_scATAC, remapClust, remapClust2, remapClust3, remapClust4, remapCol, file = "11_st2_02_count_scATACmergedPeaks_CelltypeMembership.rda")

# Quick check of the distribution by ggplot 
plot_scATAC = count_scATAC %>%
  ggplot(aes(x = "", y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  ylab("Number of scATAC peaks") +
  xlab("Differential scATAC peaks") +
  theme_bw() +
  ggtitle("scATAC") +
  scale_fill_manual(values = rev(count_scATAC$Ct_col)) +
  labs(fill = "Cell-type membership") +
  guides(fill = guide_legend(ncol = 2)) # break the legends into two columns

pdf("11_st2_02_ggplot_clusterMembership_DiffscATAC_merged_peaks.pdf", width = 5, height = 5)
plot_scATAC
dev.off()

# Observations:
sum(count_scATAC$nATACpeaks[count_scATAC$cl_mem %in% c("Endo,Glia,Neuron", "Glia", "Neuron", "Endo,Glia", "Endo,Neuron")])/sum(count_scATAC$nATACpeaks) # 2% 
sum(count_scATAC$nATACpeaks[count_scATAC$cl_mem %in% c("Endo,Glia,Neuron", "Glia", "Neuron", "Endo,Glia", "Endo,Neuron", "OPC-ODC")])/sum(count_scATAC$nATACpeaks) # 2.4% 
# Almost all diff scATAC peaks are unique to the major cell types, only 2-2.5% overlap in multiple cell clusters.

##### (3) Load and build GRanges object for bulk diffATAC peaks (up-/down-reg)

rm(list = ls())
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/DiffATAC.rda") # DiffATAC

fdr_t = 0.05; logfc_t = 0
ATACup = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC > 0), 1:3]
ATACdown = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC < 0), 1:3]

ATACup$bulk_category = "Up"; ATACdown$bulk_category = "Down"
ATACdiff = rbind(ATACup, ATACdown)
ATACdiff_gr = annoDF2GR(ATACdiff)

##### (4) Make a piechart dataframe of distribution of each category of bulk diffATAC peaks in scATAC clusters (based on major celltypes)

### Get scATAC cluster membership for each bulk DiffATAC
lnames = load("11_02_diff_scATACpeaks_by_cluster_and_remapCluster.rda") # diffSC_all, diffSC_all_hg19_df, remapClust, remapClust2, remapClust3, remapClust4
diffSC_all_hg19_gr = annoDF2GR(diffSC_all_hg19_df[,-4])

hits = as.data.frame(findOverlaps(ATACdiff_gr, diffSC_all_hg19_gr))
ATACdiff$DiffATAC = rownames(ATACdiff)
hits_df_DiffATAC = ATACdiff[hits$queryHits,]
hits_df_scATAC = diffSC_all_hg19_df[hits$subjectHits,]
hits_df = cbind(hits_df_DiffATAC[,c("DiffATAC", "bulk_category")],
                hits_df_scATAC[,c("cl")]); colnames(hits_df)[ncol(hits_df)] = "cl"

# Add DiffATAC that does not overlap with any scATAC peak
any(! ATACdiff$DiffATAC %in% hits_df$DiffATAC) # T
hits_df_0 = ATACdiff[! ATACdiff$DiffATAC %in% hits_df$DiffATAC, c("DiffATAC", "bulk_category")]
hits_df_0$cl = "None"

hits_df = rbind(hits_df, hits_df_0)

### Extend cell type information
# Load celltype and col match table
lnames = load("11_st2_02_count_scATACmergedPeaks_CelltypeMembership.rda")
# hits_scATAC_collapse2, count_scATAC, remapClust, remapClust2, remapClust3, remapClust4, remapCol

hits_df$Clst = remapClust[match(hits_df$cl, names(remapClust))]
hits_df$MajorCt = remapClust2[match(hits_df$cl, names(remapClust2))]
hits_df$stgBroadCt = remapClust3[match(hits_df$cl, names(remapClust3))]
hits_df$BroadCt = remapClust4[match(hits_df$cl, names(remapClust4))]

colnames(hits_df)
# [1] "DiffATAC" "bulk_category" "cl" "Clst" "MajorCt" "stgBroadCt" "BroadCt"    
hits_collapse = hits_df %>%
  group_by(DiffATAC) %>%
  dplyr::mutate(
    Clsts = paste0(unique(Clst)[order(unique(Clst))], collapse = ","),
    MajorCts = paste0(unique(MajorCt)[order(unique(MajorCt))], collapse = ","),
    stgBroadCts = paste0(unique(stgBroadCt)[order(unique(stgBroadCt))], collapse = ","),
    BroadCts = paste0(unique(unlist(str_split(BroadCt, ",")))[order(unique(unlist(str_split(BroadCt, ","))))], collapse = ",")
  ) %>%
  ungroup()

len_Clsts = sapply(hits_collapse$Clsts, function(x) length(unlist(str_split(x, ","))))
len_MajorCts = sapply(hits_collapse$MajorCts, function(x) length(unlist(str_split(x, ","))))
len_stgBroadCts = sapply(hits_collapse$stgBroadCts, function(x) length(unlist(str_split(x, ","))))
len_BroadCts = sapply(hits_collapse$BroadCts, function(x) length(unlist(str_split(x, ","))))
hits_collapse$Membership = ifelse(len_MajorCts == 1, hits_collapse$MajorCts, ifelse(len_stgBroadCts == 1, hits_collapse$stgBroadCts, hits_collapse$BroadCts))

hits_collapse2 = unique(hits_collapse[,c("DiffATAC", "bulk_category", "Membership")]) # 5033 rows, including all the DiffATAC peaks, great.

save(hits_collapse, hits_collapse2, remapClust, remapClust2, remapClust3, remapClust4, remapCol, file = "11_st2_04_bulkDiffATAC_DiffscATACMembership.rda")

### Table the counts of bulk DiffATAC peaks in scATAC clusters
# Get the bulk DiffATAC peaks that do not overlap those of any scATAC cluster
rm(list = ls())
lnames = load("11_st2_04_bulkDiffATAC_DiffscATACMembership.rda") # hits_collapse, hits_collapse2, remapClust, remapClust2, remapClust3, remapClust4, remapCol

count_DiffATAC = as.data.frame(table(hits_collapse2[,c("bulk_category","Membership")]))
colnames(count_DiffATAC)[2:3] = c("cl_mem", "nATACpeaks")
count_DiffATAC = count_DiffATAC %>%
  arrange(bulk_category, desc(nATACpeaks))

### Add colors
count_DiffATAC$Ct_col = remapCol[match(count_DiffATAC$cl_mem, names(remapCol))]
any(is.na(count_DiffATAC$Ct_col)) # F, great
lnames = load("11_st2_02_count_scATACmergedPeaks_CelltypeMembership.rda") # hits_scATAC_collapse2, count_scATAC, remapClust, remapClust2, remapClust3, remapClust4, remapCol
count_DiffATAC$cl_mem = factor(count_DiffATAC$cl_mem, levels = levels(count_scATAC$cl_mem))

saveRDS(count_DiffATAC, file = "11_st2_04_count_bulkDiffATAC_DiffscATACMembership.rds")

# Quick check of the distribution by ggplot 
plot_DiffATAC = count_DiffATAC %>%
  ggplot(aes(x = bulk_category, y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("ASD vs. CTL") +
  ylab("Number of ATAC peaks") +
  theme_bw() +
  scale_fill_manual(values = rev(count_scATAC$Ct_col), drop = F) +
  ggtitle("Differential ATAC peaks in bulk tissue") +
  labs(fill = "Cell-type membership of\ndifferential scATAC peaks")
plot_DiffATAC

pdf("11_st2_04_ggplot_DiffscATACMembership_of_DiffATAC.pdf", width = 5, height = 5) # width=8 is too narrow
plot_DiffATAC
dev.off()

# Observations:
# Most (69%) bulk diffATAC peaks are not found as diff. scATAC peaks.
# Among the ones that were mapped to diff scATAC peaks, only 20-25% were shared across different major cell types, others were uniquely differential in specific cell types. Can't really control the power of diff scATAC detection either.
# As Dan pointed out, celltypes that are more abundant in the brain will have more power on bulk ATAC differential analysis. So we can't say much about the relative proportion of diff scATAC in up/down-reg Bulk DiffATAC, except that we see they were represented by all the different cell types. OPC/ODC ATAC peaks tend to be down-reg rather than up.
(a = sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("None")])/nrow(hits_collapse2)) # 69%
(b = sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("Neuron","Glia","Glia,Neuron")])/nrow(hits_collapse2)) # 7.3%
b/(1-a) # 23.3%

sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("ODC", "OPC-ODC","OPC") & count_DiffATAC$bulk_category == "Up"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Up"]) # 2.4%
sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("ODC", "OPC-ODC","OPC") & count_DiffATAC$bulk_category == "Down"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Down"]) # 9.5%

sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("EXT", "Neuron") & count_DiffATAC$bulk_category == "Up"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Up"]) # 10.6%
sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("EXT", "Neuron") & count_DiffATAC$bulk_category == "Down"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Down"]) # 10.5%

sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("ASTRO") & count_DiffATAC$bulk_category == "Up"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Up"]) # 6.4%
sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("ASTRO") & count_DiffATAC$bulk_category == "Down"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Down"]) # 8%

sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("MG") & count_DiffATAC$bulk_category == "Up"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Up"]) # 1.2%
sum(count_DiffATAC$nATACpeaks[count_DiffATAC$cl_mem %in% c("MG") & count_DiffATAC$bulk_category == "Down"])/sum(count_DiffATAC$nATACpeaks[count_DiffATAC$bulk_category == "Down"]) # 0.75%

##### (5) Plot the distribution as stacked bars
rm(list = ls())
count_DiffATAC = readRDS("11_st2_04_count_bulkDiffATAC_DiffscATACMembership.rds")
lnames = load("11_st2_02_count_scATACmergedPeaks_CelltypeMembership.rda") # count_scATAC, hits_scATAC_collapse2, remapClust, remapClust2, remapClust3, remapClust4, remapCol

plot_scATAC = count_scATAC %>%
  ggplot(aes(x = "", y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  ylab("Number of scATAC peaks") +
  xlab("Differential scATAC peaks") +
  theme_bw() +
  ggtitle("scATAC") +
  scale_fill_manual(values = rev(count_scATAC$Ct_col)) +
  labs(fill = "Cell-type membership") #+
  #guides(fill = guide_legend(ncol = 2)) # break the legends into two columns
plot_scATAC

plot_DiffATAC = count_DiffATAC %>%
  ggplot(aes(x = bulk_category, y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("ASD vs. CTL") +
  ylab("Number of ATAC peaks") +
  theme_bw() +
  scale_fill_manual(values = rev(count_scATAC$Ct_col), drop = F) +
  ggtitle("Differential ATAC peaks in bulk tissue") +
  labs(fill = "Cell-type membership of\ndifferential scATAC peaks")
plot_DiffATAC

TextSize = 5
#lg = get_legend(plot_DiffATAC)

pdf("11_st2_05_ggplot_clusterMembership_DiffscATACandDiffATAC.pdf", width = 8, height = 5.5)
plot_grid(plot_scATAC + theme(legend.position = "none"), 
          plot_DiffATAC,
          nrow = 1, rel_widths = c(1,2),
          label_size = TextSize
          )
dev.off()

# Observation - from previous sections really

