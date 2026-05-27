
### Goal: Ask cell-type specificity of my bulk DiffATAC in Brie's scATAC-seq clusters.
### Note that Brie's scATAC-seq was done on GRCh38, while my bulk ATAC was on hg19.
## v2: C24 labeled as a mixture EXT-ODC

##### pipeline #####
### (1) Load GRanges objects and cell-type specificity counts of scATAC pseudobulk peaks and 
### (2) Load and build GRanges object for bulk diffATAC peaks (up-/down-reg) 
### (3) Separate bulk diffATAC peaks into categories: Fetal/Adult/Previously unknown 
### (4) Make a piechart dataframe of distribution of each category of bulk diffATAC peaks in scATAC clusters (based on major celltypes)
### (5) Add number of known adult high-confidence and fetal ATAC peaks
### (6) Plot the distribution as stacked bars
### (7) Can we link the up-reg EXT novel peaks to target genes based on promoters

rm(list = ls())
DIR = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/scATACseq/results/04_ArchR_minTSS4_minFrag1000/10_Pseudobulk_SubclusterIdent_rm3BadSamples_NR4_minTSS5/BulkDiffATAC_DistributionIn_scATACclusters/"
dir.create(DIR); setwd(DIR)

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
#library(liftOver)

##### (1) Load GRanges objects and cell-type specificity counts of scATAC pseudobulk peaks and
lnames = load("10_st3_01_scATAC_cluster_ASDandCTLmerged_peaks_hg19.rda")
# scATAC_cluster_peaks_hg19_gr, scATAC_cluster_peaks_hg19_df
scATAC_cluster_peaks_gr = scATAC_cluster_peaks_hg19_gr
scATAC_cluster_peaks_df = scATAC_cluster_peaks_hg19_df
rm(scATAC_cluster_peaks_hg19_gr, scATAC_cluster_peaks_hg19_df)

lnames = load("10_st3_02_scATAC_merged_peaks_hg19_clusterMembership.rda") # scATAC_merged_peaks_gr, merged_peaks, hits_scATAC_collapse
lnames = load("10_st3_v4_02_scATACmergedPeaks_membershipOfclusters.rda") # hits_scATAC_collapse, remapClust, remapClust2, remapClust3
rm(hits_scATAC_collapse)

remapClust["C24"] = "EXT-ODC"
remapClust2["C24"] = "EXT-ODC"
remapClust3["C24"] = "Glia,Neuron"

hits_scATAC = as.data.frame(findOverlaps(scATAC_merged_peaks_gr, scATAC_cluster_peaks_gr))
hits_scATAC$cl = scATAC_cluster_peaks_df$cl[hits_scATAC$subjectHits]

hits_scATAC$Clst = remapClust[match(hits_scATAC$cl, names(remapClust))]
hits_scATAC$MajorCt = remapClust2[match(hits_scATAC$cl, names(remapClust2))]
hits_scATAC$BroadCt = remapClust3[match(hits_scATAC$cl, names(remapClust3))]

colnames(hits_scATAC)
# [1] "queryHits" "subjectHits" "cl" "Clst" "MajorCt" "BroadCt"    
hits_scATAC_collapse = hits_scATAC %>%
  dplyr::select(queryHits, cl, Clst, MajorCt, BroadCt) %>%
  group_by(queryHits) %>%
  dplyr::mutate( # must use dplyr::mutate
    #Clsts = paste0(Clst, collapse = ",")#,
    Clsts = paste0(unique(Clst), collapse = ","),
    MajorCts = paste0(unique(MajorCt)[order(unique(MajorCt))], collapse = ","),
    BroadCts = paste0(unique(unlist(str_split(BroadCt, ",")))[order(unique(unlist(str_split(BroadCt, ","))))], collapse = ",")
  ) %>%
  ungroup()

unique(hits_scATAC_collapse$BroadCts) 
# [1] "Neuron"           "Glia"             "Glia,Neuron"      "Endo,Glia,Neuron"
# [5] "Endo"             "Endo,Neuron"      "Endo,Glia"       
# looks good

length(unique(hits_scATAC_collapse$Clsts)) # 1793 combinations of Cluster memberships (individual scATAC clusters)
length(unique(hits_scATAC_collapse$MajorCts)) # 255 combinations of Major celltype/cluster memberships (EXT/INT/MG/ASTRO/OPC/ODC/Endo/EXT-ODC from unique(remapClust2))
length(unique(hits_scATAC_collapse$BroadCts)) # 7 combinations of Broad celltype memberships (Neuron/Glia/Endo). Turns out only BroadCts needs order.

len_Clsts = sapply(hits_scATAC_collapse$Clsts, function(x) length(unlist(str_split(x, ","))))
len_MajorCts = sapply(hits_scATAC_collapse$MajorCts, function(x) length(unlist(str_split(x, ","))))
len_BroadCts = sapply(hits_scATAC_collapse$BroadCts, function(x) length(unlist(str_split(x, ","))))
hits_scATAC_collapse$Membership = ifelse(len_MajorCts == 1, hits_scATAC_collapse$MajorCts, hits_scATAC_collapse$BroadCts)

# quick check
tmp = hits_scATAC_collapse[hits_scATAC_collapse$Clsts == "EXT-ODC",] # membership is EXT-ODC in v2
unique(tmp$Membership) # "EXT-ODC", no overlap with other cell clusters ...

length(unique(hits_scATAC_collapse$Membership)) 
# old len_MajorCts <= 3: 73 combinations of cluster/celltype memberships, a bit long for legends. len_MajorCts <= 2: 39, still a bit long for legends.
# old v3 len_MajorCts == 1: 19 combinations of cluster/celltype memberships. this is fine
# now based on MajorCts: 14
(tmp = unique(hits_scATAC_collapse$Membership)) # many single cell clusters as well

save(hits_scATAC_collapse, remapClust, remapClust2, remapClust3, file = "10_st4_v2_01_scATACmergedPeaks_membershipOfclusters.rda")

### count
hits_scATAC_collapse2 = unique(hits_scATAC_collapse[,c("queryHits", "Membership")])
count_scATAC = as.data.frame(table(hits_scATAC_collapse2$Membership))
colnames(count_scATAC) = c("cl_mem", "nATACpeaks")
count_scATAC = count_scATAC %>%
  arrange(desc(nATACpeaks))

# Add a row with cl_mem = "None"
count_scATAC[1,]
count_scATAC_0 = data.frame(cl_mem = "None", nATACpeaks = 0)
count_scATAC = rbind(count_scATAC, count_scATAC_0)

unique(count_scATAC$cl_mem)
unique(count_scATAC$cl_mem)[order(unique(count_scATAC$cl_mem))]

remapCol <- c(
  "ASTRO" = "yellow",
  "Endo" = "cyan",
  "MG" = "darkorchid1",
  "ODC" = "green",
  #"ODC-like" = "darkseagreen1", # na3
  "OPC" = "darkolivegreen1",
  "EXT-ODC" = "blue",
  "EXT" = "brown1",
  "INT" = "pink",
  #"INT_LHX6" = "lightblue",
  #"INT_ADARB2" = "deepskyblue",
  
  "Glia" = "burlywood3",
  "Neuron" = "bisque",
  "Endo,Glia,Neuron" = "darkgoldenrod1",
  "Glia,Neuron" = "burlywood1",
  "Endo,Glia" = "goldenrod4",
  "Endo,Neuron" = "goldenrod3",
  
  "None" = "grey"
) # No Neuron, cuz that would be EXT,INT
count_scATAC$cl_mem = gsub(" ","", count_scATAC$cl_mem)
count_scATAC$Ct_col = remapCol[match(count_scATAC$cl_mem, names(remapCol))]
any(is.na(count_scATAC$Ct_col)) # F

count_scATAC$cl_mem = factor(count_scATAC$cl_mem, levels = rev(unique(count_scATAC$cl_mem)))

save(count_scATAC, remapCol, file = "10_st4_v2_01_count_scATACmergedPeaks_membershipOfMajorCts.rda")

##### (2) Load and build GRanges object for bulk diffATAC peaks (up-/down-reg)

rm(list = ls())
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/DiffATAC.rda") # DiffATAC

fdr_t = 0.05; logfc_t = 0
ATACup = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC > 0), 1:3]
ATACdown = DiffATAC[which(DiffATAC$ATAC_FDR < fdr_t & DiffATAC$ATAC_logFC < 0), 1:3]

ATACup$bulk_category = "Up"; ATACdown$bulk_category = "Down"
ATACdiff = rbind(ATACup, ATACdown)
ATACdiff_gr = annoDF2GR(ATACdiff)

##### (3) Separate bulk diffATAC peaks into categories: Fetal/Adult/Previously unknown 

### Load PsychEncode Wang et al adult ATAC peaks
adult_peaks = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/enhancer_database_copy/PsychENCODE/Wang_DER-03b_hg19_high_confidence_PEC_enhancers.bed")

### Load Luis fetal ATAC peaks
fetal_peaks = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/enhancer_database_copy/PsychENCODE/Luis_GSE95023_readswithinpeaks.bed")

### Overlap bulk DiffATAC peaks with fetal and adult ATAC peaks
colnames(adult_peaks)[1:3] = colnames(fetal_peaks)[1:3] = c("chr", "start", "end")
adult_peaks_gr = annoDF2GR(adult_peaks)
fetal_peaks_gr = annoDF2GR(fetal_peaks)

hits_DiffATAC_adult = as.data.frame(findOverlaps(ATACdiff_gr, adult_peaks_gr))
hits_DiffATAC_fetal = as.data.frame(findOverlaps(ATACdiff_gr, fetal_peaks_gr))

idx_const = unique(c(hits_DiffATAC_adult$queryHits, hits_DiffATAC_fetal$queryHits)) # 2799
idx_adult = setdiff(idx_const, unique(hits_DiffATAC_fetal$queryHits)) # 733
idx_fetal = setdiff(idx_const, unique(hits_DiffATAC_adult$queryHits)) # 1691
idx_novel = setdiff(1:nrow(ATACdiff), idx_const) # 2234

ATACdiff$Known = "NA"
ATACdiff$Known[idx_const] = "Constitutive"
ATACdiff$Known[idx_adult] = "Adult"
ATACdiff$Known[idx_fetal] = "Fetal"
ATACdiff$Known[idx_novel] = "Novel"
any(is.na(ATACdiff$Known)) # F
ATACdiff_gr = annoDF2GR(ATACdiff)

save(ATACdiff, ATACdiff_gr, file = "10_st4_v2_03_bulkDiffATAC_Category_ConstAdultFetalNovel.rda") # same as non v2

##### (4) Make a piechart dataframe of distribution of each category of bulk diffATAC peaks in scATAC clusters (based on major celltypes)

rm(list = ls())
lnames = load("10_st4_v2_03_bulkDiffATAC_Category_ConstAdultFetalNovel.rda") # ATACdiff, ATACdiff_gr

### Get scATAC cluster membership for each bulk DiffATAC
lnames = load("10_st3_01_scATAC_cluster_ASDandCTLmerged_peaks_hg19.rda") # scATAC_cluster_peaks_hg19_gr, scATAC_cluster_peaks_hg19_df
scATAC_cluster_peaks_gr = scATAC_cluster_peaks_hg19_gr
scATAC_cluster_peaks_df = scATAC_cluster_peaks_hg19_df
rm(scATAC_cluster_peaks_hg19_gr, scATAC_cluster_peaks_hg19_df)

hits = as.data.frame(findOverlaps(ATACdiff_gr, scATAC_cluster_peaks_gr))
ATACdiff$DiffATAC = rownames(ATACdiff)
hits_df_DiffATAC = ATACdiff[hits$queryHits,]
hits_df_scATAC = scATAC_cluster_peaks_df[hits$subjectHits,]
hits_df = cbind(hits_df_DiffATAC[,c("DiffATAC", "bulk_category", "Known")],
                hits_df_scATAC[,c("cl")]); colnames(hits_df)[ncol(hits_df)] = "cl"

### Load celltype and col match table
lnames = load("10_st4_v2_01_scATACmergedPeaks_membershipOfclusters.rda")
# hits_scATAC_collapse, remapClust, remapClust2, remapClust3

hits_df$Clst = remapClust[match(hits_df$cl, names(remapClust))]
hits_df$MajorCt = remapClust2[match(hits_df$cl, names(remapClust2))]
hits_df$BroadCt = remapClust3[match(hits_df$cl, names(remapClust3))]

colnames(hits_df)
# [1] "DiffATAC" "bulk_category" "Known" "cl" "Clst" "MajorCt" "BroadCt"    
hits_collapse = hits_df %>%
  group_by(DiffATAC) %>%
  dplyr::mutate(
    Clsts = paste0(unique(Clst), collapse = ","),
    MajorCts = paste0(unique(MajorCt)[order(unique(MajorCt))], collapse = ","),
    BroadCts = paste0(unique(unlist(str_split(BroadCt, ",")))[order(unique(unlist(str_split(BroadCt, ","))))], collapse = ",")
  ) %>%
  ungroup()

len_Clsts = sapply(hits_collapse$Clsts, function(x) length(unlist(str_split(x, ","))))
len_MajorCts = sapply(hits_collapse$MajorCts, function(x) length(unlist(str_split(x, ","))))
len_BroadCts = sapply(hits_collapse$BroadCts, function(x) length(unlist(str_split(x, ","))))
hits_collapse$Membership = ifelse(len_MajorCts == 1, hits_collapse$MajorCts, hits_collapse$BroadCts)

hits_collapse2 = unique(hits_collapse[,c("DiffATAC", "bulk_category", "Known", "Membership")])
hits_0 = ATACdiff[which(! ATACdiff$DiffATAC %in% hits_collapse2$DiffATAC), c("DiffATAC", "bulk_category", "Known")]
hits_0$Membership = "None"
hits_collapse2 = rbind(hits_collapse2, hits_0) # 5033 rows now, including all the DiffATAC peaks.

save(hits_collapse, hits_collapse2, file = "10_st4_v2_04_bulkDiffATAC_ConstAdultFetalNovel_membershipOfclusters.rda")

### Table the counts of bulk DiffATAC peaks in scATAC clusters
# Get the bulk DiffATAC peaks that do not overlap those of any scATAC cluster
rm(list = ls())
lnames = load("10_st4_v2_04_bulkDiffATAC_ConstAdultFetalNovel_membershipOfclusters.rda") # hits_collapse, hits_collapse2

count_DiffATAC = as.data.frame(table(hits_collapse2[,c("bulk_category","Known","Membership")]))
colnames(count_DiffATAC)[3:4] = c("cl_mem", "nATACpeaks")
count_DiffATAC = count_DiffATAC %>%
  arrange(bulk_category, rev(Known), desc(nATACpeaks))

### Add colors
lnames = load("10_st4_v2_01_count_scATACmergedPeaks_membershipOfMajorCts.rda") # count_scATAC, remapCol
count_DiffATAC$Ct_col = remapCol[match(count_DiffATAC$cl_mem, names(remapCol))]
#count_DiffATAC$cl_mem = factor(count_DiffATAC$cl_mem, levels = rev(unique(count_DiffATAC$cl_mem)))
any(is.na(count_DiffATAC$Ct_col)) # F
count_DiffATAC$cl_mem = factor(count_DiffATAC$cl_mem, levels = levels(count_scATAC$cl_mem))

saveRDS(count_DiffATAC, file = "10_st4_v2_04_count_bulkDiffATAC_ConstAdultFetalNove_membershipOfMajorCts.rds")

##### (5) Add number of known adult high-confidence and fetal ATAC peaks
# based on 5_9_st1_diffATAC_overlap_FetalvsAdultEnhancers.R

rm(list = ls())

### Load PsychEncode Wang et al adult ATAC peaks
adult_peaks = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/enhancer_database_copy/PsychENCODE/Wang_DER-03b_hg19_high_confidence_PEC_enhancers.bed")

### Load Luis fetal ATAC peaks
fetal_peaks = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/enhancer_database_copy/PsychENCODE/Luis_GSE95023_readswithinpeaks.bed")

### Overlap adult and fetal ATAC peaks
colnames(adult_peaks)[1:3] = colnames(fetal_peaks)[1:3] = c("chr", "start", "end")
adult_peaks_gr = annoDF2GR(adult_peaks)
fetal_peaks_gr = annoDF2GR(fetal_peaks)

hits = as.data.frame(findOverlaps(adult_peaks_gr, fetal_peaks_gr))

n_adult_const = length(unique(hits$queryHits))
n_fetal_const = length(unique(hits$subjectHits))
#n_adult_specific = length(setdiff(1:nrow(adult_peaks), unique(hits$queryHits))) # same as the following:
n_adult_specific = nrow(adult_peaks) - n_adult_const
n_fetal_specific = nrow(fetal_peaks) - n_fetal_const

df = data.frame(Tissue = rep(c("Adult", "Fetal"), each = 2), 
                Category = rep(c("Constitutive", "Specific"), 2),
                nATACpeaks = c(n_adult_const, n_adult_specific, n_fetal_const, n_fetal_specific))

saveRDS(df, file = "10_st4_v2_05_Known_AdultFetalConst_nATACpeaks.rds") # same as non v2

##### (6) Plot the distribution as stacked bars
rm(list = ls())
count_ATAC = readRDS("10_st4_v2_04_count_bulkDiffATAC_ConstAdultFetalNove_membershipOfMajorCts.rds")
lnames = load("10_st4_v2_01_count_scATACmergedPeaks_membershipOfMajorCts.rda") # count_scATAC, remapCol
count_KnownTissue = readRDS("10_st4_v2_05_Known_AdultFetalConst_nATACpeaks.rds")

plot_KnownTissue = count_KnownTissue %>%
  ggplot(aes(x = Tissue, y = nATACpeaks, fill = Category)) +
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  ylab("Number of ATAC peaks") +
  xlab("Brain tissue") +
  theme_bw() +
  ggtitle("Known ATAC peaks") +
  labs(fill = "Tissue\nspecificity")
plot_KnownTissue

plot_scATAC = count_scATAC %>%
  ggplot(aes(x = "", y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  ylab("Number of ATAC peaks") +
  xlab("Pseudobulk ATAC peaks") +
  theme_bw() +
  ggtitle("scATAC") +
  scale_fill_manual(values = rev(count_scATAC$Ct_col)) +
  labs(fill = "Cell-type membership") +
  guides(fill = guide_legend(ncol = 2)) # break the legends into two columns
plot_scATAC

plot_ATACdown = count_ATAC[count_ATAC$bulk_category == "Down",] %>%
  ggplot(aes(x = Known, y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity") +
  #labs(fill = "Enhancer tissue specificity") +
  xlab("Overlap with known ATAC peaks") +
  ylab("Number of ATAC peaks") +
  theme_bw() +
  #scale_fill_manual(values = rev(unique(count_ATAC$Ct_col))) +
  scale_fill_manual(values = rev(count_scATAC$Ct_col)) +
  ggtitle("Down-regulated in ASD bulk tissue") +
  labs(fill = "Cell-type membership")
plot_ATACdown

plot_ATACup = count_ATAC[count_ATAC$bulk_category == "Up",] %>%
  ggplot(aes(x = Known, y = nATACpeaks, fill = cl_mem)) +
  geom_bar(position = "stack", stat = "identity") +
  #labs(fill = "Enhancer tissue specificity") +
  xlab("Overlap with known ATAC peaks") +
  ylab("Number of ATAC peaks") +
  theme_bw() +
  #scale_fill_manual(values = rev(unique(count_ATAC$Ct_col))) +
  scale_fill_manual(values = rev(count_scATAC$Ct_col)) +
  ggtitle("Up-regulated in ASD bulk tissue") +
  labs(fill = "Cell-type membership")
plot_ATACup

TextSize = 5
lg = get_legend(plot_scATAC)

plot_row1 = plot_grid(plot_KnownTissue, 
                      plot_scATAC,
                      nrow = 1, 
                      rel_widths = c(1,1.5),
                      label_size = TextSize)

plot_row2 = plot_grid(plot_ATACdown + theme(legend.position = "none"), 
                      plot_ATACup + theme(legend.position = "none"),nrow = 1, 
                      label_size = TextSize)

pdf("10_st4_v2_06_DiffATAC_intersect_scATACpseudobulk.pdf", width = 7, height = 5.5)
plot_grid(plot_row1, plot_row2,
          nrow = 2#,rel_heights = c(1,2)
)
dev.off()

# Observation: Nearly all ASD down-reg ATAC peaks were found across multiple cell types, and now the ODC novel peaks faded away (because they are actually EXT-ODC cluster), together with the Glia peaks (probably become Glia,Neuron). 60% of ASD up-reg ATAC peaks were found in multiple cell types, but 30% of the novel ATAC peaks were specifically found in excitatory neurons.
sum(count_ATAC$nATACpeaks[count_ATAC$cl_mem %in% c("Glia,Neuron", "Endo,Glia,Neuron") & count_ATAC$bulk_category == "Up"])/sum(count_ATAC$nATACpeaks[count_ATAC$bulk_category == "Up"]) # 61%
sum(count_ATAC$nATACpeaks[count_ATAC$cl_mem %in% c("Neuron", "EXT", "INT") & count_ATAC$bulk_category == "Up"])/sum(count_ATAC$nATACpeaks[count_ATAC$bulk_category == "Up"]) # 33%

### Question: how do we integrate this result into my ASD manuscript?
# I think what this suggests is that Differential bulk ATAC analysis is finding signals shared across different cell types. Cell-type specific changes from scATAC can be separated into a different paper.
# The result suggests that many fetal ATAC peaks were not shut down properly during neural development. The novel EXT peaks that were up-reg could be compensatory or not - need to see whether they are promoter ATAC regions and can be linked to gene expression.

##### (7) Can we link the up-reg EXT novel peaks to target genes based on promoters
rm(list = ls())

lnames = load("10_st4_v2_04_bulkDiffATAC_ConstAdultFetalNovel_membershipOfclusters.rda") # hits_collapse, hits_collapse2

# based on Manuscripts/ASD_Dup15q_project/Rscripts_v3_diffIS/Figure1_SupFigure3_L45terms_ghCorrected.R
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda") # sig_CorATAC, sig_CorATAC_uniq, promATAC

ATACup_Novel_EXT = hits_collapse2$DiffATAC[which(hits_collapse2$bulk_category == "Up" & hits_collapse2$Known == "Novel" & hits_collapse2$Membership == "EXT")] # 454 ATAC peaks
idx = which(promATAC$ATACpeak %in% ATACup_Novel_EXT) # only 8
promATACup_Novel_EXT = promATAC[idx,] # 2/8 genes are significantly down-reg (FDR < 0.05, or 0.1 the same results), the other 6 show no gene expression change (FDR > 0.8).

# Conclusion: Hard to interpret the ATAC changes.

### Next step, perform diff analysis on scATAC-seq and check cell-type specificity of bulk DiffATAC.


