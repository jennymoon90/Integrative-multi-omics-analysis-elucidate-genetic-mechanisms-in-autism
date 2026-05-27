rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(VennDiagram)
library(UpSetR)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/BINDetect")

load("ATACupORdown_fdr01logfc02_overlap_TFoiBound.rda") # ATACup ATACdown
ATACdown$ATACpeak = paste0(ATACdown$chr, "_", ATACdown$start, "_", ATACdown$end)
ATACup$ATACpeak = paste0(ATACup$chr, "_", ATACup$start, "_", ATACup$end)
TFs_lowerBinding = c("CTCF", "ETV5", "ELK4", "SOX9", "NR2F1")

## 1) Venn Diagram of # Down-reg ATAC peaks co-bound by the 5 TFs of interest
CTCF_boundATACdown = ATACdown$ATACpeak[which(ATACdown[,grepl("CTCF", colnames(ATACdown))] == 1)] # 744
ETV5_boundATACdown = ATACdown$ATACpeak[which(ATACdown[,grepl("ETV5", colnames(ATACdown))] == 1)] # 92
ELK4_boundATACdown = ATACdown$ATACpeak[which(ATACdown[,grepl("ELK4", colnames(ATACdown))] == 1)] # 87
SOX9_boundATACdown = ATACdown$ATACpeak[which(ATACdown[,grepl("SOX9", colnames(ATACdown))] == 1)] # 281
NR2F1_boundATACdown = ATACdown$ATACpeak[which(ATACdown[,grepl("NR2F1", colnames(ATACdown))] == 1)] # 169

venn.diagram(x = list(CTCF_boundATACdown, ETV5_boundATACdown, ELK4_boundATACdown, SOX9_boundATACdown, NR2F1_boundATACdown), category.names = TFs_lowerBinding, filename = "VennDiagram_OverlapOfNumATACdownBroundBy5TFs.png", fill = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF", "gold"), cat.pos = c(0, -20, -120, 90, 20), cat.fontface = "bold", cex = 1.2, lwd = c(3, rep(1, 4))) # cat.pos isn't about left and right. It's vector giving the position (in degrees) of each category name along the circle, with 0 at 12 o'clock

# May be clear use upset
pdf("UpSetR_OverlapOfNumATACdownBroundBy5TFs.pdf", width = 7, height = 5, onefile = F)
upset(fromList(list(CTCF = CTCF_boundATACdown, ETV5 = ETV5_boundATACdown, ELK4 = ELK4_boundATACdown, SOX9 = SOX9_boundATACdown, NR2F1 = NR2F1_boundATACdown)), order.by = "freq")
dev.off()
# Observation:
# Uniquely bound by CTCF > SOX9 > NR2F1 > ELK4/ETV5 (same TF cluster and share more than unique)
# Many down-reg ATACs are co-bound by CTCF and other factors.

# VennDiagram and UpSetR automatically takes care of the duplicated entries.

#length(ATACdown$ATACpeak[which(ATACdown$`CTCF-MA0139.1_overview_ASCorCTLbound.txt` == 1 & ATACdown$`ELK4-MA0076.1_overview_ASCorCTLbound.txt` == 1)]) # 35 co-bound by CTCF and ELK4.  It's safer using unique(). 
length(unique(ATACdown$ATACpeak[which(ATACdown$`CTCF-MA0139.1_overview_ASCorCTLbound.txt` == 1 & ATACdown$`ELK4-MA0076.1_overview_ASCorCTLbound.txt` == 1)])) # 35 co-bound by CTCF and ELK4
length(unique(ATACdown$ATACpeak[which(ATACdown$`CTCF-MA0139.1_overview_ASCorCTLbound.txt` == 1 & ATACdown$`ETV5-MA0765.1_overview_ASCorCTLbound.txt` == 1)])) # 40 co-bound by CTCF and ETV5
length(unique(ATACdown$ATACpeak[which(ATACdown$`CTCF-MA0139.1_overview_ASCorCTLbound.txt` == 1 & ATACdown$`NR2F1-MA0017.2_overview_ASCorCTLbound.txt` == 1)])) # 82 co-bound by CTCF and NR2F1
length(unique(ATACdown$ATACpeak[which(ATACdown$`CTCF-MA0139.1_overview_ASCorCTLbound.txt` == 1 & ATACdown$`SOX9-MA0077.1_overview_ASCorCTLbound.txt` == 1)])) # 65 co-bound by CTCF and SOX9

## 2) Venn Diagram of # Up-reg ATAC peaks co-bound by the 5 TFs of interest
# TF names
TFs_higherBinding = 
  setdiff(colnames(ATACdown)[4:ncol(ATACdown)], 
          colnames(ATACdown)[grepl(paste(TFs_lowerBinding, collapse = "|"), colnames(ATACdown))])
TFs_higherBinding = sub("-MA.*", "", TFs_higherBinding)
TFs_higherBinding = sub("_HUMAN.*", "", TFs_higherBinding) # 34
TFs_higherBinding = sub("var.*", "", TFs_higherBinding)
TFs_higherBinding = unique(TFs_higherBinding) # 28

# List of TF bound up-reg ATAC peaks
TF_higherBinding_ATACpeaks = vector("list", length(TFs_higherBinding))
names(TF_higherBinding_ATACpeaks) = TFs_higherBinding

for (i in 1:length(TFs_higherBinding)) {
  TF = TFs_higherBinding[i]
  ColNames = which(grepl(paste(c(paste0(TF,"_HUMAN"), paste0(TF,"-MA"), paste0(TF,"var")), collapse = "|"), colnames(ATACup)))
  idx = which(apply(ATACup[,ColNames, drop = F], 1, function(x) any(x == 1)))
  TF_higherBinding_ATACpeaks[[i]] = unique(ATACup$ATACpeak[idx])
}

save(TF_higherBinding_ATACpeaks, file = "UpSetR_ATACupPeaks_boundbyTFhigherBinding.rda")

# UpSetR
pdf("UpSetR_OverlapOfNumATACupBroundByTFs.pdf", width = 14, height = 10, onefile = F)
upset(fromList(TF_higherBinding_ATACpeaks), order.by = "freq", nsets = length(TFs_higherBinding)) # , (default) nintersects = 40 meaning only top 40 intersections will be plotted
dev.off()
# Observation:
# RFX2-5 share the largest number of dys-reg ATAC peaks
# followed by JUN/FOS cluster 7 TFs (JUNB, JUN, FOSL2, FOS, JUND, FOSL1, BATF, JDP2, NFE2)
# followed by CPEB1 uniquely bound
# followed by JUN/FOS cluster 5 TFs (no JDP2 or NFE2)

## Conclusions:
# The representative TF in the JUN/FOS cluster would be FOSL1 (all bound dys-reg ATAC are shared with other cluster members)
# The representative TF in the RFX cluster would be RFX2 (all bound dys-reg ATAC are shared with other cluster members). RFX5 has 5 uniquely bound, while RFX4 has 4 uniquely bound.
# No representative TF in the MEF2 cluster. Just use the 10+3 ATAC peaks co-bound by all 4 members. 
# No representative TF in the DBP cluster (DBP, HLF, TEF). Just use the 10+2 ATAC peaks co-bound by all 3 members and NFIL3. 

