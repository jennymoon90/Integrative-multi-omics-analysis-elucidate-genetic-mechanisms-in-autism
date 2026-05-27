# Hypothesis:
# In ASD, the correlation between ATAC peaks at Hi-C loop ends drops, in accompany to drop of Hi-C contact, especially in neurons.

# Steps:
# 1) Confirm known findings: ATAC pairs linked by Hi-C loops have significantly higher correlation
# 2) In ASD, such correlation drops
# 3) The drop of correlation mainly happens at down-regulated Hi-C loops
# 4) pSI/GO enrichment of the genes at these down-regulated Hi-C loops

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)
library(ggplot2)
library(reshape2)
library(plyr)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/6_4_v1_ATACcorAtHiCloopends/")

### 1) ATAC peaks at Hi-C loop ends show significantly higher correlation
## load CQN normalized ATAC peaks
lnames = load("../../5_DiffATAC/CQN.rda") # PRKM.cqn

## Load Bulk Hi-C consensus loops
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/PL_by_sample_datasets/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_rmOutlierTechRep.rda") # Bulk_logCPM, NeuNp_logCPM, NeuNn_logCPM

rm(list = setdiff(ls(), c("RPKM.cqn", "Bulk_logCPM")))

## Assign ATAC peaks to 10kb bins
ATAC_location = as.data.frame(str_split_fixed(rownames(RPKM.cqn), "_", 3))
ATAC_location$V2 = as.integer(ATAC_location$V2);ATAC_location$V3 = as.integer(ATAC_location$V3)
ATAC_location$mid = (ATAC_location$V2 + ATAC_location$V3)/2
ATAC_location$fragmentMid = floor(ATAC_location$mid/10e3)*10e3 + 5e3
ATAC_location$mid = paste0(ATAC_location$V1, "_", ATAC_location$fragmentMid)

## ATAC pairs linked by Hi-C
Loop_location = as.data.frame(str_split_fixed(rownames(Bulk_logCPM), "_", 3))
Loop_location$V2 = as.integer(Loop_location$V2);Loop_location$V3 = as.integer(Loop_location$V3)
Loop_location$mid1 = paste0(Loop_location$V1, "_", Loop_location$V2)
Loop_location$mid2 = paste0(Loop_location$V1, "_", Loop_location$V3)

# Filter for loops >= 20kb
Loop_location$distance = Loop_location$V3 - Loop_location$V2
range(Loop_location$distance) # 20kb - 4.9Mb, already filtered

ATAC_location$ATACpeak = rownames(RPKM.cqn) # to match back to RPKM.cqn
Loop_location$Loop = rownames(Bulk_logCPM)

ATAC1 = ATAC2 = ATAC_location
colnames(ATAC1)[4] = "mid1"; colnames(ATAC2)[4] = "mid2"
colnames(ATAC1)[6] = "ATACpeak1"; colnames(ATAC2)[6] = "ATACpeak2"

ATAC_pairs = left_join(Loop_location[,c("Loop", "distance", "mid1", "mid2")], ATAC1[,c("mid1", "ATACpeak1")]) # 49065 rows
ATAC_pairs = left_join(ATAC_pairs, ATAC2[,c("mid2", "ATACpeak2")]) # 70767 rows
ATAC_pairs = unique(ATAC_pairs[complete.cases(ATAC_pairs),]) # 58462 rows
# Total 58462 ATAC pairs linked by Hi-C

## Randomly select equal number of ATAC pairs on the same chromosome but not linked by Hi-C
ATAC_pairs_nonHiC = matrix(ncol = 2, nrow = 0)
colnames(ATAC_pairs_nonHiC) = c("ATACpeak1", "ATACpeak2")
n = 1
set.seed(123)
while (n <= nrow(ATAC_pairs)) {
  if (n %% 5000 == 0) {print(n)}
  idx1 = sample(1:nrow(ATAC_location), 1)
  ATACpeak_idx1 = ATAC_location$ATACpeak[idx1]
  ATACpeak_paired = c(ATAC_pairs$ATACpeak2[which(ATAC_pairs$ATACpeak1 == ATACpeak_idx1)], ATAC_pairs$ATACpeak1[which(ATAC_pairs$ATACpeak2 == ATACpeak_idx1)])
  pool = which(ATAC_location$V1 == ATAC_location$V1[idx1] & ! ATAC_location$ATACpeak %in% c(ATACpeak_idx1, ATACpeak_paired))
  idx2 = sample(pool, 1)
  
  ATAC_pairs_nonHiC = rbind(ATAC_pairs_nonHiC, c(ATAC_location$ATACpeak[c(min(idx1, idx2), max(idx1, idx2))])) # order the ATAC pairs, so that it's unique
  ATAC_pairs_nonHiC = unique(ATAC_pairs_nonHiC)
  n = nrow(ATAC_pairs_nonHiC)
}
# 2022/8/19 3:39-4:29 PM
# This is taking too long, next time modify the script to find a non-Hi-C linked ATAC peak on the same chromosome for each ATAC peak, instead of sampling one by one and unique()

ATAC_pairs_nonHiC = as.data.frame(ATAC_pairs_nonHiC)
save(ATAC_pairs, ATAC_pairs_nonHiC, file = "st1_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiC_vs_NotLinked.rda")

## Calculate correlation of ATAC pairs 
colnames(ATAC_pairs)[5:6] # "ATACpeak1" "ATACpeak2"
ATAC_pairs$cor = apply(ATAC_pairs, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[5]),], RPKM.cqn[which(rownames(RPKM.cqn) == x[6]),])) # took a full minute.

colnames(ATAC_pairs_nonHiC) # "ATACpeak1" "ATACpeak2"
ATAC_pairs_nonHiC$cor = apply(ATAC_pairs_nonHiC, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[1]),], RPKM.cqn[which(rownames(RPKM.cqn) == x[2]),])) # took a full minute.

## Violin Plot
ATAC_pairs$type = "Linked_by_HiC"
ATAC_pairs_nonHiC$type = "Not_linked_by_HiC"
df = rbind(ATAC_pairs[,c("ATACpeak1", "ATACpeak2", "cor", "type")], ATAC_pairs_nonHiC[,c("ATACpeak1", "ATACpeak2", "cor", "type")])

save(df, ATAC_pairs, ATAC_pairs_nonHiC, file = "st1_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiC_vs_NotLinked.rda")

# Wilcox test
unique(df$type)
df$type = factor(df$type)
(wx_res = wilcox.test(df$cor[df$type == "Linked_by_HiC"], df$cor[df$type == "Not_linked_by_HiC"]))
p = wx_res$p.value # p<2.2e-16

pdf("st1_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiC_vs_NotLinked.pdf", width = 3, height = 4)
df %>%
  ggplot(aes(x = type, y = cor)) +
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  ylab("Correlation of ATAC pairs") +
  ylim(min(df$cor), 1.1) +
  annotate("text", x = 1.5, y = 1.1, label = "p < 2.2e-16")
dev.off()

# Observation: obvious and significant difference! GREAT!

### 2) In ASD, does such correlation drop? - No
lnames = load("../../5_DiffATAC/CQN.rda") # PRKM.cqn
rm(list = setdiff(ls(), c("RPKM.cqn", "Covariates")))
lnames = load("st1_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiC_vs_NotLinked.rda")

## Calculate correlation of ATAC pairs in ASD vs. CTL
idx_ASD = which(colnames(RPKM.cqn) %in% Covariates$sample_name[Covariates$Diagnosis == "ASD"]) # 16
idx_CTL = which(colnames(RPKM.cqn) %in% Covariates$sample_name[Covariates$Diagnosis == "CTL"]) # 22

colnames(ATAC_pairs)[5:6] # "ATACpeak1" "ATACpeak2"
ATAC_pairs$cor_ASD = apply(ATAC_pairs, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[5]),idx_ASD], RPKM.cqn[which(rownames(RPKM.cqn) == x[6]),idx_ASD])) # took a full minute.
ATAC_pairs$cor_CTL = apply(ATAC_pairs, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[5]),idx_CTL], RPKM.cqn[which(rownames(RPKM.cqn) == x[6]),idx_CTL])) # took a full minute.

colnames(ATAC_pairs_nonHiC)[1:2] # "ATACpeak1" "ATACpeak2"
ATAC_pairs_nonHiC$cor_ASD = apply(ATAC_pairs_nonHiC, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[1]),idx_ASD], RPKM.cqn[which(rownames(RPKM.cqn) == x[2]),idx_ASD])) # took a full minute.
ATAC_pairs_nonHiC$cor_CTL = apply(ATAC_pairs_nonHiC, 1, function(x) cor(RPKM.cqn[which(rownames(RPKM.cqn) == x[1]),idx_CTL], RPKM.cqn[which(rownames(RPKM.cqn) == x[2]),idx_CTL])) # took a full minute.

## Violin Plot
df = rbind(ATAC_pairs[,c("ATACpeak1", "ATACpeak2", "cor", "type", "cor_CTL", "cor_ASD")], ATAC_pairs_nonHiC[,c("ATACpeak1", "ATACpeak2", "cor", "type", "cor_CTL", "cor_ASD")])
df_melt = melt(df, id = c("ATACpeak1", "ATACpeak2", "cor", "type"), variable.name = "Diagnosis", value.name = "cor_Diagnosis")

save(df, df_melt, ATAC_pairs, ATAC_pairs_nonHiC, file = "st2_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiCvsNot_inASDvsCTL.rda")

# ANOVA test
aov_res = aov(cor_Diagnosis ~ type + Diagnosis, data = df_melt)
(aov_summary = summary(aov_res)) # Only type has significant effect (p<2e-16)

# Wilcoxon test
unique(df$type)
(wx_res = wilcox.test(df$cor_ASD[df$type == "Linked_by_HiC"], df$cor_CTL[df$type == "Linked_by_HiC"]))
p_linked = round(wx_res$p.value,2) # p=0.21, no significant change - makes sense as most loops did not change

wx_res = wilcox.test(df$cor_ASD[df$type == "Not_linked_by_HiC"], df$cor_CTL[df$type == "Not_linked_by_HiC"])
p_unlinked = round(wx_res$p.value,2) # p=0.01, significantly increased in ASD - may be false positive
mean(df$cor_ASD[df$type == "Not_linked_by_HiC"]) # 0.0493
mean(df$cor_CTL[df$type == "Not_linked_by_HiC"]) # 0.0434

df_melt$Diagnosis = gsub("cor_", "", df_melt$Diagnosis)
df_melt$Diagnosis = factor(df_melt$Diagnosis, levels = c("CTL", "ASD"))

pdf("st2_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiCvsNot_inASDvsCTL.pdf", width = 5, height = 4)
df_melt %>%
  ggplot(aes(x = type, y = cor_Diagnosis)) +
  geom_violin(aes(fill = Diagnosis)) +
  geom_boxplot(aes(fill = Diagnosis), width = 0.1, position=position_dodge(0.9)) +
  theme_bw() +
  ylab("Correlation of ATAC pairs") +
  ylim(min(df$cor), 1.22) +
  #annotate("text", x = 1, y = 1.1, label = paste0("p = ",p_linked)) +
  #annotate("text", x = 2, y = 1.1, label = paste0("p = ",p_unlinked))
  annotate("text", x = 1.5, y = 1.22, label = paste0("p < 2e-16")) +
  geom_segment(aes(x = 0.8, y = 1.05, xend = 1.2, yend = 1.05)) +
  geom_segment(aes(x = 1.8, y = 1.05, xend = 2.2, yend = 1.05)) +
  geom_segment(aes(x = 1, y = 1.15, xend = 2, yend = 1.15)) +
  geom_segment(aes(x = 1, y = 1.05, xend = 1, yend = 1.15)) +
  geom_segment(aes(x = 2, y = 1.05, xend = 2, yend = 1.15))
dev.off()

pdf("st2_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiC_inASDvsCTL.pdf", width = 3, height = 4)
df_melt[df_melt$type == "Linked_by_HiC",] %>%
  ggplot(aes(x = Diagnosis, y = cor_Diagnosis)) +
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  ylab("Correlation of ATAC pairs linked by Hi-C loops") +
  ylim(min(df$cor), 1.1) +
  annotate("text", x = 1.5, y = 1.1, label = paste0("p = ",p_linked))
dev.off()

# Observation: No significant change in correlation of ATAC pairs linked by HiC between ASD and CTL

## Get significantly correlated ATAC peaks
load("st2_Correlation_of_ATACpeaksCqnNormalizedLogRPKM_LinkedByHiCvsNot_inASDvsCTL.rda")
ATAC_pairs_nonHiC = ATAC_pairs_nonHiC %>%
  arrange(cor)
ATAC_pairs$cor_p = sapply(ATAC_pairs$cor, function(x) length(which(ATAC_pairs_nonHiC$cor > x))/nrow(ATAC_pairs_nonHiC))
ATAC_pairs$cor_fdr = p.adjust(ATAC_pairs$cor_p, method = "fdr")
length(which(ATAC_pairs$cor_fdr < 0.05)) # 42 pairs
length(which(ATAC_pairs$cor_fdr < 0.1)) # 112 pairs
length(which(ATAC_pairs$cor_p < 0.05)) # 10780 pairs
length(which(ATAC_pairs$cor_p < 0.01)) # 3774 pairs
mean(ATAC_pairs_nonHiC$cor) # 0.05
mean(ATAC_pairs$cor[which(ATAC_pairs$cor_fdr < 0.1)]) # 0.93
mean(ATAC_pairs$cor[which(ATAC_pairs$cor_p < 0.05)]) # 0.68 - decent
mean(ATAC_pairs$cor[which(ATAC_pairs$cor_p < 0.01)]) # 0.78 - decent

Correlated_ATAC_pairs = ATAC_pairs[which(ATAC_pairs$cor_p < 0.05), ] # 10780, but they do not necessarily link gene promoter with distal enhancers
save(df, df_melt, ATAC_pairs, ATAC_pairs_nonHiC, Correlated_ATAC_pairs, file = "st2_plus_CorrelatedATACpairsP005.rda")


  