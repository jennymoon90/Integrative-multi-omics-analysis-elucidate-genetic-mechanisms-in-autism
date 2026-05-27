rm(list = ls())
options(scipen = 0)

library(stringr)
library(tidyverse)
library(GenomicRanges)
library(Repitools)

#setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/04_mQTLs/")
setwd("/Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/04_mQTLs/")

# ---- Procedures ----
### (1) Make data frame of SNP genomic location, effect allele, effect size on probe methylation, effect size on ASD risk from GWAS, CpG island of the target probe, whether the CpG island is in my DMR_up/down/ns list, whether the CpG island overlaps TAD boundaries.
### (2) Plot the mQTL effect on DNA methylation vs. GWAS risk, for those at TAD boundaries
### (3) Plot the mQTL effect on DNA methylation vs. GWAS risk, highlighting 3 categories: DMRs at weakened TAD boundaries, at all TAD boundaries, all mQTLs
# --------------------

### (1) Make data frame of SNP genomic location, effect allele, effect size on probe methylation, effect size on ASD risk from GWAS, CpG island of the target probe, whether the CpG island is in my DMR_up/down/ns list, whether the CpG island overlaps TAD boundaries.

#lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/40_Hannon_2016_Nat Neurosci_mQTL fetal brain/Supplementary_Table_2_mQTLs_hg19.rda") # fetal_mQTLs
lnames = load("/Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Database_download/40_Hannon_2016_Nat Neurosci_mQTL fetal brain/Supplementary_Table_2_mQTLs_hg19.rda") # fetal_mQTLs

df = data_frame(SNP = fetal_mQTLs$SNP_ID, methyEffectAllele = fetal_mQTLs$Fetal_Allele, methyBeta = fetal_mQTLs$Fetal_Regression_coefficient, TargetProbe = fetal_mQTLs$Probe_ID)

## Add CpG island and DMR info
lnames = load("../02_redo_DMG_analysis/08_DMR_fdr005.rda") # DMR, HM450k_final, Probe_Island_grouped

df$TargetCpGisland = HM450k_final$UCSC_CpG_Islands_Name[match(df$TargetProbe, HM450k_final$Probe)]
df$DMR_Beta = DMR$magnitude[match(df$TargetCpGisland, rownames(DMR))]
df$DMR_FDR = DMR$fdr[match(df$TargetCpGisland, rownames(DMR))]

## Add TADb info
ASD_TAD_dir = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/6_TADs/IS_GENOVA/differential_IS/v4_afterGRC_TADbOnly_includeTransCisInLm/"
lnames = load(paste0(ASD_TAD_dir, "p_magnitude_Bulk_lmeDiagnosisAgeSexBatch2.TemporalABNNICHDvalid.interactioncis.longRange_includingTechRepsAndOutliers.rda")) # p, magnitude, diffIS

weakened_TADb = rownames(diffIS)[which(diffIS$FDR < 0.05 & diffIS$Magnitude > 0)] # 754 weakened TADb 
TADb_df = as.data.frame(str_split_fixed(rownames(diffIS), "_", 3))
TADb_df$V2 = as.integer(TADb_df$V2); TADb_df$V3 = as.integer(TADb_df$V3)
colnames(TADb_df) = c("chr", "start", "end")
TADb_gr = annoDF2GR(TADb_df)

TargetCpG = unique(df$TargetCpGisland)
TargetCpG = TargetCpG[!is.na(TargetCpG)] # 1068 unique CpG islands
TargetCpG_df = as.data.frame(str_split_fixed(TargetCpG, ":|-", 3))
TargetCpG_df$V2 = as.integer(TargetCpG_df$V2); TargetCpG_df$V3 = as.integer(TargetCpG_df$V3)
colnames(TargetCpG_df) = c("chr", "start", "end")
TargetCpG_gr = annoDF2GR(TargetCpG_df)

hits = as.data.frame(findOverlaps(TargetCpG_gr, TADb_gr))
hits_df = data_frame(TargetCpGisland = TargetCpG[hits$queryHits], TADb = rownames(diffIS)[hits$subjectHits], TADb_Magnitude = diffIS$Magnitude[hits$subjectHits], TADb_FDR = diffIS$FDR[hits$subjectHits])

df$TADb = hits_df$TADb[match(df$TargetCpGisland, hits_df$TargetCpGisland)]
df$TADb_Magnitude = hits_df$TADb_Magnitude[match(df$TargetCpGisland, hits_df$TargetCpGisland)]
df$TADb_FDR = hits_df$TADb_FDR[match(df$TargetCpGisland, hits_df$TargetCpGisland)]

save(df, file = "2_01_FetalBrain_mQTL_DMR_TADb_dataframe.rda")

## Load GWAS info (Grove 2019 ASD)
ASD_GWAS = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/GWAS_summary_statistics/ASD_Grove_2018/ASD_Grove_2019.sumstats.gz", header = T)

df_wGWAS = left_join(df, ASD_GWAS[,1:4]) # 16809 rows
df_wGWAS = df_wGWAS[complete.cases(df_wGWAS[,which(colnames(df_wGWAS) %in% c("A1", "A2"))]),] # 15910 rows
# Not all methyEffectAllele are in A1 or A2. For instance, the first SNP, rs645279, show up as T>A,C,G in dbSNP. The methyEffectAllele is one of the three minor alleles, A1 is T (the ref allele) and A2 is C (one of the other minor alleles). So this won't be useful. -> I shall delete all the SNPs that have methyEffectAllele not in A1 or A2. -> But if the methyEffectAllele is in A1 and A2, how do I know the non-Effect allele is in A1 or A2? No way.
idx = which(apply(df_wGWAS[,c("methyEffectAllele", "A1", "A2")], 1, function(x) x[1] %in% c(x[2], x[3]))) # 8057 SNPs left

df_wGWAS_validAlleles = df_wGWAS[idx,]

## Change methyBeta to the same direction as DMR_Beta.
idx = which(df_wGWAS_validAlleles$DMR_Beta * df_wGWAS_validAlleles$methyBeta < 0) # 1461
for (i in idx) {
  if (df_wGWAS_validAlleles$methyEffectAllele[i] == df_wGWAS_validAlleles$A1[i]) {
    df_wGWAS_validAlleles$methyEffectAllele[i] = df_wGWAS_validAlleles$A2[i]
    df_wGWAS_validAlleles$methyBeta[i] = -df_wGWAS_validAlleles$methyBeta[i]
  } else if (df_wGWAS_validAlleles$methyEffectAllele[i] == df_wGWAS_validAlleles$A2[i]) {
    df_wGWAS_validAlleles$methyEffectAllele[i] = df_wGWAS_validAlleles$A1[i]
    df_wGWAS_validAlleles$methyBeta[i] = -df_wGWAS_validAlleles$methyBeta[i]
  } else {
   print(i) 
  }
}

## Match A1 to the same as methyEffectAllele
idx = which(df_wGWAS_validAlleles$A1 != df_wGWAS_validAlleles$methyEffectAllele) # 3821 SNPs
df_wGWAS_alladjusted = df_wGWAS_validAlleles
df_wGWAS_alladjusted$A2[idx] = df_wGWAS_alladjusted$A1[idx]
df_wGWAS_alladjusted$A1[idx] = df_wGWAS_alladjusted$methyEffectAllele[idx]
df_wGWAS_alladjusted$Z[idx] = -df_wGWAS_alladjusted$Z[idx]

save(df_wGWAS_alladjusted, file = "2_01_FetalBrain_mQTL_DMR_TADb_GWAS_alladjusted.rda")

### (2) Plot the mQTL effect on DNA methylation vs. GWAS risk

rm(list = ls())
lnames = load("2_01_FetalBrain_mQTL_DMR_TADb_GWAS_alladjusted.rda") # df_wGWAS_alladjusted
all(df_wGWAS_alladjusted$methyEffectAllele == df_wGWAS_alladjusted$A1) # T

## all SNPs at TAD boundaries
idx = which(!is.na(df_wGWAS_alladjusted$TADb)) # 781 SNPs
plot(df_wGWAS_alladjusted$methyBeta[idx], df_wGWAS_alladjusted$Z[idx])
# looks good, hyper-methylation increases ASD risk effect

res_TADb = summary(lm(df_wGWAS_alladjusted$Z[idx] ~ 0 + df_wGWAS_alladjusted$methyBeta[idx]))
lm_coef_TADb = round(res_TADb$coefficients[1], 1) # 3
lm_p_TADb = signif(res_TADb$coefficients[4], 1) # 0.003
# Great, significant positive slope! But I don't have a way to determine which effect allele to use for mQTLs that are not associated with any CpG island. That's not true. Because I found TAD boundaries based on CpG islands. So if an mQTL does not associate with a CpG islands, it won't be associated with a TAD boundary. Yay!
cor(df_wGWAS_alladjusted$Z[idx], df_wGWAS_alladjusted$methyBeta[idx]) # 0.1055

pdf("2_02_mQTLs_atAllTADb_effect_on_MethyCpGislands_ASDrisk.pdf", height = 5, width = 8)
df_wGWAS_alladjusted[idx,] %>%
  ggplot(aes(methyBeta, Z)) +
  geom_point() +
  geom_abline(slope = lm_coef_TADb, intercept = 0, col = "red") +
  theme_bw() +
  xlab("Effect on DNA methylation") +
  ylab("Effect on ASD risk") +
  annotate("text", x = 0.25, y = 4, label = "DNA hyper-\nmethylation in ASD") +
  annotate("segment", x = 0.22, xend = 0.24, y = 2.9, yend = 3.2, arrow = arrow(length = unit(2,"mm")), size = 1) +
  #annotate("text", x = -0.25, y = 4, label = paste0("cor = ", cor_all, "\np = ", cor_p_all), col = "red") # cor is too little
  annotate("text", x = -0.25, y = 4, label = paste0("slope = ", lm_coef_TADb, "\np = ", lm_p_TADb), col = "red") +
  ggtitle("mQTLs at all TAD boundaries")
dev.off()

save(df_wGWAS_alladjusted, lm_coef_TADb, lm_p_TADb, file = "2_02_mQTLs_atAllTADb_effect_on_MethyCpGislands_ASDrisk.rda")

### (3) Plot the mQTL effect on DNA methylation vs. GWAS risk - 3 categories

rm(list = ls())

## Load previous lm results for the mQTLs
lnames = load("03_mQTL_effect_on_MethyCpGislands_ASDrisk.rda") # df_wGWAS_alladjusted, lm_coef_all, lm_coef_TADb, lm_p_all, lm_p_TADb, cor_res_all, cor_res_TADb
lm_coef_DMR = lm_coef_TADb
lm_p_DMR = lm_p_TADb
rm(df_wGWAS_alladjusted,cor_res_all,cor_res_TADb,lm_coef_TADb,lm_p_TADb)

## Load current dataset and define the 3 mQTL categories
lnames = load("2_02_mQTLs_atAllTADb_effect_on_MethyCpGislands_ASDrisk.rda")

df_wGWAS_alladjusted$Category = ifelse(df_wGWAS_alladjusted$TADb_FDR < 0.05 & df_wGWAS_alladjusted$TADb_Magnitude > 0 & df_wGWAS_alladjusted$DMR_FDR < 0.05, "DMRs at weakened TAD boundaries", ifelse(!is.na(df_wGWAS_alladjusted$TADb), "At TAD boundaries", "Not at TAD boundaries"))
df_wGWAS_alladjusted$Category = factor(df_wGWAS_alladjusted$Category, levels = c("DMRs at weakened TAD boundaries", "At TAD boundaries", "Not at TAD boundaries"))
df_wGWAS_alladjusted$Category[is.na(df_wGWAS_alladjusted$Category)] = "Not at TAD boundaries"
table(df_wGWAS_alladjusted$Category) # 46 DMRs xxx, 735 At TADb, 2098->7276 Not at TADb. Consistent with previous results.
df_wGWAS_alladjusted = df_wGWAS_alladjusted[order(df_wGWAS_alladjusted$Category, decreasing = T),]

pdf("2_03_mQTL_3Categories_effect_on_MethyCpGislands_ASDrisk.pdf", height = 5, width = 8)
df_wGWAS_alladjusted %>% # [!is.na(df_wGWAS_alladjusted$DMR_Beta),]
  ggplot(aes(methyBeta, Z, col = Category)) +
  geom_point() +
  geom_abline(slope = lm_coef_all, intercept = 0, col = "grey20") +
  geom_abline(slope = lm_coef_TADb, intercept = 0, col = "dodgerblue") +
  geom_abline(slope = lm_coef_DMR, intercept = 0, col = "red") + # #F8766D too faint
  theme_bw() +
  xlab("Effect on DNA methylation") +
  ylab("Effect on ASD risk") +
  #labs(col = "Category") +
  scale_color_manual(values = c("red", "dodgerblue", "grey")) +
  annotate("text", x = 0.25, y = 4, label = "DNA hyper-\nmethylation in ASD") +
  annotate("segment", x = 0.22, xend = 0.24, y = 2.8, yend = 3.2, arrow = arrow(length = unit(2,"mm")), size = 1) +
  annotate("text", x = -0.25, y = 1, label = paste0("slope = ", lm_coef_all, "\np = ", lm_p_all), col = "grey20") +
  annotate("text", x = -0.25, y = 2.5, label = paste0("slope = ", lm_coef_TADb, "\np = ", lm_p_TADb), col = "dodgerblue") +
  annotate("text", x = -0.25, y = 4, label = paste0("slope = ", lm_coef_DMR, "\np = ", lm_p_DMR), col = "red")  # lm_p is significant
dev.off()

save(df_wGWAS_alladjusted, lm_coef_TADb, lm_p_TADb, lm_coef_all, lm_p_all, lm_coef_DMR, lm_p_DMR, file = "2_03_mQTL_3Categories_effect_on_MethyCpGislands_ASDrisk.rda")

## Do mQTLs that increase DNA methylation tend to weaken TAD boundaries in general? Yes, based on the dodgerblue dots (TADb are attached via CpG islands associated with the mQTLs, so the effect alleles are adjusted already.)
