
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/ind_BINDetect_Jaspar2022/")

load("Jaspar2022hs_TFmotifs_mean_score_by_sample.rda") # df
colnames(df) = sub("_mean.score", "", colnames(df))
colnames(df) = sub("\\.", "-", colnames(df))
colnames(df)[2:5] = sub("X", "", colnames(df)[2:5])

## Filter for the motifs of interest
lnames = load("DiffBS_lm_Volcano.rda") # DiffBS_lm_Volcano
df = df[df$motif_id %in% rownames(DiffBS_lm_Volcano)[DiffBS_lm_Volcano$TFlabel != ""],] # 45 TFmotifs
df$TFmotif = DiffBS_lm_Volcano$TFmotif[match(df$motif_id, rownames(DiffBS_lm_Volcano))]
df = df[,-1]

## Melt motif binding score by sample dataframe (df -> df_melt), and attach sample meta data (Age, Diagnosis)
df_melt = melt(df, variable.name = "sample_name", value.name = "mean_binding_score")

lnames = load("../../3_QC/BiolTechCov_and_QC.rda") # df_wAllCov
colnames(df_wAllCov)
datMeta = df_wAllCov[match(df_melt$sample_name, df_wAllCov$sample_name), c("sample_name", "Diagnosis", "Age", "Sex", "Cortex", "Region")]
df_melt = cbind(df_melt, datMeta[,-1])

## Plot
df_melt$change = DiffBS_lm_Volcano$col[match(df_melt$TFmotif, DiffBS_lm_Volcano$TFmotif)]
df_melt = df_melt %>%
  arrange(change, TFmotif)
#df_melt$TFmotif = factor(df_melt$TFmotif, levels = rev(levels(df_melt$TFmotif)))
df_melt$TFmotif = factor(df_melt$TFmotif, levels = unique(df_melt$TFmotif))

pdf("DiffTFbylme_MeanBindingScore_bySampleAge.pdf", width = 20, height = 12)
df_melt %>%
#pdf("DiffTFbylme_MeanBindingScore_bySampleAge_rmB5000.pdf", width = 20, height = 12)
#df_melt[! grepl("B5000", df_melt$sample_name),] %>%
  ggplot(aes(x = Age, y = mean_binding_score, col = Diagnosis)) +
  geom_point() +
  geom_smooth(method = "loess", span = 2, se = F) +
  theme_bw() +
  facet_wrap(~ TFmotif, scales = "free_y")
dev.off()
# Observation:
# TFs of the same family almost show exactly the same binding scores.
# For TF motif up in ASD, an ASD sample at age 27-ish (B5000B and B5000C) show extremely high TF binding score, affecting the smooth line.
# After removing B5000, the TF motifs up in ASD appear much less different from CTL. The TF motifs down in ASD are not affected.

# Cluster by TF family
df_melt$cluster = DiffBS_lm_Volcano$cluster[match(df_melt$TFmotif, DiffBS_lm_Volcano$TFmotif)]

TFclusters = unique(df_melt[,c("TFmotif", "cluster", "change")])
TFclusters = TFclusters %>%
  arrange(change, cluster)
TFclusters$TF_shortname = DiffBS_lm_Volcano$TF_shortname[match(TFclusters$TFmotif, DiffBS_lm_Volcano$TFmotif)]

TFclusters$TFs = TFclusters$TF_shortname
for (i in 1:nrow(TFclusters)) {
  c = TFclusters$cluster[i]
  idx = which(TFclusters$cluster == c)
  if (length(idx) > 1) {
    TFclusters$TFs[i] = paste0(unique(TFclusters$TF_shortname[idx]), collapse = "/")
  }
}

df_melt$TFs = TFclusters$TFs[match(df_melt$TFmotif, TFclusters$TFmotif)]
unique(df_melt$TFs)
df_melt$TFs[which(df_melt$TFs == "BACH1/BACH2/BATF/BATF::JUN/BATF3/BNC2/FOS/FOS::JUN/FOS::JUNB/FOS::JUND/FOSB::JUNB/FOSL1/FOSL1::JUN/FOSL1::JUNB/FOSL1::JUND/FOSL2/FOSL2::JUN/FOSL2::JUNB/FOSL2::JUND/JDP2/JUN::JUNB/JUNB/JUND/NFE2")] = "BACH(1/2)/BATF(3)/BNC2/FOS(L1/2)/JUN(B/D)/JDP2/NFE2"
df_melt$TFs[which(df_melt$TFs == "MEF2A/MEF2B/MEF2C/MEF2D")] = "MEF2A/B/C/D"
df_melt$TFs[which(df_melt$TFs == "BHLHA15/BHLHE23/OLIG1/OLIG2/OLIG3")] = "BHLHA15/BHLHE23/OLIG(1/2/3)"

idx_dup = which(duplicated(TFclusters[,c("cluster", "TFs")]))
TF_dups = TFclusters$TFmotif[idx_dup] # 33
df_melt2 = df_melt[-which(df_melt$TFmotif %in% TF_dups),]
df_melt2 = df_melt2 %>%
  arrange(change, TFs)
df_melt2$TFs = factor(df_melt2$TFs, levels = unique(df_melt2$TFs))

pdf("DiffTFbylme_MeanBindingScore_bySampleAge_ReducedToTFCluster.pdf", width = 12, height = 6)
df_melt2 %>%
#pdf("DiffTFbylme_MeanBindingScore_bySampleAge_rmB5000_ReducedToTFCluster.pdf", width = 12, height = 6)
#df_melt2[! grepl("B5000", df_melt2$sample_name),] %>%
  ggplot(aes(x = Age, y = mean_binding_score, col = Diagnosis)) +
  geom_point() +
  geom_smooth(method = "loess", span = 2, se = F) +
  theme_bw() +
  facet_wrap(~ TFs, scales = "free_y")
dev.off()
# Observations:
# Yes, much cleaner to view, bu clustering the TFs
# CTCF, YY2, ZBTB33 show lower binding score in ASD, while FOS/JUN/etc, BHLH/OLIG, DBP, MEF2, NEUROG1, ZNF211 show higher binding score in ASD.

## Check TF gene expression
library(readxl)
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

DiffBS_lm_Volcano$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]
DiffBS_lm_Volcano$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]

DiffBS_lm_Volcano$WholeCortex_dup15q_logFC = Jill_DEG$WholeCortex_dup15q_logFC[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]
DiffBS_lm_Volcano$WholeCortex_dup15q_FDR = Jill_DEG$WholeCortex_dup15q_FDR[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]

TFclusters = left_join(TFclusters, DiffBS_lm_Volcano[,c("TFmotif", "magnitude", "p", "sig", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR", "WholeCortex_dup15q_logFC", "WholeCortex_dup15q_FDR")])

plot(TFclusters$WholeCortex_ASD_logFC, TFclusters$magnitude, xlab = "Differential gene expression", ylab = "Differential binding score", pch = 19, col = alpha("black", 0.5))
plot(TFclusters$WholeCortex_dup15q_logFC, TFclusters$magnitude, xlab = "Differential gene expression", ylab = "Differential binding score", pch = 19, col = alpha("black", 0.5))
TFclusters$GE_sig = ifelse(TFclusters$WholeCortex_ASD_FDR < 0.05, "significant", "n.s.") # use FDR < 0.05 rather than 0.1! If 0.1, one more point in Q2.
TFclusters$GE_sig[is.na(TFclusters$GE_sig)] = "n.s."
TFclusters$GE_sig = factor(TFclusters$GE_sig, levels = c("significant", "n.s."))
TFclusters$GE_sig_dup15q = ifelse(TFclusters$WholeCortex_dup15q_FDR < 0.05, "significant", "n.s.") # use FDR < 0.05 rather than 0.1! If 0.1, one more point in Q2.
TFclusters$GE_sig_dup15q[is.na(TFclusters$GE_sig_dup15q)] = "n.s."
TFclusters$GE_sig_dup15q = factor(TFclusters$GE_sig_dup15q, levels = c("significant", "n.s."))

TFclusters$GE_label = ifelse(TFclusters$GE_sig == "significant" & TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0, TFclusters$TF_shortname, "")
TFclusters$GE_label = ifelse(TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0, TFclusters$TF_shortname, "")
TFclusters$GE_label_dup15q = ifelse(TFclusters$GE_sig_dup15q == "significant" & TFclusters$magnitude * TFclusters$WholeCortex_dup15q_logFC > 0, TFclusters$TF_shortname, "")
TFclusters$GE_label_dup15q = ifelse(TFclusters$magnitude * TFclusters$WholeCortex_dup15q_logFC > 0, TFclusters$TF_shortname, "")

# Split xx::xx to get DGE of each TF [to include JUN, optional]
# TFclusters1 = TFclusters[-which(grepl("::", TFclusters$TF_shortname)),]
# TFclusters3 = TFclusters2 = TFclusters[grepl("::", TFclusters$TF_shortname),]
# 
# TFclusters1$Gene_name = TFclusters1$TF_shortname
# TFclusters2$Gene_name = sapply(TFclusters2$TF_shortname, function(x) unlist(str_split(x, "::"))[1])
# TFclusters3$Gene_name = sapply(TFclusters2$TF_shortname, function(x) unlist(str_split(x, "::"))[2])
# 
# TFclusters2$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(TFclusters2$Gene_name, Jill_DEG$external_gene_name)]
# TFclusters2$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(TFclusters2$Gene_name, Jill_DEG$external_gene_name)]
# TFclusters3$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(TFclusters3$Gene_name, Jill_DEG$external_gene_name)]
# TFclusters3$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(TFclusters3$Gene_name, Jill_DEG$external_gene_name)]
# 
# TFclusters = rbind(TFclusters1, TFclusters2, TFclusters3)
# TFclusters$GE_sig = ifelse(TFclusters$WholeCortex_ASD_FDR < 0.05, "significant", "n.s.") # use FDR < 0.05 rather than 0.1! 
# TFclusters$GE_sig = factor(TFclusters$GE_sig, levels = c("significant", "n.s."))
# TFclusters$GE_label = ifelse(TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0, TFclusters$Gene_name, "")
# tmp = TFclusters[,c("Gene_name", "magnitude", "WholeCortex_ASD_logFC")]
# tmp$magnitude = round(tmp$magnitude, 4)
# idx = which(duplicated(tmp))
# TFclusters$GE_label[idx] = ""
# DiffTFbylme_DGE_scatterplot_includeJUN.pdf

# scatter plot
#pdf("DiffTFbylme_DGE_scatterplot_includeJUN.pdf", width = 8, height = 5)
pdf("DiffTFbylme_DGE_scatterplot.pdf", width = 8, height = 5)
#pdf("DiffTFbylme_DGE_scatterplot_allTFnameLabeled.pdf", width = 8, height = 5)
TFclusters %>%
  ggplot(aes(x = WholeCortex_ASD_logFC, y = magnitude, col = GE_sig, label = GE_label)) +
  #ggplot(aes(x = WholeCortex_ASD_logFC, y = magnitude, col = GE_sig, label = TF_shortname)) + # _allTFnameLabeled
  geom_point(alpha = 0.7) +
  geom_text_repel(box.padding = 0.6, max.overlaps = Inf) +
  labs(col = "TF differential\ngene expression") +
  theme_bw() +
  xlab("LogFC of differential TF gene expression (ASD vs CTL)") +
  ylab("Difference in TF binding score (ASD vs CTL)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
dev.off()
# Observation:
# TFs significantly change GE in the same dir as biding score:
# a. both up: BACH1, NFIL3, NFE2, BATF, FOSL1 (No JUN?)
# b. both down: None
# Note that CTCF expression dropped, though n.s. (FDR> 0.1)

pdf("DiffTFbylme_DGE_scatterplot_Dup15q.pdf", width = 8, height = 5)
#pdf("DiffTFbylme_DGE_scatterplot_allTFnameLabeled.pdf", width = 8, height = 5)
TFclusters %>%
  ggplot(aes(x = WholeCortex_dup15q_logFC, y = magnitude, col = GE_sig_dup15q, label = GE_label_dup15q)) +
  #ggplot(aes(x = WholeCortex_ASD_logFC, y = magnitude, col = GE_sig, label = TF_shortname)) + # _allTFnameLabeled
  geom_point(alpha = 0.7) +
  geom_text_repel(box.padding = 0.6, max.overlaps = Inf) +
  labs(col = "TF differential\ngene expression") +
  theme_bw() +
  xlab("LogFC of differential TF gene expression (Dup15q vs CTL)") +
  ylab("Difference in TF binding score (ASD vs CTL)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
dev.off()

save(DiffBS_lm_Volcano, TFclusters, file = "DiffBS_lm_Volcano_wCluster_DGE_ASDorDup15q.rda")

