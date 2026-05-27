rm(list = ls())
options(stringsAsFactors = F)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_10_OverlapATACwTFbound/")

## Load Promoter and Distal ATAC-GE pairs
lnames = load("../10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda") # "sig_CorATAC"      "sig_CorATAC_uniq" "promATAC"

## Load TFoi bound regions
lnames = load("../10_TOBIAS/ind_BINDetect_Jaspar2022/TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda") # TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF
TFoi = TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF
TFoi = TFoi[-which(grepl("NFIL3", TFoi))]

tf_bound = vector("list")
for (tf in c("BATF", "FOSL1", "BACH1", "NFE2")) {
  tfoi = TFoi[grepl(tf, TFoi)]
  bound = read.table(paste0("../10_TOBIAS/ind_BINDetect_Jaspar2022/TFoi_bound_ATACpeaks/",tfoi, "_BoundInAnySample.txt"), header = T)
  bound$peak = paste0(bound$peak_chr, "_", bound$peak_start, "_", bound$peak_end)
  #assign(paste0(tf,"_bound"), unique(bound$peak))
  tf_bound[[tf]] = unique(bound$peak)
}

tf_bound2 = vector("list")
tfoi = TFoi[grepl("CTCF", TFoi)]
for (tfoi_cur in tfoi) {
  bound = read.table(paste0("../10_TOBIAS/ind_BINDetect_Jaspar2022/TFoi_bound_ATACpeaks/",tfoi_cur, "_BoundInAnySample.txt"), header = T)
  bound$peak = paste0(bound$peak_chr, "_", bound$peak_start, "_", bound$peak_end)
  #assign(paste0(tf,"_bound"), unique(bound$peak))
  tf_bound2[[tfoi_cur]] = unique(bound$peak)
}

## Label TF-bound in the ATAC-GE pair dataframe
promATAC$BATF = ifelse(promATAC$ATACpeak %in% tf_bound[["BATF"]], "yes", "no")
promATAC$FOSL1 = ifelse(promATAC$ATACpeak %in% tf_bound[["FOSL1"]], "yes", "no")
promATAC$BACH1 = ifelse(promATAC$ATACpeak %in% tf_bound[["BACH1"]], "yes", "no")
promATAC$NFE2 = ifelse(promATAC$ATACpeak %in% tf_bound[["NFE2"]], "yes", "no")
promATAC$CTCF_MA0139.1 = ifelse(promATAC$ATACpeak %in% tf_bound2[["CTCF-MA0139.1"]], "yes", "no")
promATAC$CTCF_MA1929.1 = ifelse(promATAC$ATACpeak %in% tf_bound2[["CTCF-MA1929.1"]], "yes", "no")
promATAC$CTCF_MA1930.1 = ifelse(promATAC$ATACpeak %in% tf_bound2[["CTCF-MA1930.1"]], "yes", "no")

sig_CorATAC_uniq$BATF = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound[["BATF"]], "yes", "no")
sig_CorATAC_uniq$FOSL1 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound[["FOSL1"]], "yes", "no")
sig_CorATAC_uniq$BACH1 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound[["BACH1"]], "yes", "no")
sig_CorATAC_uniq$NFE2 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound[["NFE2"]], "yes", "no")
sig_CorATAC_uniq$CTCF_MA0139.1 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound2[["CTCF-MA0139.1"]], "yes", "no")
sig_CorATAC_uniq$CTCF_MA1929.1 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound2[["CTCF-MA1929.1"]], "yes", "no")
sig_CorATAC_uniq$CTCF_MA1930.1 = ifelse(sig_CorATAC_uniq$distal_ATAC %in% tf_bound2[["CTCF-MA1930.1"]], "yes", "no")

pdf("DiffATAC_DEG_logFC_PromoterATAC_eaTFbound.pdf", height = 4, width = 6)
idx_all = which(promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05)
for (i in 8:14) {
  tf = colnames(promATAC)[i]
  idx_plot = which(promATAC[,i] == "yes" & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05)
  
  lm_res_prom = summary(lm(WholeCortex_ASD_logFC ~ 0 + ATAC_logFC, data = promATAC[idx_plot, ]))
  lm_coef_prom = lm_res_prom$coefficients[1]
  lm_p_prom = lm_res_prom$coefficients[4]
  lm_r2_prom = lm_res_prom$r.squared
  
  p_cur = promATAC[idx_plot, ] %>%
    ggplot(aes(x = ATAC_logFC, y = WholeCortex_ASD_logFC)) +
    geom_point() + # col = "grey"
    geom_smooth(method='lm', formula= y~x) +
    #geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
    theme_bw() + 
    ggtitle(paste0("Promoter ATAC bound by ", tf)) +
    xlim(range(promATAC$ATAC_logFC[idx_all])[1], range(promATAC$ATAC_logFC[idx_all])[2]) +
    ylim(range(promATAC$WholeCortex_ASD_logFC[idx_all])[1], range(promATAC$WholeCortex_ASD_logFC[idx_all])[2]) +
    annotate("text", x = 0, y = range(promATAC$WholeCortex_ASD_logFC[idx_plot], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2_prom),2), "\n", ifelse(lm_p_prom < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p_prom, format = "e", digits = 1)))), col = "red")
  
  print(p_cur)
}
dev.off()
# Observations:
# All TFoi-bound promoter ATAC peaks show significant correlation with DEG logFC.
# Many differential promoter ATAC are co-bound by BATF, FOSL1, and BACH1.
# NFE2 show the highest correlation (cor = 0.41, p = 3e-3)
# Many more differential promoter ATAC are bound by CTCF, cor = 0.17 ~ 0.21. But CTCF may play a bigger role in distal ATAC-GE logFC cor

pdf("DiffATAC_DEG_logFC_DistalATAC_eaTFbound.pdf", height = 4, width = 6)
idx_all = which(sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05)
for (i in 7:13) {
  tf = colnames(sig_CorATAC_uniq)[i]
  idx_plot = which(sig_CorATAC_uniq[,i] == "yes" & sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05)
  
  lm_res = summary(lm(WholeCortex_ASD_logFC ~ 0 + distalATAC_logFC, data = sig_CorATAC_uniq[idx_plot, ]))
  lm_coef = lm_res$coefficients[1]
  lm_p = lm_res$coefficients[4]
  lm_r2 = lm_res$r.squared
  
  p_cur = sig_CorATAC_uniq[idx_plot, ] %>%
    ggplot(aes(x = distalATAC_logFC, y = WholeCortex_ASD_logFC)) +
    geom_point() + # col = "grey"
    geom_smooth(method='lm', formula= y~x) +
    #geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
    theme_bw() + 
    ggtitle(paste0("Distal ATAC bound by ", tf)) +
    xlim(range(sig_CorATAC_uniq$distalATAC_logFC[idx_all])[1], range(sig_CorATAC_uniq$distalATAC_logFC[idx_all])[2]) +
    ylim(range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_all])[1], range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_all])[2]) +
    annotate("text", x = 0, y = range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_plot], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2),2), "\n", ifelse(lm_p < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p, format = "e", digits = 1)))), col = "red")
  
  print(p_cur)
}
dev.off()
# Observations:
# No significant correlation between ATAC-DEG logFC for BATF/FOSL1/BACH1/NFE2-bound distal ATAC peaks. Only CTCF-bound distal ATAC show significant correlation with DEG logFC.
# CTCF-bound distal ATAC get significant and higher logFC cor (cor = 0.24 ~ 0.31) than all distal ATAC or CTCF-bound promoter ATAC

# What if I do not subset for diffATAC or DEG?
pdf("ATAC_GE_logFC_PromoterATAC_eaTFbound.pdf", height = 4, width = 6)
for (i in 8:14) {
  tf = colnames(promATAC)[i]
  idx_plot = which(promATAC[,i] == "yes")
  
  lm_res_prom = summary(lm(WholeCortex_ASD_logFC ~ 0 + ATAC_logFC, data = promATAC[idx_plot, ]))
  lm_coef_prom = lm_res_prom$coefficients[1]
  lm_p_prom = lm_res_prom$coefficients[4]
  lm_r2_prom = lm_res_prom$r.squared
  
  p_cur = promATAC[idx_plot, ] %>%
    ggplot(aes(x = ATAC_logFC, y = WholeCortex_ASD_logFC)) +
    geom_point() + # col = "grey"
    geom_smooth(method='lm', formula= y~x) +
    #geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
    theme_bw() + 
    ggtitle(paste0("Promoter ATAC bound by ", tf)) +
    #xlim(range(promATAC$ATAC_logFC, na.rm = T)[1], range(promATAC$ATAC_logFC, na.rm = T)[2]) +
    ylim(range(promATAC$WholeCortex_ASD_logFC, na.rm = T)[1], range(promATAC$WholeCortex_ASD_logFC, na.rm = T)[2]) +
    annotate("text", x = range(promATAC$ATAC_logFC[idx_plot], na.rm = T)[1] * 0.7, y = range(promATAC$WholeCortex_ASD_logFC, na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2_prom),2), "\n", ifelse(lm_p_prom < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p_prom, format = "e", digits = 1)))), col = "red")
  
  print(p_cur)
}
dev.off()
# Observation:
# cor = 0.1 (p<2e-7) for BATF/FOSL1/BACH1/NFE2-bound promoter ATAC-GE pairs. 

## What are the DEGs with diff. promoter ATAC bound by BATF/FOSL1/BACH1/NFE2
library(readxl)
Jill_DEG = readxl::read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")
promATAC$external_gene_name = Jill_DEG$external_gene_name[match(promATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
promATAC = promATAC[,c(1:2, ncol(promATAC), 3:(ncol(promATAC) - 1))]
sig_CorATAC_uniq$external_gene_name = Jill_DEG$external_gene_name[match(sig_CorATAC_uniq$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
sig_CorATAC_uniq = sig_CorATAC_uniq[,c(1:2, ncol(sig_CorATAC_uniq), 3:(ncol(sig_CorATAC_uniq) - 1))]

save(promATAC, sig_CorATAC_uniq, sig_CorATAC, tf_bound, tf_bound2, file = "DiffATAC_DEG_logFC_PromoterAndDistalATAC_eaTFbound.rda")

lnames = load("DiffATAC_DEG_logFC_PromoterAndDistalATAC_eaTFbound.rda")
Jill_DEG = readxl::read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

DEGup_promboundby_BatfFosl1Bach1 = unique(promATAC$external_gene_name[which((promATAC$BATF == "yes" | promATAC$FOSL1 == "yes" | promATAC$BACH1 == "yes") & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05 & promATAC$WholeCortex_ASD_logFC > 0 & promATAC$ATAC_logFC > 0)]) # 33 DEG up with promoter ATAC bound by BATF/FOSL1/BACH1

DEGup_promboundby_NFE2 = unique(promATAC$external_gene_name[which(promATAC$NFE2 == "yes" & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05 & promATAC$WholeCortex_ASD_logFC > 0 & promATAC$ATAC_logFC > 0)]) # 19 DEG up with promoter ATAC bound by NFE2

DEGup_promboundby_BatfFosl1Bach1Nfe2 = unique(c(DEGup_promboundby_BatfFosl1Bach1, DEGup_promboundby_NFE2)) # 34 DEGs, including NFIL3!

# GO enrichment on level 4-5 terms
BiocManager::install("clusterProfiler")
library(clusterProfiler) # T Wu etc. The Innovation. 2021
#?groupGO

# GO enrichment
# install.packages("gprofiler2")
#library(gprofiler2)
DEGup_promboundby_BatfFosl1Bach1Nfe2_geneid = Jill_DEG$ensembl_gene_id[match(DEGup_promboundby_BatfFosl1Bach1Nfe2, Jill_DEG$external_gene_name)]

#gost_DEGup_promboundby_BatfFosl1Bach1Nfe2 = gost(query = DEGup_promboundby_BatfFosl1Bach1Nfe2_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
#GO_DEGup_promboundby_BatfFosl1Bach1Nfe2 = gost_DEGup_promboundby_BatfFosl1Bach1Nfe2$result # apoptotic process
#GO_DEGup_promboundby_BatfFosl1Bach1Nfe2$minusLog10Pval = -log10(GO_DEGup_promboundby_BatfFosl1Bach1Nfe2$p_value)

gost_DEGup_promboundby_BatfFosl1Bach1Nfe2 = enrichGO(gene = DEGup_promboundby_BatfFosl1Bach1Nfe2_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
# gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_q005 = gost_DEGup_promboundby_BatfFosl1Bach1Nfe2@result
# gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_q005 = gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_q005[gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_q005$qvalue < 0.05,] # NULL
gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_filtered = gofilter(gost_DEGup_promboundby_BatfFosl1Bach1Nfe2, level = 4)
GO_DEGup_promboundby_BatfFosl1Bach1Nfe2 = gost_DEGup_promboundby_BatfFosl1Bach1Nfe2_filtered@result
GO_DEGup_promboundby_BatfFosl1Bach1Nfe2$minusLog10Pval = -log10(GO_DEGup_promboundby_BatfFosl1Bach1Nfe2$pvalue)

GO_df = GO_DEGup_promboundby_BatfFosl1Bach1Nfe2[1:10,]
GO_df$term_name = factor(GO_df$Description, levels = rev(GO_df$Description))

pdf("DEGpromoterATACbothUp_promoterATACboundbyBatfFosl1Bach1Nfe2_GOenrichmentL4.pdf", width = 8, height = 5)
GO_df[1:10,] %>%
  ggplot() +
  #geom_point(aes(x = term_name, y = minusLog10Pval), alpha = 0.8) +
  geom_bar(aes(x = term_name, y = minusLog10Pval), stat="identity", fill = "dodgerblue") +
  geom_hline(aes(yintercept = -log10(0.05), linetype = "-log10(0.05)")) +
  theme_bw() +
  coord_flip() + 
  scale_linetype_manual(name = "Theshold", values = 2) +
  ylab("-log10(p-value)") +
  xlab("GO terms") +
  ggtitle("GO enrichment of ASD up-reg genes with promoter ATAC up-reg\nand bound by BATF/FOSL1/BACH1/NFE2")
dev.off()
# gprofiler2 Observed: apoptotic process, cell death, regulation of signaling
# clusterProfiler Observed: necrotic cell death, cellular senescence, protein folding in ER, other terms are not legitimate

# Con't: DEGs with distal ATAC bound by CTCF
DEGdown_distalboundby_CTCF = unique(sig_CorATAC_uniq$external_gene_name[which((sig_CorATAC_uniq$CTCF_MA0139.1 == "yes" | sig_CorATAC_uniq$CTCF_MA1929.1 == "yes" | sig_CorATAC_uniq$CTCF_MA1930.1 == "yes") & sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_logFC < 0 & sig_CorATAC_uniq$distalATAC_logFC < 0)]) # 108 DEG down with distal ATAC bound by CTCF

DEGdown_promUPboundby_CTCF = unique(promATAC$external_gene_name[which((promATAC$CTCF_MA0139.1 == "yes" | promATAC$CTCF_MA1929.1 == "yes" | promATAC$CTCF_MA1930.1 == "yes") & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05 & promATAC$WholeCortex_ASD_logFC < 0 & promATAC$ATAC_logFC > 0)]) # 108 DEG down with distal ATAC bound by CTCF

DEGdown_promDOWNboundby_CTCF = unique(promATAC$external_gene_name[which((promATAC$CTCF_MA0139.1 == "yes" | promATAC$CTCF_MA1929.1 == "yes" | promATAC$CTCF_MA1930.1 == "yes") & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05 & promATAC$WholeCortex_ASD_logFC < 0 & promATAC$ATAC_logFC < 0)]) # 108 DEG down with distal ATAC bound by CTCF

# GO enrichment?
# install.packages("gprofiler2")
# library(gprofiler2)
DEGdown_distalboundby_CTCF_geneid = Jill_DEG$ensembl_gene_id[match(DEGdown_distalboundby_CTCF, Jill_DEG$external_gene_name)] # 108

# gost_DEGdown_distalboundby_CTCF = gost(query = DEGdown_distalboundby_CTCF_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
# GO_DEGdown_distalboundby_CTCF = gost_DEGdown_distalboundby_CTCF$result # nervous system development, synapse organization
# GO_DEGdown_distalboundby_CTCF$minusLog10Pval = -log10(GO_DEGdown_distalboundby_CTCF$p_value)
gost_DEGdown_distalboundby_CTCF = enrichGO(gene = DEGdown_distalboundby_CTCF_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
gost_DEGdown_distalboundby_CTCF_filtered = gofilter(gost_DEGdown_distalboundby_CTCF, level = 4)
GO_DEGdown_distalboundby_CTCF = gost_DEGdown_distalboundby_CTCF_filtered@result
GO_DEGdown_distalboundby_CTCF$minusLog10Pval = -log10(GO_DEGdown_distalboundby_CTCF$pvalue)

#GO_df = GO_DEGdown_distalboundby_CTCF[GO_DEGdown_distalboundby_CTCF$term_size < 6000,]
#GO_df = GO_df[1:10,]
#GO_df$term_name = factor(GO_df$term_name, levels = rev(GO_df$term_name))
GO_df = GO_DEGdown_distalboundby_CTCF[1:10,]
GO_df$term_name = factor(GO_df$Description, levels = rev(GO_df$Description))

save(GO_df, file = "DEGdistalATACbothDown_distalATACboundbyCTCF_GOenrichmentL4.rda")

pdf("DEGdistalATACbothDown_distalATACboundbyCTCF_GOenrichmentL4.pdf", width = 8, height = 5)
GO_df[1:10,] %>%
  ggplot() +
  #geom_point(aes(x = term_name, y = minusLog10Pval), alpha = 0.8) +
  geom_bar(aes(x = term_name, y = minusLog10Pval), stat="identity", fill = "dodgerblue") +
  geom_hline(aes(yintercept = -log10(0.05), linetype = "-log10(0.05)")) +
  theme_bw() +
  coord_flip() + 
  scale_linetype_manual(name = "Theshold", values = 2) +
  ylab("-log10(p-value)") +
  xlab("GO terms") +
  ggtitle("GO enrichment of ASD down-reg genes\nwith distal ATAC down-reg and bound by CTCF")
dev.off()
# Expected: synaptic terms

# Any SFARI genes?
SFARI = read.csv("~/Documents/Documents/Geschwind_lab/LAB/Database_download/SFARI_human_genes/SFARI-Gene_genes_01-11-2022release_03-03-2022export.csv")
SFARI_genes = SFARI$gene.symbol[SFARI$gene.score %in% c("S", "1", "2")] # 554 genes

intersect(DEGup_promboundby_BatfFosl1Bach1, SFARI_genes) # 0
intersect(DEGup_promboundby_NFE2, SFARI_genes) # 0
intersect(DEGdown_distalboundby_CTCF, SFARI_genes) # 6 genes: "SLC6A1" "ELAVL3" "TSC2"   "ANK3"   "NACC1"  "RAI1". 
intersect(DEGdown_promUPboundby_CTCF, SFARI_genes) # 6 genes: "PRICKLE1" "GPHN"     "KIF5C"    "GABRB2"   "ICA1"     "CAMK2B"
intersect(DEGdown_promDOWNboundby_CTCF, SFARI_genes) # 3 genes: "MAP1A"  "TSC2"   "ELAVL3"
# SFARI genes TSC2 and ELAVL3 have both promoter and distal ATAC down and gene expression down.

promATAC$SFARI = ifelse(promATAC$external_gene_name %in% SFARI_genes, "yes", "no")
sig_CorATAC_uniq$SFARI = ifelse(sig_CorATAC_uniq$external_gene_name %in% SFARI_genes, "yes", "no")
promATAC$SFARI = factor(promATAC$SFARI, levels = c("yes", "no"))
sig_CorATAC_uniq$SFARI = factor(sig_CorATAC_uniq$SFARI, levels = c("yes", "no"))

pdf("DiffATAC_DEG_logFC_PromoterATAC_eaTFbound_colSFARI.pdf", height = 4, width = 6)
idx_all = which(promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05)
for (i in 9:15) {
  tf = colnames(promATAC)[i]
  idx_plot = which(promATAC[,i] == "yes" & promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05)
  
  lm_res_prom = summary(lm(WholeCortex_ASD_logFC ~ 0 + ATAC_logFC, data = promATAC[idx_plot, ]))
  lm_coef_prom = lm_res_prom$coefficients[1]
  lm_p_prom = lm_res_prom$coefficients[4]
  lm_r2_prom = lm_res_prom$r.squared
  
  p_cur = promATAC[idx_plot, ] %>%
    ggplot(aes(x = ATAC_logFC, y = WholeCortex_ASD_logFC)) +
    geom_point(aes(col = SFARI)) + # col = "grey"
    geom_smooth(method='lm', formula= y~x) +
    #geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
    theme_bw() + 
    ggtitle(paste0("Promoter ATAC bound by ", tf)) +
    xlim(range(promATAC$ATAC_logFC[idx_all])[1], range(promATAC$ATAC_logFC[idx_all])[2]) +
    ylim(range(promATAC$WholeCortex_ASD_logFC[idx_all])[1], range(promATAC$WholeCortex_ASD_logFC[idx_all])[2]) +
    annotate("text", x = 0, y = range(promATAC$WholeCortex_ASD_logFC[idx_plot], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2_prom),2), "\n", ifelse(lm_p_prom < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p_prom, format = "e", digits = 1)))), col = "black")
  
  print(p_cur)
}
dev.off()

pdf("DiffATAC_DEG_logFC_DistalATAC_eaTFbound_colSFARI.pdf", height = 4, width = 6)
idx_all = which(sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05)
for (i in 8:14) {
  tf = colnames(sig_CorATAC_uniq)[i]
  idx_plot = which(sig_CorATAC_uniq[,i] == "yes" & sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05)
  
  lm_res = summary(lm(WholeCortex_ASD_logFC ~ 0 + distalATAC_logFC, data = sig_CorATAC_uniq[idx_plot, ]))
  lm_coef = lm_res$coefficients[1]
  lm_p = lm_res$coefficients[4]
  lm_r2 = lm_res$r.squared
  
  p_cur = sig_CorATAC_uniq[idx_plot, ] %>%
    ggplot(aes(x = distalATAC_logFC, y = WholeCortex_ASD_logFC)) +
    geom_point(aes(col = SFARI)) + # col = "grey"
    geom_smooth(method='lm', formula= y~x) +
    #geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
    theme_bw() + 
    ggtitle(paste0("Distal ATAC bound by ", tf)) +
    xlim(range(sig_CorATAC_uniq$distalATAC_logFC[idx_all])[1], range(sig_CorATAC_uniq$distalATAC_logFC[idx_all])[2]) +
    ylim(range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_all])[1], range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_all])[2]) +
    annotate("text", x = 0, y = range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_plot], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2),2), "\n", ifelse(lm_p < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p, format = "e", digits = 1)))), col = "black")
  
  print(p_cur)
}
dev.off()
# Observations:
# For BATF/FOSL1/BACH1/NFE2, 0-1 SFARI gene has distal ATAC and DEG down (Q1), 1-3 SFARI genes have distal ATAC and DEG down (Q3), and 1 SFARI gene in Q2 and Q4 each.
# For CTCF, 4-6 SFARI genes have distal ATAC and DEG down (Q3), 2-3 SFARI genes have distal ATAC and DEG up (Q1), 3 SFARI genes have distal ATAC up while DEG down (Q4), 0-2 SFARI gene have distal ATAC down while DEG up (Q2).

## For CTCF-bound differential distal ATAC-DEG pairs, how do Hi-C loop change? 
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_10_OverlapATACwTFbound/")
lnames = load("DiffATAC_DEG_logFC_PromoterAndDistalATAC_eaTFbound.rda")
library(readxl)
Jill_DEG = readxl::read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

sig_CorATAC$CTCF_MA0139.1 = ifelse(sig_CorATAC$distal_ATAC %in% tf_bound2[["CTCF-MA0139.1"]], "yes", "no")
sig_CorATAC$CTCF_MA1929.1 = ifelse(sig_CorATAC$distal_ATAC %in% tf_bound2[["CTCF-MA1929.1"]], "yes", "no")
sig_CorATAC$CTCF_MA1930.1 = ifelse(sig_CorATAC$distal_ATAC %in% tf_bound2[["CTCF-MA1930.1"]], "yes", "no")

sig_CorATAC$external_gene_name = Jill_DEG$external_gene_name[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
sig_CorATAC = sig_CorATAC[,c(1:2, ncol(sig_CorATAC), 3:(ncol(sig_CorATAC) - 1))]

DEG_distalATAC_bothDown_boundby_CTCF = sig_CorATAC[sig_CorATAC$distalATAC_logFC < 0 & sig_CorATAC$WholeCortex_ASD_logFC < 0 & sig_CorATAC$distalATAC_FDR < 0.05 & sig_CorATAC$WholeCortex_ASD_FDR < 0.05 & (sig_CorATAC$CTCF_MA0139.1 == "yes" | sig_CorATAC$CTCF_MA1929.1 == "yes" | sig_CorATAC$CTCF_MA1930.1 == "yes"), ] # 176 rows

# Add Hi-C loop logFC
lnames = load("../../../../../Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/diffloop_limma_results/Limma_ASD_Bulk_lmDiagnosisBatchValidTransAgeScaledSex.rda") # top.table_ASD

DEG_distalATAC_bothDown_boundby_CTCF$Loop_Bulk_logFC = top.table_ASD$logFC[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_ASD))]
DEG_distalATAC_bothDown_boundby_CTCF$Loop_Bulk_P = top.table_ASD$P.Value[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_ASD))]

lnames = load("../../../../../Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/diffloop_limma_results/Limma_Dup15q_NeuNp_lmDiagnosisTransAgeScaledValid.rda") # top.table_Dup15q_NeuNp

DEG_distalATAC_bothDown_boundby_CTCF$Loop_NeuNp_logFC = top.table_Dup15q_NeuNp$logFC[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_Dup15q_NeuNp))]
DEG_distalATAC_bothDown_boundby_CTCF$Loop_NeuNp_P = top.table_Dup15q_NeuNp$P.Value[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_Dup15q_NeuNp))]

lnames = load("../../../../../Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/diffloop_limma_results/Limma_Dup15q_NeuNn_lmDiagnosisTransDuplicateValid.rda") # table_Dup15q_NeuNn

DEG_distalATAC_bothDown_boundby_CTCF$Loop_NeuNn_logFC = top.table_Dup15q_NeuNn$logFC[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_Dup15q_NeuNn))]
DEG_distalATAC_bothDown_boundby_CTCF$Loop_NeuNn_P = top.table_Dup15q_NeuNn$P.Value[match(DEG_distalATAC_bothDown_boundby_CTCF$Loop, rownames(top.table_Dup15q_NeuNn))]

# Bargraph showing percent Hi-C up/down separated by cell-type
# a. Nominal Hi-C logFC change?
df = DEG_distalATAC_bothDown_boundby_CTCF
DEG_distalATAC_bothDown_boundby_CTCF_summary = 
  data_frame(Celltype = rep(c("Bulk", "NeuNp", "NeuNn"), each = 3), HiC_dir = rep(c("NS", "Up", "Down"), 3), N_Loop = 0)
DEG_distalATAC_bothDown_boundby_CTCF_summary
Loop_PT = 0.1 # short for Loop_p_threshold
DEG_distalATAC_bothDown_boundby_CTCF_summary$N_Loop = 
  c(length(which(df$Loop_Bulk_P > Loop_PT)), 
    length(which(df$Loop_Bulk_P < Loop_PT & df$Loop_Bulk_logFC > 0)),
    length(which(df$Loop_Bulk_P < Loop_PT & df$Loop_Bulk_logFC < 0)),
    length(which(df$Loop_NeuNp_P > Loop_PT)), 
    length(which(df$Loop_NeuNp_P < Loop_PT & df$Loop_NeuNp_logFC > 0)),
    length(which(df$Loop_NeuNp_P < Loop_PT & df$Loop_NeuNp_logFC < 0)),
    length(which(df$Loop_NeuNn_P > Loop_PT)), 
    length(which(df$Loop_NeuNn_P < Loop_PT & df$Loop_NeuNn_logFC > 0)),
    length(which(df$Loop_NeuNn_P < Loop_PT & df$Loop_NeuNn_logFC < 0))
    )
DEG_distalATAC_bothDown_boundby_CTCF_summary$HiC_dir = factor(DEG_distalATAC_bothDown_boundby_CTCF_summary$HiC_dir, levels = c("Up", "NS", "Down"))

AllHiC_summary = 
  data_frame(Celltype = rep(c("Bulk", "NeuNp", "NeuNn"), each = 3), HiC_dir = rep(c("NS", "Up", "Down"), 3), N_Loop = 0)
AllHiC_summary$N_Loop = 
  c(length(which(top.table_ASD$P.Value > Loop_PT)), 
    length(which(top.table_ASD$P.Value < Loop_PT & top.table_ASD$logFC > 0)),
    length(which(top.table_ASD$P.Value < Loop_PT & top.table_ASD$logFC < 0)),
    length(which(top.table_Dup15q_NeuNp$P.Value > Loop_PT)), 
    length(which(top.table_Dup15q_NeuNp$P.Value < Loop_PT & top.table_Dup15q_NeuNp$logFC > 0)),
    length(which(top.table_Dup15q_NeuNp$P.Value < Loop_PT & top.table_Dup15q_NeuNp$logFC < 0)),
    length(which(top.table_Dup15q_NeuNn$P.Value > Loop_PT)), 
    length(which(top.table_Dup15q_NeuNn$P.Value < Loop_PT & top.table_Dup15q_NeuNn$logFC > 0)),
    length(which(top.table_Dup15q_NeuNn$P.Value < Loop_PT & top.table_Dup15q_NeuNn$logFC < 0))
  )
AllHiC_summary$HiC_dir = factor(AllHiC_summary$HiC_dir, levels = c("Up", "NS", "Down"))

# plot DEGdistalATACbothDown_distalATACboundbyCTCF_HiCloopSignifChanges.pdf after c. statistical test

# b. Any Hi-C logFC change?
DEG_distalATAC_bothDown_boundby_CTCF_summary2 = 
  data_frame(Celltype = rep(c("Bulk", "NeuNp", "NeuNn"), each = 2), HiC_dir = rep(c("Up", "Down"), 3), N_Loop = 0)
DEG_distalATAC_bothDown_boundby_CTCF_summary2
DEG_distalATAC_bothDown_boundby_CTCF_summary2$N_Loop = 
  c(length(which(df$Loop_Bulk_logFC > 0)),
    length(which(df$Loop_Bulk_logFC < 0)),
    length(which(df$Loop_NeuNp_logFC > 0)),
    length(which(df$Loop_NeuNp_logFC < 0)),
    length(which(df$Loop_NeuNn_logFC > 0)),
    length(which(df$Loop_NeuNn_logFC < 0))
  )
DEG_distalATAC_bothDown_boundby_CTCF_summary2$HiC_dir = factor(DEG_distalATAC_bothDown_boundby_CTCF_summary2$HiC_dir, levels = c("Up", "Down"))

AllHiC_summary2 = 
  data_frame(Celltype = rep(c("Bulk", "NeuNp", "NeuNn"), each = 2), HiC_dir = rep(c("Up", "Down"), 3), N_Loop = 0)
AllHiC_summary2$N_Loop = 
  c(length(which(top.table_ASD$logFC > 0)),
    length(which(top.table_ASD$logFC < 0)),
    length(which(top.table_Dup15q_NeuNp$logFC > 0)),
    length(which(top.table_Dup15q_NeuNp$logFC < 0)),
    length(which(top.table_Dup15q_NeuNn$logFC > 0)),
    length(which(top.table_Dup15q_NeuNn$logFC < 0))
  )
AllHiC_summary2$HiC_dir = factor(AllHiC_summary2$HiC_dir, levels = c("Up", "Down"))

# plot after statistical test

# c. Statistical test - bootstrap
ls(pattern = ".*_summary.*")
df_p = data_frame(Celltype = rep(c("Bulk", "NeuNp", "NeuNn"), each = 2), HiCchangeType = rep(c("Any", "Sig"), 3), HiCdown_p = 1)

for (i in 1:nrow(df_p)) {
  ct = df_p$Celltype[i]
  if (ct == "Bulk") {
    hic_actual = top.table_ASD
  } else if (ct == "NeuNp") {
    hic_actual = top.table_Dup15q_NeuNp
  } else {
    hic_actual = top.table_Dup15q_NeuNn
  }
  
  hic_change = df_p$HiCchangeType[i]
  if (hic_change == "Sig") {
    df_actual = DEG_distalATAC_bothDown_boundby_CTCF_summary
    df_bg = AllHiC_summary
  } else {
    df_actual = DEG_distalATAC_bothDown_boundby_CTCF_summary2
    df_bg = AllHiC_summary2
  }
  HiCdown_actual = df_actual$N_Loop[df_actual$Celltype == ct & df_actual$HiC_dir == "Down"]
  HiC_total = sum(df_actual$N_Loop[df_actual$Celltype == ct])
  
  set.seed(284)
  HiCdown_bg_c = c()
  for (B in 1:100) { # B up to 500 does not change the resulting p-value
    idx_random = sample(1:nrow(hic_actual), HiC_total)
    if (hic_change == "Sig") {
      HiCdown_bg = length(which(hic_actual$logFC[idx_random] < 0 & hic_actual$P.Value[idx_random] < Loop_PT))
    } else {
      HiCdown_bg = length(which(hic_actual$logFC[idx_random] < 0))
    }
    HiCdown_bg_c = c(HiCdown_bg_c, HiCdown_bg)
  }
  df_p$HiCdown_p[i] = length(which(HiCdown_bg_c > HiCdown_actual))/100
}
df_p
# Bulk and NeuNn Hi-C loop of CTCF-bound diff. distal ATAC-DEG pairs are significantly more likely to be down-regulated (hic p-value < 0.1) in ASD.
# NeuNp Hi-C loop of CTCF-bound diff. distal ATAC-DEG pairs are significantly more likely to be have reduced intensity (no hic p-val threshold) in ASD.
# Add p-values back to the plots

save(DEG_distalATAC_bothDown_boundby_CTCF_summary, AllHiC_summary, DEG_distalATAC_bothDown_boundby_CTCF_summary2, AllHiC_summary2, df_p, file = "DEGdistalATACbothDown_distalATACboundbyCTCF_HiCloopSignifOrAnyChanges.rda")

# plot a.
pdf("DEGdistalATACbothDown_distalATACboundbyCTCF_HiCloopSignifChanges.pdf", width = 6, height = 5)
DEG_distalATAC_bothDown_boundby_CTCF_summary %>%
  ggplot(aes(x = Celltype, y = N_Loop, fill = HiC_dir)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=N_Loop), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5) +
  xlab("Cell-type of Hi-C interactions") +
  ylab("Number of promoter-distal ATAC interactions") +
  scale_fill_discrete(name = "Hi-C change (ASD vs. CTL)\nsignificance threshold: p<0.1") +
  theme_bw() +
  ggtitle("CTCF-bound distal ATAC-GE pairs\nwith down-regulated ATAC and gene expression") +
  theme(legend.position="bottom") +
  annotate("text", x = 1.33, y = 20, label = paste0("p = ", df_p$HiCdown_p[df_p$Celltype == "Bulk" & df_p$HiCchangeType == "Sig"])) +
  annotate("text", x = 2.3, y = 20, label = paste0("p = ", df_p$HiCdown_p[df_p$Celltype == "NeuNn" & df_p$HiCchangeType == "Sig"]))

AllHiC_summary %>%
  ggplot(aes(x = Celltype, y = N_Loop, fill = HiC_dir)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=N_Loop), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5) +
  xlab("Cell-type of Hi-C interactions") +
  ylab("Number of Hi-C interactions") +
  scale_fill_discrete(name = "Hi-C change (ASD vs. CTL)\nsignificance threshold: p<0.1") +
  theme_bw() +
  ggtitle(" \nAll Hi-C interactions") +
  theme(legend.position="bottom")

dev.off()
# Observation: Mostly no nominal Hi-C loop changes. More norminally down-reg loops than up-reg loops in Bulk and NeuNn (loop p < 0.1).

# plot b.
pdf("DEGdistalATACbothDown_distalATACboundbyCTCF_HiCloopAnyChanges.pdf", width = 6, height = 5)
DEG_distalATAC_bothDown_boundby_CTCF_summary2 %>%
  ggplot(aes(x = Celltype, y = N_Loop, fill = HiC_dir)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=N_Loop), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5) +
  xlab("Cell-type of Hi-C interactions") +
  ylab("Number of promoter-distal ATAC interactions") +
  scale_fill_discrete(name = "Hi-C change (ASD vs. CTL)\nno p-value threshold") +
  theme_bw() +
  ggtitle("CTCF-bound distal ATAC-GE pairs\nwith down-regulated ATAC and gene expression") +
  theme(legend.position="bottom") +
  annotate("text", x = 3.23, y = 45, label = paste0("p = ", df_p$HiCdown_p[df_p$Celltype == "NeuNp" & df_p$HiCchangeType == "Any"]))

AllHiC_summary2 %>%
  ggplot(aes(x = Celltype, y = N_Loop, fill = HiC_dir)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=N_Loop), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5) +
  xlab("Cell-type of Hi-C interactions") +
  ylab("Number of Hi-C interactions") +
  scale_fill_discrete(name = "Hi-C change (ASD vs. CTL)\nno p-value threshold") +
  theme_bw() +
  ggtitle(" \nAll Hi-C interactions") +
  theme(legend.position="bottom")

dev.off()
# Observation: 
# Background Hi-C loops in Bulk, NeuNp, NeuNn have similar number of up and down-reg Hi-C interactions
# More down-reg distalATAC-DEG pairs have down-reg (logFC < 0, no p threshold) Hi-C interactions in Dup15q NeuNp and NeuNn than up-reg Hi-C interactions. No difference in Bulk. 

# Any SFARI genes have down-reg Hi-C?
DEG_distalATAC_bothDown_boundby_CTCF$SFARI = ifelse(DEG_distalATAC_bothDown_boundby_CTCF$external_gene_name %in% SFARI_genes, "yes", "no")

tmp = DEG_distalATAC_bothDown_boundby_CTCF[which(DEG_distalATAC_bothDown_boundby_CTCF$SFARI == "yes"),] # 10 rows, ANK3 promoter-distal ATAC interaction chr10_62115000_62145000 nominally down-reg

save(list = ls(), file = "save_all.rda")

# ---- after meeting with Dan on 12/07/2022 ------
## Dan suggests apoptotic process enrichment in up-reg TF bound DEGs is wierd. Look at them one by one.
lnames = load("save_all.rda")

DEGup_promboundby_BatfFosl1Bach1Nfe2 
# [1] "CACHD1"     "BCAR3"      "PEA15"      "ALDH18A1"   "ADD3"       "ST5"       
# [7] "MDK"        "UBE2L6"     "YBX3"       "VDR"        "PHF11"      "HEATR5A"   
# [13] "SAMD4A"     "STON2"      "PDIA3"      "SPPL2A"     "NFATC1"     "CSRNP1"    
# [19] "BOC"        "ACTRT3"     "LRRC34"     "SGMS2"      "USP53"      "OTUD4"     
# [25] "ZFR"        "MAP3K5"     "HEBP2"      "COG5"       "CLU"        "GADD45G"   
# [31] "NFIL3"      "NEK6"       "TIMP1"      "AC084082.3"
# Note that NFIL3 is up. 

# Any in Parikshak's modules M9/19/20?
library(readxl)
Parikshak_GM = read_excel("~/Documents/Documents/Geschwind_lab/LAB/Literature/ASD_omics_wide_changes/Parikshak_2016_SupTable2_GeneModuleMembership.xlsx", skip = 1)
#M9_genes = Parikshak_GM$`HGNC Symbol`[Parikshak_GM$`WGCNA Module Label` == 9]
#M9_genes[duplicated(M9_genes)] # NA no HGNC symbol
M9_genes = unique(Parikshak_GM$`HGNC Symbol`[Parikshak_GM$`WGCNA Module Label` == 9]) # 507 genes, enriched in astrocytes
M19_genes = unique(Parikshak_GM$`HGNC Symbol`[Parikshak_GM$`WGCNA Module Label` == 19]) # 274 genes, enriched in microglia
M20_genes = unique(Parikshak_GM$`HGNC Symbol`[Parikshak_GM$`WGCNA Module Label` == 20]) # 333 genes, enriched in neurons/astrocytes

intersect(DEGup_promboundby_BatfFosl1Bach1Nfe2, M9_genes) # 11 genes, wow! - test significant overlap
intersect(DEGup_promboundby_BatfFosl1Bach1Nfe2, M19_genes) # 1 gene: YBX3
intersect(DEGup_promboundby_BatfFosl1Bach1Nfe2, M20_genes) # 2 genes: CSRNP1, ZFR

# test significant overlap b/w DEGup_promboundby_BatfFosl1Bach1Nfe2 and M9_genes
all_genes = intersect(Jill_DEG$external_gene_name, Parikshak_GM$`HGNC Symbol`) # 14196 genes
M9_genes = intersect(M9_genes, all_genes) # 504 genes in Jill's dataset

both = length(intersect(DEGup_promboundby_BatfFosl1Bach1Nfe2, M9_genes)) # 11 genes
M9only = length(setdiff(M9_genes, DEGup_promboundby_BatfFosl1Bach1Nfe2)) # 493 genes
TFonly = length(setdiff(DEGup_promboundby_BatfFosl1Bach1Nfe2, M9_genes)) # 23 genes
neither = length(setdiff(all_genes, union(M9_genes, DEGup_promboundby_BatfFosl1Bach1Nfe2))) # 13672
fisher_res = fisher.test(matrix(c(both, TFonly, M9only, neither), nrow = 2))
fisher_p = formatC(fisher_res$p.value, digits = 2) # 1.38e-08
#fisher_OD = (both * neither)/(M9only * TFonly) # 13.3
fisher_OD = fisher_res$estimate # 13.3

# Venn diagram
library(VennDiagram)
library(gridExtra)

save(DEGup_promboundby_BatfFosl1Bach1Nfe2, M9_genes, both, fisher_res, fisher_p, fisher_OD, file = "VennDiagram_Overlap_DEGupPromboundbyBatfFosl1Bach1Nfe2_ParikshakM9.rda")

png("VennDiagram_Overlap_DEGupPromboundbyBatfFosl1Bach1Nfe2_ParikshakM9.png", height = 500, width = 500)
grid.newpage()
my_venndiag <- draw.pairwise.venn(length(DEGup_promboundby_BatfFosl1Bach1Nfe2), length(M9_genes), both, fill = c("orange","chartreuse"), alpha = 0.6, category = c("Up-regulated \ngenes with \n promoter ATAC \n bound by BATF,\nFOSL1, BACH1 \nor NFE2", "Parikshak 2016 M9 module genes"), cat.pos = c(0,0), cat.dist = c(0.08, 0.02), cex = 2, cat.cex = 1.5, fontfamily = "Arial", cat.fontfamily = "Arial")
grid.arrange(gTree(children = my_venndiag), # Add title & subtitle
             #top = "My Main Title",
             bottom = textGrob(paste0("Significant overlap, p = ", fisher_p), gp = gpar(fontsize=20)), heights = c(30,1))
dev.off()
dev.off()

## How do JDP2/NFIL3/JUN/JUNB/FOX/BATF3/TEF/BACH1 bound ATAC target genes overlap with scRNA regulons? 
# Run script 10_3_11_overlap_TFboundATACtargetGenes_with_LucyRegulonScRNA.R


