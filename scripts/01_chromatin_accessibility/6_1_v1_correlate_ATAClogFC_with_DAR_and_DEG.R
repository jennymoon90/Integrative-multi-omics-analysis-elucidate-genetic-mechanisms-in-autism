rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/")

## Load differential ATAC results
load("../5_DiffATAC/DiffATAC.rda")
DiffATAC_location = DiffATAC[,1:3]

## Load differential H3K27Ac results
lnames = load("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/32_Gokul_H3K27Ac_datasets/processed_data/differential_acetylation/Convergent.Rdata") # DAR

DAR_location = as.data.frame(str_split_fixed(rownames(DAR), "-", 3))
colnames(DAR_location) = c("chr", "start", "end")
DAR_location$start = as.integer(DAR_location$start)
DAR_location$end = as.integer(DAR_location$end)
DAR = cbind(DAR_location, DAR)
rownames(DAR_location) = rownames(DAR)

unique(DiffATAC$chr) # chr1-Y
unique(DAR$chr) # chr1-Y, chrM and xxx_random
DAR = DAR[-which(grepl("_random", DAR$chr)),] # 56451 H3K27ac peaks
DAR = DAR[DAR$chr != "chrM",] # 56449 H3K27ac peaks

## Load DEG results
library(readxl)
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

## Peak size histograms
DiffATAC_location$distance = DiffATAC_location$end - DiffATAC_location$start
DAR_location = DAR_location[rownames(DAR_location) %in% rownames(DAR),]
DAR_location$distance = DAR_location$end - DAR_location$start

meanATAClen = round(mean(DiffATAC_location$distance)) # 1184 bp
meanH3K27AClen = round(mean(DAR_location$distance)) # 2949 bp

pdf("ATAC_and_H3K27Ac_peaks_size_distribution.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
hist(DiffATAC_location$distance, breaks = 50, xlab = "peak length", 
     main = paste0(nrow(DiffATAC), " ATAC peaks\nmean length: ",meanATAClen," bp"))
hist(DAR_location$distance, breaks = 100, xlab = "peak length",
     main = paste0(nrow(DAR), " H3K27ac peaks\nmean length: ",meanH3K27AClen," bp"))
dev.off()

## Correlation of logFC of differential ATAC and H3K27ac peaks
# Filter for FDR < 0.05 of all profiles
DiffATAC = DiffATAC[!is.na(DiffATAC$ATAC_FDR),] # 127971 ATAC peaks
DAR = DAR[!is.na(DAR$beta.condition),] # 56168 H3K27ac peaks
sigATAC = DiffATAC[DiffATAC$ATAC_FDR < 0.05,] # 5033 significant ATAC peaks
sigH3K27ac = DAR[DAR$p.condition.fdr < 0.05,] # 8158 significant H3K27ac peaks

# Intersect the two significant peak sets
sigATAC_gr = annoDF2GR(sigATAC)
sigH3K27ac_gr = annoDF2GR(sigH3K27ac)
hits = as.data.frame(findOverlaps(sigATAC_gr, sigH3K27ac_gr))
hits$ATACpeak = rownames(sigATAC)[hits$queryHits]
hits$H3K27ACpeak = rownames(sigH3K27ac)[hits$subjectHits]
hits$ATAC_logFC = sigATAC$ATAC_logFC[hits$queryHits]
hits$H3K27ac_logFC = sigH3K27ac$beta.condition[hits$subjectHits]

# Venn diagram
length(unique(hits$ATACpeak)) # 818
length(unique(hits$H3K27ACpeak)) # 801
# make venn diagram using ppt. statistics below:

# Out of all accessible regions, how many bp belong to differential ATAC peaks and H3K27ac peaks? Is it significant?
write.table(DiffATAC[,1:3], file = "venn_diagram/ATACpeaks.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(DAR[,1:3], file = "venn_diagram/H3K27ACpeaks.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(sigATAC[,1:3], file = "venn_diagram/sigATAC.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(sigH3K27ac[,1:3], file = "venn_diagram/sigH3K27ac.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# Come back after running 6_1_bedtools_overlap_ATAC_H3K27ac.sh

all_accessible = read.table("venn_diagram/all_accessible.merged.bed")
all_accessible$len = all_accessible$V3 - all_accessible$V2 + 1
total_bp = sum(all_accessible$len) # 242,513,186 bp

sig_intersect = read.table("venn_diagram/sig_intersect.bed")
sig_intersect$len = sig_intersect$V3 - sig_intersect$V2 + 1
intersect_bp = sum(sig_intersect$len) # 1,292,719 bp

InSigatacNotSigh3k27ac = sigATAC[! rownames(sigATAC) %in% hits$ATACpeak,] # 4215 peaks
InSigatacNotSigh3k27ac$len = InSigatacNotSigh3k27ac$end - InSigatacNotSigh3k27ac$start + 1
InSigatacNotSigh3k27ac_bp = sum(InSigatacNotSigh3k27ac$len) # 5,570,806 bp

InSigh3k27acNotSigatac = sigH3K27ac[! rownames(sigH3K27ac) %in% hits$H3K27ACpeak,] # 7357 peaks
InSigh3k27acNotSigatac$len = InSigh3k27acNotSigatac$end - InSigh3k27acNotSigatac$start + 1
InSigh3k27acNotSigatac_bp = sum(InSigh3k27acNotSigatac$len) # 22,546,475 bp

sig_union = read.table("venn_diagram/sig_union.merged.bed")
sig_union$len = sig_union$V3 - sig_union$V2 + 1
sig_union_bp = sum(sig_union$len)
out_bp = total_bp - sig_union_bp # 211,481,821

fisher_res = fisher.test(matrix(c(intersect_bp, InSigatacNotSigh3k27ac_bp, InSigh3k27acNotSigatac_bp, out_bp), nrow = 2))
fisher_res$p.value # 0, p-value < 2.2e-16
fisher_res$estimate # OR = 2.2

save(hits, sigATAC, sigH3K27ac, intersect_bp, fisher_res, file = "venn_diagram/VennDiagram_DAR_overlap.rda")

# failed using bedr
# ---------
# install.packages("bedr")
# library(bedr)
# check.binary(x = "bedtools", verbose = TRUE) # F

# ATAC = paste0(DiffATAC$chr, ":", DiffATAC$start, "-", DiffATAC$end)
# is.a.valie = is.valid.region(ATAC) # all PASS
# H3K27ac = paste0(DAR$chr, ":", DAR$start, "-", DAR$end)
# is.a.valie = is.valid.region(H3K27ac) # all PASS
# 
# is.sorted = is.sorted.region(ATAC) # TRUE
# is.sorted = is.sorted.region(H3K27ac) # FALSE
# H3K27ac.sorted = bedr.sort.region(H3K27ac)
# 
# # total bp
# All = c(ATAC, H3K27ac.sorted)
# All.sorted = bedr.sort.region(All)
# All.merge = bedr.merge.region(All.sorted) # ERROR: missing binary/executable bedtools * Collapsing
# ---------

# Plot correlation of ATAC and H3K27ac logFCs
hit_cor = round(cor(hits$ATAC_logFC, hits$H3K27ac_logFC),2) # 0.82
tmp = cor.test(hits$ATAC_logFC, hits$H3K27ac_logFC)
hit_corsig = formatC(tmp$p.value, 1) # 4e-197

df = data_frame(x = hits$ATAC_logFC, y = hits$H3K27ac_logFC)
model = lm(y ~ 0 + x, df)
model_res = summary(model)
lm_coef = round(model_res$coefficients[1,1], digits = 2) # 1.26
r_coef = round(sqrt(model_res$r.squared), digits = 2) # 0.89
p_coef = formatC(model_res$coefficients[1,4], digits = 1) # 1e-275

pdf("LogFC_of_ASD_significant_ATAC_and_H3K27Ac_peaks.pdf", height = 4, width = 6)
plot(hits$ATAC_logFC, hits$H3K27ac_logFC, xlab = "ATAC peak logFC", ylab = "H3K27ac peak logFC", pch = 19, col = alpha("black",0.5), main = "Correlation between logFC of differential ATAC and H3K27ac") 
#abline(a = 0, b = hit_cor, col = "red", lwd = 1.5)
#text(x = 0.4, y = -0.1, paste0("Pearson's correlation\ncoef = ", hit_cor, "\np = ", hit_corsig), col = "red")
abline(a = 0, b = lm_coef, col = "red", lwd = 1.5)
text(x = 0.4, y = 0, paste0("R^2 = ", round(r_coef^2, 1), "\np = ", p_coef), col = "red")
dev.off()
# Great, it correlates so well! Finally I learnt how to do differential analysis.

# save the variables for the logFC cor
save(hits, lm_coef, r_coef, p_coef, file = "LogFC_of_ASD_significant_ATAC_and_H3K27Ac_peaks_20221006.rda")

## Correlation of logFC of differential promoter ATAC and DEG
# get promoter regions (tss +- 2000)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
genes = genes(txdb)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("external_gene_name","ensembl_gene_id","entrezgene_id"), values = names(genes), filters = 'entrezgene_id', mart = mart)
save(bm, file = "BioMart_Entrezgene_id_2_Ensembl_gene_id_and_name.rda")

promoter = as.data.frame(promoters(genes, upstream = 2000, downstream=2000))
colnames(promoter)[which(colnames(promoter) == "gene_id")] = "entrezgene_id"
promoter$external_gene_name = bm$external_gene_name[match(promoter$entrezgene_id,bm$entrezgene_id)]
promoter$ensembl_gene_id = bm$ensembl_gene_id[match(promoter$entrezgene_id,bm$entrezgene_id)]
colnames(promoter)[1] = "chr"
save(promoter, file = "promoter2kb.rda")

# Link diff ATAC peaks to DEGs
sigDEG = Jill_DEG[Jill_DEG$WholeCortex_ASD_FDR < 0.05,] # 4219 genes
sigDEG_promoter = promoter[promoter$ensembl_gene_id %in% sigDEG$ensembl_gene_id,-4]
sigDEG_promoter$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(sigDEG_promoter$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
sigDEG_promoter_gr = annoDF2GR(sigDEG_promoter)

hits2 = as.data.frame(findOverlaps(sigATAC_gr, sigDEG_promoter_gr))
hits2$ATACpeak = rownames(sigATAC)[hits2$queryHits]
hits2$gene_id = sigDEG_promoter$ensembl_gene_id[hits2$subjectHits]
hits2$ATAC_logFC = sigATAC$ATAC_logFC[hits2$queryHits]
hits2$DEG_logFC = sigDEG_promoter$WholeCortex_ASD_logFC[hits2$subjectHits]

# Plot correlation of ATAC and H3K27ac logFCs
hit2_cor = round(cor(hits2$ATAC_logFC, hits2$DEG_logFC),2) # 0.2 -> 0.11 (2kb -> 3kb)
tmp = cor.test(hits2$ATAC_logFC, hits2$DEG_logFC)
hit2_corsig = formatC(tmp$p.value, 1) # 0.002

df2 = data_frame(x = hits2$ATAC_logFC, y = hits2$DEG_logFC)
model = lm(y ~ 0 + x, df2)
model_res = summary(model)
lm_coef = round(model_res$coefficients[1,1], digits = 2) # 0.21
r_coef = round(sqrt(model_res$r.squared), digits = 2) # 0.21
p_coef = formatC(model_res$coefficients[1,4], digits = 1) # 0.002 -> 0.04

pdf("LogFC_of_ASD_significant_2kbpromoterATAC_and_DEG.pdf", height = 4, width = 6)
plot(hits2$ATAC_logFC, hits2$DEG_logFC, xlab = "Promoter ATAC peak logFC (TSS +- 2kb)", ylab = "DEG  logFC", pch = 19, col = alpha("black",0.5), main = "Correlation between logFC of differential\nATAC and gene expression") 
abline(a = 0, b = lm_coef, col = "red", lwd = 1.5)
text(x = 0, y = 0.5, paste0("R^2 = ", round(r_coef^2, 2), "\np = ", p_coef), col = "red")
dev.off()
# Significant correlation, but many dots in Q4 - increased promoter ATAC but reduced GE.

## pSI enrichment for the 4 quadrants
hits2$gene_name = sigDEG_promoter$external_gene_name[hits2$subjectHits]
hits2$Quadrants = case_when(
  hits2$ATAC_logFC > 0 & hits2$DEG_logFC > 0 ~ "Q1",
  hits2$ATAC_logFC < 0 & hits2$DEG_logFC > 0 ~ "Q2",
  hits2$ATAC_logFC < 0 & hits2$DEG_logFC < 0 ~ "Q3",
  hits2$ATAC_logFC > 0 & hits2$DEG_logFC < 0 ~ "Q4"
)

# Load pSI data
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/Zhang_2015_pSI/Zhang2015pSI_usingMatureOfAllAges.rda") 
colnames(pSI_res)

# pSI enrichment function
pSI_enrichment = function(df) {
  idx_Astro = which(df$gene_name %in% rownames(pSI_res)[which(pSI_res$Astrocyte<0.05)])
  idx_Neuron = which(df$gene_name %in% rownames(pSI_res)[which(pSI_res$Neuron<0.05)]) 
  idx_Oligo = which(df$gene_name %in% rownames(pSI_res)[which(pSI_res$Oligodendrocyte<0.05)])
  idx_Microglia = which(df$gene_name %in% rownames(pSI_res)[which(pSI_res$Microglia<0.05)])
  idx_Endoth = which(df$gene_name %in% rownames(pSI_res)[which(pSI_res$Endothelial<0.05)])
  
  df$SI_Astro = 
    df$SI_Neuron = 
    df$SI_Oligo = 
    df$SI_Microglia = 
    df$SI_Endoth = 
    df$SI_Glia = "non-specific"
  
  df$SI_Astro[idx_Astro] = "specific"
  df$SI_Neuron[idx_Neuron] = "specific"
  df$SI_Oligo[idx_Oligo] = "specific"
  df$SI_Microglia[idx_Microglia] = "specific"
  df$SI_Endoth[idx_Endoth] = "specific"
  df$SI_Glia[unique(c(idx_Astro, idx_Oligo, idx_Microglia))] = "specific"
  
  p_Fisher = OR_Fisher = as.data.frame(matrix(nrow = 4, ncol = 6))
  Quadrants = c("Q1", "Q2","Q3", "Q4")
  colnames(df)
  k = which(colnames(df) == "SI_Glia") - 1
  
  for (i in 1:length(Quadrants)) { # Quadrants
    for (j in 1:ncol(p_Fisher)) { # Celltypes
      Q_Cell = length(which(df$Quadrants == Quadrants[i] & df[,k + j] == "specific"))
      Q_nonCell = length(which(df$Quadrants == Quadrants[i] & df[,k + j] != "specific"))
      nonQ_Cell = length(which(df$Quadrants != Quadrants[i] & df[,k + j] == "specific"))
      nonQ_nonCell = length(which(df$Quadrants != Quadrants[i] & df[,k + j] != "specific"))
      (Res = fisher.test(matrix(c(Q_Cell,Q_nonCell,nonQ_Cell,nonQ_nonCell), nrow = 2)))
      p_Fisher[i,j] =  format(Res$p.value, scientific = T, digits = 2)
      OR_Fisher[i,j] = format(Res$estimate, scientific = F, digits = 2)
    } 
  }
  
  OR_Fisher = as.data.frame(apply(OR_Fisher, 2, as.numeric))
  p_Fisher = as.data.frame(apply(p_Fisher, 2, as.numeric))
  colnames(p_Fisher) = colnames(OR_Fisher) = sapply(colnames(df)[(k+1):(k+6)], function(x) substr(x, 4, nchar(x)))
  rownames(p_Fisher) = rownames(OR_Fisher) = Quadrants
  
  OR_Fisher$Quadrant = rownames(OR_Fisher)
  library(reshape2)
  OR_Fisher_melt = reshape2::melt(OR_Fisher, value.name = "Odds_ratio", variable.name = "Celltype")
  p_Fisher$Quadrant = rownames(p_Fisher)
  p_Fisher_melt = reshape2::melt(p_Fisher, value.name = "P_value", variable.name = "Celltype")
  Enrichment_df = left_join(OR_Fisher_melt, p_Fisher_melt)
  Enrichment_df$FDR = p.adjust(Enrichment_df$P_value, method = "fdr") # All 1
  Enrichment_df$FDR_chr = format(Enrichment_df$FDR, scientific = T, digits = 2)
  Enrichment_df$label_text = paste0(Enrichment_df$Odds_ratio, "\n(", Enrichment_df$FDR_chr, ")")
  Enrichment_df$label_text[which(Enrichment_df$Odds_ratio <= 1 | Enrichment_df$FDR > 0.1)] = NA
  
  Enrichment_df$Celltype = factor(Enrichment_df$Celltype, 
                                  levels = c("Glia","Endoth", "Microglia", "Oligo", "Astro", "Neuron"))
  return(Enrichment_df)
}

# run the function
Enrichment_df = pSI_enrichment(hits2)

# plot
pdf("CorLogFC_sigATACsigDEG_2kbpromoter_Quadrants_pSIenrichment.pdf", width = 6, height = 4)
Enrichment_df[Enrichment_df$Celltype != "Glia",] %>% 
  ggplot(aes(x = Celltype, y = Quadrant)) +
  geom_tile(aes(fill = Odds_ratio)) + # , color = "black", size = 2
  scale_fill_gradient2(low = "white", high = "red", name = "Odds ratio", midpoint = 1, na.value = "red") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"), #panel.background = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank(),
        axis.text = element_text(size = 12)) +
  geom_text(aes(label=label_text)) +
  ggtitle("Enrichment of cell-type specifically expressed genes\nin differential promoter ATAC-DEG pairs")
dev.off()

# Observation:
# 2kb promoter: Q1 is significantly enriched for microglia genes, Q3 for oligodendrocyte genes, Q4 for neuron genes.

# For the neuron genes in Q4, is there any other genes that share the ATAC peak?
tmp = hits2[hits2$Quadrants == "Q4",] # 76 pairs
tmp = tmp[which(tmp$gene_name %in% rownames(pSI_res)[!is.na(pSI_res$Neuron)]),] # 22 pairs

tmp2 = hits2[hits2$ATACpeak %in% tmp$ATACpeak & hits2$Quadrants != "Q4",] # 0 neuronal ATAC peak overlap a different DEG. What about non-DEGs, any genes?

promoter_gr = annoDF2GR(promoter[,-4])
hits3 = as.data.frame(findOverlaps(sigATAC_gr, promoter_gr))
hits3$ATACpeak = rownames(sigATAC)[hits3$queryHits]
hits3$gene_id = promoter$ensembl_gene_id[hits3$subjectHits]
hits3$ATAC_logFC = sigATAC$ATAC_logFC[hits3$queryHits]
hits3$gene_name = promoter$external_gene_name[hits3$subjectHits]
tmp3 = hits3[hits3$ATACpeak %in% tmp$ATACpeak & ! (hits3$gene_id %in% tmp$gene_id),] # 3 neuronal ATAC peaks overlap a different gene (not DEG): DPYSL2 and GOLM1
pSI_res[tmp3$gene_name,] # all NAs, these other genes are not specifically expressed in any cell-types

# This agrees with Gokul's paper in that neuronal DARs compensate DEG.


