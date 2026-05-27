## Differential definition: FDR < 0.05
## Proximal: TSS - 2kb to + 100bp, Distal: ATAC peaks at the distal end of promoter Hi-C loops
## Distal ATAC target genes: significant promoter-distal ATAC cor. from 10_3_8_LinkATACtoGene_v5RPKMcqnATACcorInPLs_PromoterATACTSSminus2kbplus100bp.R

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)
library(plyr)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/5_6_v3_v5_correct_PromoterTssMinus2kbPlus100bp_DistalTargetGeneByPromoterDistalATACsigCor/")

###################################################################
## Classify diffATAC into Promoter/Distal_HiC/Distal_nonHiC, and ##
## Quantify up- vs. down-reg ATAC peaks in each category         ##
###################################################################

## Get promoter regions (tss - 2kb to +100bp)
# done in 5_6_v3_number_of_ATAC_peaks_upORdown_proximalVShicdistal_v4relativeAll4cat_CorrectPromoterMinus2kbToPlus100bp.R
load("../promoter_tss_minus2kb_to_plus100bp.rda") # promoter
promoter_gr = annoDF2GR(promoter[,-4])

## Get promoter Hi-C loops
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/PL_by_sample_datasets/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_rmOutlierTechRep.rda") # Bulk_logCPM, NeuNp_logCPM, NeuNn_logCPM

promoter$tss = ifelse(promoter$strand == "+", promoter$end - 100, promoter$start + 100)
table(promoter$strand) # 11737 +, 11319 -, 0 *
promoter$mid = ifelse(promoter$strand == "+", floor(promoter$tss/10e3)*10e3 + 5e3, floor(promoter$tss/10e3)*10e3 - 5e3)
promoter$mid = paste0(promoter$chr, "_", promoter$mid)

get_promoter_HiCloop_location = function(logCPM) {
  logCPM = as.data.frame(str_split_fixed(rownames(logCPM), "_", 3))
  logCPM$V2 = as.integer(logCPM$V2); logCPM$V3 = as.integer(logCPM$V3)
  logCPM$mid1 = paste0(logCPM$V1, "_", logCPM$V2)
  logCPM$mid2 = paste0(logCPM$V1, "_", logCPM$V3)
  
  # intersect each loop end with DEG promoter
  idx1 = which(logCPM$mid1 %in% promoter$mid)
  idx2 = which(logCPM$mid2 %in% promoter$mid)
  
  df1 = logCPM[idx1,]; colnames(df1)[4:5] = c("mid","distal_mid")
  df2 = logCPM[idx2,]; colnames(df2)[4:5] = c("distal_mid","mid")
  
  df1 = left_join(df1, promoter[,c("mid", "ensembl_gene_id", "external_gene_name")])
  df2 = left_join(df2, promoter[,c("mid", "ensembl_gene_id", "external_gene_name")])
  
  # add distal end start and end information
  df1$start = df1$V3 - 5e3; df1$end = df1$V3 + 5e3
  df2$start = df2$V2 - 5e3; df2$end = df2$V2 + 5e3
  
  df2 = df2[,colnames(df1)]
  df = unique(rbind(df1, df2))
  colnames(df)[1:3] = c("chr", "fragmentMid1", "fragmentMid2")
  
  return(df)
}

Bulk_logCPM = get_promoter_HiCloop_location(Bulk_logCPM)
Bulk_HiCdistal_gr = annoDF2GR(Bulk_logCPM)

## Proximal and distal ATAC peak % up or down reg plot
for (i in c(1,2)) {
  ## Load differential ATAC results
  if (i == 1) {
    fdr_t = 0.05; logfc_t = 0
  } else {
    fdr_t = 0.1; logfc_t = 0.2
  }
  lnames = load(paste0("../../5_DiffATAC/MAdf_FDRthreshold",fdr_t,"_logFCthreshold",logfc_t,".rda")) # MA_df
  
  ## filter for differential ATAC peaks
  sigATAC = MA_df[which(MA_df$col != "lightgrey" & !is.na(MA_df$ATAC_logFC)),]
  sigATAC_gr = annoDF2GR(sigATAC)
  
  ## Proximal and distal ATAC peaks
  # Promoter diff ATAC
  hits = as.data.frame(findOverlaps(sigATAC_gr, promoter_gr))
  idx_proximal = unique(hits$queryHits)
  sigATAC$type = "Distal_nonHiC"
  sigATAC$type[idx_proximal] = "Promoter"
  sigATAC_promoter = sigATAC[idx_proximal,]
  
  # Which genes are at promoter diff ATAC?
  hits$ATACpeak = rownames(sigATAC)[hits$queryHits]
  hits = cbind(hits, sigATAC[hits$queryHits,4:6])
  hits = cbind(hits, promoter[hits$subjectHits, 8:7])
  
  # promoter Hi-C diff ATAC
  hits2 = as.data.frame(findOverlaps(sigATAC_gr, Bulk_HiCdistal_gr))
  length(unique(hits2$queryHits)) # 908 distal ATAC peaks at bulk promoter Hi-C loop end
  distal_ATAC = unique(rownames(sigATAC)[hits2$queryHits]) # 908
  idx_distalHiC = which(rownames(sigATAC) %in% distal_ATAC)
  idx_distalHiC_nonPromoter = setdiff(idx_distalHiC, idx_proximal) # 639. meaning 269 ATAC peaks are both promoter and distal Hi-C end ATAC peak
  #sigATAC$type[idx_distalHiC] = "Distal_HiC"
  sigATAC$type[idx_distalHiC_nonPromoter] = "Distal_HiC"
  sigATAC_distalHiC = sigATAC[idx_distalHiC,]; sigATAC_distalHiC$type = "Distal_HiC" # include ATAC that are both promoter and distal (dual ATAC)
  
  # Which genes are at promoter Hi-C loop distal diff ATAC?
  hits2$ATACpeak = rownames(sigATAC)[hits2$queryHits] # includes dual ATAC
  hits2 = cbind(hits2, sigATAC[hits2$queryHits,4:6])
  #hits2 = cbind(hits2, loop = rownames(Bulk_logCPM)[hits2$subjectHits], Bulk_logCPM[hits2$subjectHits, 6:7])
  hits2 = cbind(hits2, loop = paste0(Bulk_logCPM$chr[hits2$subjectHits], "_", Bulk_logCPM$fragmentMid1[hits2$subjectHits], "_", Bulk_logCPM$fragmentMid2[hits2$subjectHits]), Bulk_logCPM[hits2$subjectHits, 6:7])
  # the distal genes right now are wo considering proximal-distal ATAC correlation
  
  # some ATAC peaks may be both promoter and Distal HiC
  sigATAC_labelType = sigATAC[sigATAC$type == "Distal_nonHiC",]
  sigATAC_labelType = rbind(sigATAC_labelType, sigATAC_promoter)
  sigATAC_labelType = rbind(sigATAC_labelType, sigATAC_distalHiC)
  
  length(unique(rownames(sigATAC_labelType))) # 5302, duplicated rownames were added 1 to the rownames. this equals to nrow(sigATAC_labelType)
  sigATAC_labelType$ATAC_peak = paste0(sigATAC_labelType$chr, "_", sigATAC_labelType$start, "_", sigATAC_labelType$end)
  length(unique(sigATAC_labelType$ATAC_peak)) # 5033, good
  
  promoter_sigATAC = hits; promoterHiCdistal_sigATAC = hits2
  save(promoter_sigATAC, promoterHiCdistal_sigATAC, sigATAC_labelType, file = paste0("promoter_minus2kbplus100bp_proximal_HiCdistal_diffATAC_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".rda"))
  # In 5_7_GO_diffATAC.R, I will print out the sorted diffATAC bed files
  
  ## Plot percentage up/down-reg
  counts = ddply(sigATAC_labelType, .(sigATAC_labelType$type, sigATAC_labelType$col), nrow) # ATAC as both promoter and Distal_HiC are counted twice
  colnames(counts) = c("ATAC_type", "Dysregulation_in_ASD", "Freq")
  counts$Dysregulation_in_ASD[which(counts$Dysregulation_in_ASD == "blue")] = "Down-reg"
  counts$Dysregulation_in_ASD[which(counts$Dysregulation_in_ASD == "red")] = "Up-reg"
  counts$Percent = NA
  counts$Percent[1] = counts$Freq[1]/sum(counts$Freq[1:2])
  counts$Percent[2] = 1 - counts$Percent[1]
  counts$Percent[3] = counts$Freq[3]/sum(counts$Freq[3:4])
  counts$Percent[4] = 1 - counts$Percent[3]
  counts$Percent[5] = counts$Freq[5]/sum(counts$Freq[5:6])
  counts$Percent[6] = 1 - counts$Percent[5]
  counts$Dysregulation_in_ASD = factor(counts$Dysregulation_in_ASD, levels = c("Up-reg", "Down-reg"))
  
  p = counts %>% 
    ggplot(aes(x = ATAC_type, y = Percent, fill = Dysregulation_in_ASD)) +
    geom_col(position = position_stack()) +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5)) +
    theme_bw() +
    coord_flip() +
    labs(fill = "") + # "H3K27Ac peaks"
    xlab("") +
    scale_y_continuous("", breaks = c(0,1), labels = c(0,sum(counts$Freq[1:2])), sec.axis = sec_axis(~.* sum(counts$Freq[5:6]), breaks = c(0,sum(counts$Freq[5:6])), name = "Number of ATAC peaks dysregulated in ASD")) #+
    # theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(t = -1.2, unit = "cm") #, 
    #       #plot.title = element_text(hjust = -20)
    # )
  
  pdf(paste0("Number_of_UpDownReg_ATACpeaks_ProximalVsHiCDistal_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,"_promoterTssMinus2kbToPlus100bp.pdf"), height = 3.5, width = 7)
  print(p)
  dev.off()
  
  save(counts, file = paste0("Number_of_UpDownReg_ATACpeaks_ProximalVsHiCDistal_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,"_promoterTssMinus2kbToPlus100bp.rda"))
}

# Observation:
# Proximal ATAC peaks are 1.9-2.6 x more likely to be up-regulated.
# Distal_HiC ATAC peaks are 1.95-2.75 x more likely to be down-regulated.

##########################################################################
## Re-define distal ATAC target genes based on promoter-distal ATAC cor ##
##########################################################################
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/5_6_v3_v5_correct_PromoterTssMinus2kbPlus100bp_DistalTargetGeneByPromoterDistalATACsigCor/")

lnames = load("../../10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda") # sig_CorATAC, sig_CorATAC_uniq, promATAC
library(readxl)
Jill_DEG = readxl::read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

sig_CorATAC_uniq$external_gene_name = Jill_DEG$external_gene_name[match(sig_CorATAC_uniq$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

for (i in c(1,2)) {
  ## Load promoter diff ATAC and promoter loop distal diff ATAC results
  if (i == 1) {
    fdr_t = 0.05; logfc_t = 0
  } else {
    fdr_t = 0.1; logfc_t = 0.2
  }
  lnames = load(paste0("promoter_minus2kbplus100bp_proximal_HiCdistal_diffATAC_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".rda")) # promoter_sigATAC, promoterHiCdistal_sigATAC, sigATAC_labelType
  distal_sigATAC = promoterHiCdistal_sigATAC; rm(promoterHiCdistal_sigATAC)
  length(unique(distal_sigATAC$ATACpeak)) # 908 for fdr < 0.05
  length(unique(promoter_sigATAC$ATACpeak)) # 844 for fdr < 0.05
  # ATAC peak as both promoter and distal_hic are counted twice.
  
  colnames(distal_sigATAC)[3] = "distal_ATAC"
  # distal_sigATAC_true = left_join(distal_sigATAC[,3:6], sig_CorATAC_uniq[,c(1,2,7)]) # 4287, this includes target genes within 30kb of the ATAC. I just want Hi-C loop target genes with significant proximal-distal ATAC cor for the downstream enrichment analysis.
  distal_sigATAC_true = inner_join(distal_sigATAC, sig_CorATAC_uniq[,c(1,2,7)]) # 284, 14.5% of the original dataframe. Great
  
  save(promoter_sigATAC, distal_sigATAC_true, sigATAC_labelType, file = paste0("DiffATAC_PromoterOrHicdistal_wTargetGenes_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,"_PromoterTssMinus2kbToPlus100bp_DistalTargetBasedOnATACsigCor.rda"))
}

#################################################################
## Celltype, SFARI, and GO enrichment of up- vs. down-reg ATAC ##
## at Promoter vs. Distal_HiC                                  ##
#################################################################

rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/5_6_v3_v5_correct_PromoterTssMinus2kbPlus100bp_DistalTargetGeneByPromoterDistalATACsigCor/")

# Load pSI data
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/Zhang_2015_pSI/Zhang2015pSI_usingMatureOfAllAges.rda") 
colnames(pSI_res)

# pSI enrichment function - modified for type instead of Quadrants
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
  
  p_Fisher = OR_Fisher = as.data.frame(matrix(nrow = length(type), ncol = 6))
  type = unique(df$type)
  colnames(df)
  k = which(colnames(df) == "SI_Glia") - 1
  
  for (i in 1:length(type)) { # type
    for (j in 1:ncol(p_Fisher)) { # Celltypes
      Q_Cell = length(which(df$type == type[i] & df[,k + j] == "specific"))
      Q_nonCell = length(which(df$type == type[i] & df[,k + j] != "specific"))
      nonQ_Cell = length(which(df$type != type[i] & df[,k + j] == "specific"))
      nonQ_nonCell = length(which(df$type != type[i] & df[,k + j] != "specific"))
      (Res = fisher.test(matrix(c(Q_Cell,Q_nonCell,nonQ_Cell,nonQ_nonCell), nrow = 2)))
      p_Fisher[i,j] =  format(Res$p.value, scientific = T, digits = 2)
      OR_Fisher[i,j] = format(Res$estimate, scientific = F, digits = 2)
    } 
  }
  
  OR_Fisher = as.data.frame(apply(OR_Fisher, 2, as.numeric))
  p_Fisher = as.data.frame(apply(p_Fisher, 2, as.numeric))
  colnames(p_Fisher) = colnames(OR_Fisher) = sapply(colnames(df)[(k+1):(k+6)], function(x) substr(x, 4, nchar(x)))
  rownames(p_Fisher) = rownames(OR_Fisher) = type
  
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

# SFARI enrichment function - modified
SFARI = read.csv("~/Documents/Documents/Geschwind_lab/LAB/Database_download/SFARI_human_genes/SFARI-Gene_genes_01-11-2022release_03-03-2022export.csv")
SFARI_enrichment = function(df) {
  type = unique(df$type)
  p_Fisher = OR_Fisher = as.data.frame(matrix(nrow = length(type), ncol = 1))
  
  for (i in 1:length(type)) { # type
    Q_Cell = length(which(df$type == type[i] & ! is.na(df$SFARI_score)))
    Q_nonCell = length(which(df$type == type[i] & is.na(df$SFARI_score)))
    nonQ_Cell = length(which(df$type != type[i] & ! is.na(df$SFARI_score)))
    nonQ_nonCell = length(which(df$type != type[i] & is.na(df$SFARI_score)))
    (Res = fisher.test(matrix(c(Q_Cell,Q_nonCell,nonQ_Cell,nonQ_nonCell), nrow = 2)))
    p_Fisher[i,1] =  format(Res$p.value, scientific = T, digits = 2)
    OR_Fisher[i,1] = format(Res$estimate, scientific = F, digits = 2)
  }
  
  OR_Fisher = as.data.frame(apply(OR_Fisher, 2, as.numeric))
  p_Fisher = as.data.frame(apply(p_Fisher, 2, as.numeric))
  colnames(p_Fisher) = colnames(OR_Fisher) = c("SFARI gene") 
  rownames(p_Fisher) = rownames(OR_Fisher) = type
  
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
  
  return(Enrichment_df)
}

for (i in c(1,2)) {
  ## Load promoter diff ATAC and promoter loop distal diff ATAC results
  if (i == 1) {
    fdr_t = 0.05; logfc_t = 0
  } else {
    fdr_t = 0.1; logfc_t = 0.2
  }
  lnames = load(paste0("DiffATAC_PromoterOrHicdistal_wTargetGenes_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,"_PromoterTssMinus2kbToPlus100bp_DistalTargetBasedOnATACsigCor.rda")) # promoter_sigATAC, distal_sigATAC_true, sigATAC_labelType
  
  ## Genes at promoter up/down-regulated ATAC peaks
  promoter_up_genes = unique(promoter_sigATAC[promoter_sigATAC$ATAC_logFC > 0, 7:8]) # 718 genes
  promoter_down_genes = unique(promoter_sigATAC[promoter_sigATAC$ATAC_logFC < 0, 7:8]) # 279 genes
  
  ## Genes at promoter loop distal up/down-regulated ATAC peaks
  distal_up_genes = unique(distal_sigATAC_true[distal_sigATAC_true$ATAC_logFC > 0, 8:9]) # 386 genes
  distal_down_genes = unique(distal_sigATAC_true[distal_sigATAC_true$ATAC_logFC < 0, 8:9]) # 982 genes
  
  ## How many genes intersect?
  overlap_genes = intersect(promoter_up_genes$external_gene_name, distal_down_genes$external_gene_name) # Empty
  print(paste0("fdr < ", fdr_t, ", |logFC| > ",logfc_t, ", number of overlapped genes: ", length(overlap_genes)))
  
  ## Are these two gene sets enriched in any cell type?
  promoter_up_genes$type = "promoter_up"
  promoter_down_genes$type = "promoter_down"
  distal_up_genes$type = "distal_up"
  distal_down_genes$type = "distal_down"
  df = rbind(promoter_up_genes, promoter_down_genes, distal_up_genes, distal_down_genes)
  colnames(df)[which(colnames(df) == "external_gene_name")] = "gene_name"
  
  df = df[complete.cases(df),] # 1167 -> 1165 rows with fdr < 0.05
  
  # run the function
  Enrichment_df = pSI_enrichment(df)
  save(Enrichment_df, file = paste0("pSIenrichment_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".rda"))
  
  # pSI plot
  pl = Enrichment_df[Enrichment_df$Celltype != "Glia",] %>% 
    ggplot(aes(x = Celltype, y = Quadrant)) +
    geom_tile(aes(fill = Odds_ratio)) + # , color = "black", size = 2
    scale_fill_gradient2(low = "white", high = "red", name = "Odds ratio", midpoint = 1, na.value = "red") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"), #panel.background = element_blank(), 
          axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text = element_text(size = 12)) +
    geom_text(aes(label=label_text)) +
    ggtitle("Enrichment of cell-type specifically expressed genes\nin differential ATAC peaks")
  
  pdf(paste0("PromoterDiffATAC_promoterHiCdistalDiffATAC_pSIenrichment_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".pdf"), width = 6, height = 3)
  print(pl)
  dev.off()
  
  ## SFARI enrichment
  df$SFARI_score = SFARI$gene.score[match(df$ensembl_gene_id, SFARI$ensembl.id)]
  Enrichment_df = SFARI_enrichment(df)
  
  save(Enrichment_df, file = paste0("SFARIenrichment_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".rda"))
  
  pl2 = Enrichment_df %>% 
    ggplot(aes(x = Celltype, y = Quadrant)) +
    geom_tile(aes(fill = Odds_ratio)) + # , color = "black", size = 2
    scale_fill_gradient2(low = "white", high = "red", name = "Odds ratio", midpoint = 1, na.value = "red") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"), #panel.background = element_blank(), 
          axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text = element_text(size = 12)) +
    geom_text(aes(label=label_text)) +
    ggtitle("Enrichment of SFARI genes\nin differential ATAC peaks")
  
  pdf(paste0("PromoterDiffATAC_promoterHiCdistalDiffATAC_SFARIenrichment_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,".pdf"), width = 5, height = 3)
  print(pl2)
  dev.off()
}

# Observation:
# In either case (fdr<0.05 or [fdr<0.1 & |logFC|>0.2]), distal_up ATAC peaks are enriched for SFARI genes, only significant with [fdr<0.1 & |logFC|>0.2].

## Go enrichment
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)
library(gprofiler2)
library(ggplot2)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/5_6_v3_v5_correct_PromoterTssMinus2kbPlus100bp_DistalTargetGeneByPromoterDistalATACsigCor/")

fdr_t = 0.05; logfc_t = 0
lnames = load(paste0("DiffATAC_PromoterOrHicdistal_wTargetGenes_ATACfdrthreshold",fdr_t,"_logfcthreshold",logfc_t,"_PromoterTssMinus2kbToPlus100bp_DistalTargetBasedOnATACsigCor.rda")) # promoter_sigATAC, distal_sigATAC_true, sigATAC_labelType

# Quick check
table(sigATAC_labelType$type) # 844 Promoter, 908 Distal_HiC, 3550 Distal_nonHiC
colnames(distal_sigATAC_true)[3] = "ATACpeak"
length(unique(promoter_sigATAC$ATACpeak)) # 844
length(unique(distal_sigATAC_true$ATACpeak)) # 209 - this is the subset of Distal_HiC ATAC peaks that have target genes linked by significant proximal-distal ATAC pairs

Promoter_down_geneid = unique(promoter_sigATAC$ensembl_gene_id[promoter_sigATAC$ATAC_logFC < 0 & ! is.na(promoter_sigATAC$ensembl_gene_id)]) # 245 genes
Promoter_up_geneid = unique(promoter_sigATAC$ensembl_gene_id[promoter_sigATAC$ATAC_logFC > 0 & ! is.na(promoter_sigATAC$ensembl_gene_id)]) # 864 genes
DistalHiC_down_geneid = unique(distal_sigATAC_true$ensembl_gene_id[distal_sigATAC_true$ATAC_logFC < 0 & ! is.na(distal_sigATAC_true$ensembl_gene_id)]) # 158 genes
DistalHiC_up_geneid = unique(distal_sigATAC_true$ensembl_gene_id[distal_sigATAC_true$ATAC_logFC > 0 & ! is.na(distal_sigATAC_true$ensembl_gene_id)]) # 78 genes

# Check the overlap across the 4 gene sets
library(VennDiagram)
venn.diagram(
  x = list(Promoter_down_geneid, Promoter_up_geneid, DistalHiC_down_geneid, DistalHiC_up_geneid),
  category.names = c("Promoter_\nDowm" , "Promoter_\nUp" , "DistalHiC_Down", "DistalHiC_Up"),
  filename = 'DiffATAC_PromoterOrDistalHiC_geneid_VennDiagram.png',
  output=TRUE,
  cat.fontfamily = "sans"
)
# Distal_down and Promoter_down share 21, Distal_up and Promoter_up share 11. Otherwise no overlap. Great!

## GO analysis and plot
library(readxl)
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

gost_PromoterDown_genes = gost(query = Promoter_down_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_PromoterDown_genes = gost_PromoterDown_genes$result # axon ensheathment in central nervous system,etc

gost_PromoterUp_genes = gost(query = Promoter_up_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_PromoterUp_genes = gost_PromoterUp_genes$result # response to stimulus, etc

gost_DistalHiCDown_genes = gost(query = DistalHiC_down_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_DistalHiCDown_genes = gost_DistalHiCDown_genes$result # nervous system development, neuron projection morphogenesis, etc

gost_DistalHiCUp_genes = gost(query = DistalHiC_up_geneid, organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_DistalHiCUp_genes = gost_DistalHiCUp_genes$result # negative regulation of neuron death, etc

# Plot
GO_PromoterDown_genes$ATAC_type = "Promoter_Down"
GO_PromoterUp_genes$ATAC_type = "Promoter_Up"
GO_DistalHiCDown_genes$ATAC_type = "DistalHiC_Down"
GO_DistalHiCUp_genes$ATAC_type = "DistalHiC_Up"

GO_df = rbind(GO_PromoterDown_genes, GO_PromoterUp_genes, GO_DistalHiCDown_genes, GO_DistalHiCUp_genes)
length(unique(GO_df$term_name)) # 741 unique GO term names
hist(GO_df$term_size, breaks = 100) # mostly under 6000
abline(v = 4500)
length(unique(GO_df$term_name[GO_df$term_size < 6e3])) # 726 unique GO term names
length(unique(GO_df$term_name[GO_df$term_size < 5000])) # 705 unique GO term names
length(unique(GO_df$term_name[GO_df$term_size < 4500])) # 699 unique GO term names

GO_df$minusLog10Pval = -log10(GO_df$p_value)
GO_df$fdr = p.adjust(GO_df$p_value, method = "fdr")
GO_df$minusLog10FDR = -log10(GO_df$fdr)
#Terms = unique(GO_df$term_name[GO_df$p_value < 1e-8 & GO_df$term_size < 6e3]) # 21 terms

# Take top 7 terms (limit to FDR < 0.05) of each 4 GO term set after removing term.size > 5000 terms
require(data.table)
d = data.table(GO_df[GO_df$term_size < 4.5e3 & GO_df$fdr < 0.05,], key = "ATAC_type")
GO_df_selected = d[, head(.SD, 8), by = ATAC_type] # tried top 5-15 terms of each ATAC_type
(Terms = unique(GO_df_selected$term_name)) # term_size < 6000 & top 15 of each cat: 46 terms, term_size < 5000 & top 10 of each cat: 36 terms, term_size < 4500 & top 8 of each cat: 30 terms (nice terms), term_size < 4500 & top 5 of each cat: 20 terms (nice size for plot, but a bit restricted), yes all typical terms for each category are included!
GO_df_selectedTerms = GO_df[GO_df$term_name %in% Terms & GO_df$fdr < 0.05,] # 95 -> 75 -> 61 -> 40 rows, including p-val of all 4 categories for the union of top GO terms
GO_df_selectedTerms$minFDR = sapply(GO_df_selectedTerms$term_name, function(x) min(GO_df_selectedTerms$fdr[GO_df_selectedTerms$term_name == x]))

GO_df_selectedTerms = GO_df_selectedTerms %>%
  # arrange(desc(ATAC_type), desc(minusLog10FDR))
  arrange(minFDR)
# unique(GO_df_selectedTerms$ATAC_type) # "Promoter_Up"    "Promoter_Down"  "DistalHiC_Up"   "DistalHiC_Down"
GO_df_selectedTerms$term_name = factor(GO_df_selectedTerms$term_name, levels = rev(unique(GO_df_selectedTerms$term_name)))

save(GO_df, GO_df_selectedTerms, Terms, file = "DiffATAC_PromoterOrDistalHiC_geneGOenrichment.rda") # Terms_levels

pdf("DiffATAC_PromoterOrDistalHiC_geneGOenrichment.pdf", width = 8, height = 5)
GO_df_selectedTerms %>%
  ggplot() +
  geom_point(aes(x = term_name, y = minusLog10FDR, col = ATAC_type), alpha = 0.8) + # , size = term_size
  geom_hline(aes(yintercept = -log10(0.05), linetype = "FDR = 0.05")) + # linetype = "-log10(0.05)"
  theme_bw() +
  coord_flip() + 
  scale_linetype_manual(name = "Theshold", values = 2) +
  ylab("-log10(FDR)") +
  #xlab("GO terms") +
  ggtitle("GO enrichment of proximal and distal genes of differential ATAC peaks") +
  theme(plot.title = element_text(hjust = 0.65),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(color = rev(c(rep("black", length(i1)), rep("blue", length(i2)), rep("black", length(i3) + 1), rep("blue", length(i4)), rep("black", length(i5) + 2), "blue")))
        #legend.position = "bottom"
        )
dev.off()
# Promoter_Up genes show very significant p-val for around 50% of the GO terms, leading terms include response to chemical, etc. 
# All categories are enriched for cell projection/neuron terms at comparable levels. 
# Promoter_Down genes show specific enrichment for central nervous system myelination, neuron projection morphogenesis, etc. 
# DistalHiC_Down enriched terms are mostly related to neuron/cell projection. 
# DistalHiC_Up show specific enrichment for negative regulation of neuron death.
# May Need to plot the top terms separately for each category, or maybe not. Maybe better to separate top half and bottom half to separate scale

GO_df_selectedTerms = GO_df_selectedTerms %>%
  arrange(fdr)
colnames(GO_df_selectedTerms)
idx = which(duplicated(GO_df_selectedTerms[,c("term_name", "minFDR")]))
tmp = GO_df_selectedTerms[-idx,]
idx = which(tmp$ATAC_type != "Promoter_Up")
Terms1 = as.character(tmp$term_name[1:(min(idx)-1)])
Terms2 = as.character(tmp$term_name[min(idx):nrow(tmp)])
GO_df_selectedTerms$TermHalf = ifelse(GO_df_selectedTerms$term_name %in% Terms1, 1, 2)

save(GO_df_selectedTerms, file = "DiffATAC_PromoterOrDistalHiC_geneGOenrichment_TwoScales.rda")

pdf("DiffATAC_PromoterOrDistalHiC_geneGOenrichment_TwoScales.pdf", width = 8, height = 5)
GO_df_selectedTerms %>%
  ggplot() +
  geom_point(aes(x = term_name, y = minusLog10FDR, col = ATAC_type), alpha = 0.8) + # , size = term_size
  geom_hline(aes(yintercept = -log10(0.05), linetype = "FDR = 0.05")) + # linetype = "-log10(0.05)"
  theme_bw() +
  coord_flip() + 
  scale_linetype_manual(name = "Theshold", values = 2) +
  ylab("-log10(FDR)") +
  #xlab("GO terms") +
  ggtitle("GO enrichment of proximal and distal genes of differential ATAC peaks") +
  theme(plot.title = element_text(hjust = 0.65),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(color = rev(c(rep("black", length(i1)), rep("blue", length(i2)), rep("black", length(i3) + 1), rep("blue", length(i4)), rep("black", length(i5) + 2), "blue")))
        #legend.position = "bottom"
  ) +
  facet_wrap(TermHalf~., scales = "free", nrow = 2) +
  theme(strip.text.x = element_blank())
dev.off()

