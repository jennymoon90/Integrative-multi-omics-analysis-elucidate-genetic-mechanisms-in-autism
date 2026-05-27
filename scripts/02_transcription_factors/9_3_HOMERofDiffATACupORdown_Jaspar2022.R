
### 1) Prepare bed files for HOMER TF enrichment
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/")

DIR="~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/"
lnames = load(paste0(DIR,"DiffATAC.rda")) # DiffATAC

for (i in c(1,2)) {
  if (i == 1) {
    fdr_t = 0.05; logfc_t = 0
  } else {
    fdr_t = 0.1; logfc_t = 0.2
  }
  
  ## Filter to differential ATAC peaks
  ATAC_up = DiffATAC[which(DiffATAC$ATAC_logFC > logfc_t & DiffATAC$ATAC_FDR < fdr_t),]
  ATAC_down = DiffATAC[which(DiffATAC$ATAC_logFC < -logfc_t & DiffATAC$ATAC_FDR < fdr_t),]
  
  ## Write to bed files
  write.table(ATAC_up, file = paste0("upATAC_fdrthreshold",fdr_t,"logfcthreshold",logfc_t,".bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(ATAC_down, file = paste0("downATAC_fdrthreshold",fdr_t,"logfcthreshold",logfc_t,".bed"), quote = F, row.names = F, col.names = F, sep = "\t")
}

## Check the number of bed files
files = list.files()
length(files) # 4, correct.

### Run HOMER: 9_4_HOMER_MotifScan_in_DiffATACupORdown_Jaspar2022.sh
# For motif, used jaspar2022_BrExpTF.motifs generated using 14_10_TFBS_enrichment_in_4dif_categories_of_gene_BasedOn_Nloop_TPM_Jaspar2022.R and 14_11_HOMER_scan_motif_in_promoterATAC_4categories_Jaspar2022.sh (this is the JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt -> restrict to human brain expressed TFs)

### TF binding motif enrichment plot
rm(list = ls())
#setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/Jaspar2022_results/")
setwd("/Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/Jaspar2022_results/")

(motif_files = list.files(pattern = "*_knownResults.txt")) # 2

for (motif_file in motif_files) {
  (new_name = substr(motif_file,1,nchar(motif_file) - 4))
  
  new_name = gsub("fdrthreshold0.05logfcthreshold0", "threshold1", new_name)
  new_name = gsub("fdrthreshold0.1logfcthreshold0.2", "threshold2", new_name)
  new_name = gsub("_Jaspar2022", "", new_name)
  
  new_motif = read.table(paste0(motif_file)) # "HOMER_TFenrichment_results/", 
  colnames(new_motif) = c("Consensus", "P_value", "logP_value", "q_value_Benjamini", 
                          "Number_of_Target_Sequences_with_Motif","Percent_of_Target_Sequences_with_Motif",
                          "Number_of_Background_Sequences_with_Motif","Percent_of_Background_Sequences_with_Motif")
  new_motif$Motif_name = rownames(new_motif)
  #new_motif = unique(new_motif)
  assign(new_name, new_motif)
}
# See that in fdr < 0.05 Q3 is enriched for RFX4 and fdr < 0.1 & |logFC| > 0.2 Q1 is enriched for RORC (q < 0.1). RFX4 is an astrocyte-marker expressed gene.

# Generate a matrix comparing the q-value of each motif.
ref = downATAC_threshold1_knownResults
#downATAC_threshold2_knownResults = downATAC_threshold2_knownResults[match(rownames(ref), rownames(downATAC_threshold2_knownResults)),]
upATAC_threshold1_knownResults = upATAC_threshold1_knownResults[match(rownames(ref), rownames(upATAC_threshold1_knownResults)),]
#upATAC_threshold2_knownResults = upATAC_threshold2_knownResults[match(rownames(ref), rownames(upATAC_threshold2_knownResults)),]

Motif_comparison_main = 
  data_frame(Motif_name = ref$Motif_name,
             Consensus = ref$Consensus,
             
             pval_down_t1 = downATAC_threshold1_knownResults$P_value,
             pval_up_t1 = upATAC_threshold1_knownResults$P_value,
             #pval_down_t2 = downATAC_threshold2_knownResults$P_value,
             #pval_up_t2 = upATAC_threshold2_knownResults$P_value,

             qval_down_t1 = downATAC_threshold1_knownResults$q_value_Benjamini,
             qval_up_t1 = upATAC_threshold1_knownResults$q_value_Benjamini,
             #qval_down_t2 = downATAC_threshold2_knownResults$q_value_Benjamini,
             #qval_up_t2 = upATAC_threshold2_knownResults$q_value_Benjamini,
  )

idx = which(apply(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "qval"))], 1, function(x) any(x < 0.05)))
Motif_comparison_main = 
  Motif_comparison_main[idx,] # only FDR<0.05 are shown

Motif_comparison_main$minus_log_qval_down_t1 = -log10(Motif_comparison_main$qval_down_t1)
#Motif_comparison_main$minus_log_qval_down_t2 = -log10(Motif_comparison_main$qval_down_t2)
Motif_comparison_main$minus_log_qval_up_t1 = -log10(Motif_comparison_main$qval_up_t1)
#Motif_comparison_main$minus_log_qval_up_t2 = -log10(Motif_comparison_main$qval_up_t2)

Motif_comparison_main$Motif_abbr = sapply(Motif_comparison_main$Motif_name, function(x) substr(x, 1, which(strsplit(x, "")[[1]] == "/")[1] -1))

library(reshape2)
Motif_comparison_main_melt = melt(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "minus_log_qval") | colnames(Motif_comparison_main) == "Motif_abbr")])
colnames(Motif_comparison_main_melt)[3] = c("significance")
Motif_comparison_main_melt$variable = as.character(Motif_comparison_main_melt$variable)
Motif_comparison_main_melt$variable = gsub("minus_log_qval_", "", Motif_comparison_main_melt$variable)
Motif_comparison_main_melt$variable = gsub("t1", "threshold1", Motif_comparison_main_melt$variable)
#Motif_comparison_main_melt$variable = gsub("t2", "threshold2", Motif_comparison_main_melt$variable)

# for visibility in the plot
#Motif_comparison_main_melt$significance[which(Motif_comparison_main_melt$significance == 0)] = 0.01

Motif_comparison_main$sig_n = apply(Motif_comparison_main[,which(grepl("^qval", colnames(Motif_comparison_main)))], 1, function(x) length(which(x < 0.05)))
table(Motif_comparison_main$sig_n) # 3 motifs are enriched in both up and down-reg ATAC peaks, 35 motifs are uniquely enriched in either up or down-reg ATAC peaks. old: 2 motifs are enriched in all ATAC categories, 4 motifs are enriched in 3 ATAC categories， 25 motifs are enriched in two ATAC categories, 16 motifs are enriched in 1 ATAC category.
Motif_comparison_main$sig_m = Motif_comparison_main$sig_n
#Motif_comparison_main$sig_m[Motif_comparison_main$sig_m == 2] = 1

Motif_comparison_main = Motif_comparison_main %>%
  arrange(dplyr::desc(sig_m), dplyr::desc(minus_log_qval_down_t1), minus_log_qval_up_t1, Motif_name)

duplicated(Motif_comparison_main$Motif_abbr) # BHLHE22 is duplicated as it has two motifs, but one is enriched in up-reg ATAC, while the other in down-reg ATAC. -> delete BHLHE22
Motif_comparison_main = Motif_comparison_main[Motif_comparison_main$Motif_abbr != "BHLHE22",]

Motif_comparison_main_melt$Motif_abbr =
  factor(Motif_comparison_main_melt$Motif_abbr,
         levels = Motif_comparison_main$Motif_abbr) 
Motif_comparison_main_melt = Motif_comparison_main_melt[!is.na(Motif_comparison_main_melt$Motif_abbr),]

save(Motif_comparison_main_melt, file = "Motif_Enrichment_in_DiffATACupORdown_HOMERq005_Jaspar2022.rda")

# plot
# pdf("Motif_Enrichment_in_DiffATACupORdown_HOMERq005_Jaspar2022.pdf", 
#     width = 10, height = 4.5)
# Motif_comparison_main_melt %>%
#   ggplot(aes(x = Motif_abbr, y = significance, fill = variable)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   ggtitle("Enrichment of TFBS in ASD differential ATAC peaks (Only TF motifs with FDR<0.05 are shown)") +
#   geom_hline(aes(yintercept=-log10(0.05), linetype="FDR = 0.05"), color = "black", size=0.3) +
#   ylab("-log10(FDR)") +
#   xlab("Transcription factors") +
#   labs(fill = "ASD differential ATAC peaks") +
#   theme_bw() +
#   theme(#axis.text=element_text(size=20),
#     axis.text.x=element_text(angle = 45, hjust = 1),
#     #axis.title=element_text(size=20),
#     #legend.text = element_text(size=20),
#     #title =element_text(size=20),
#     #plot.title = element_text(hjust = 0.5)
#     legend.position="top"
#   ) +
#   scale_linetype_manual(name = "Enrichment", values = c("FDR = 0.05" = 2)) + #, "FDR = 0.2" = 3
#   scale_fill_manual(values = c("down_threshold1" = "deepskyblue", "down_threshold2" = "dodgerblue", "up_threshold1" = "#F8766D", "up_threshold2" = "pink")) # purple #C77CFF
# dev.off()

