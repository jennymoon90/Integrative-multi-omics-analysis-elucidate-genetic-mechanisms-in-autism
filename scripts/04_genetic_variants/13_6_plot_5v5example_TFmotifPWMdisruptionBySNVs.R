rm(list = ls())
library(ggplot2)
library(tidyverse)
# install.packages("ggseqlogo")
library(ggseqlogo)
library(MotifDb)
library(cowplot)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/13_SFARIgeneEnrichment/13_5_v4_TFBS_disrupted_by_SNVs/")

#### 1) Load the two examples

lnames = load("SNVs_LabelMotifbreakrEffect_ProbandOrSibling_within500bpOfDownregDistATAC_targetGenes_SFARIgeneLabeled.rda") # hits_Proband, hits_Sibling, df_proband4plotgardner, df_sibling4plotgardner
df1 = hits_Proband[hits_Proband$external_gene_name == "ADSL" & hits_Proband$alleleDiff < 0,]
df2 = hits_Proband[hits_Proband$external_gene_name == "MAP2K2" & hits_Proband$alleleDiff < 0,]

# ----
### How to plot TF binding sequence logo?
# Refer to https://omarwagih.github.io/ggseqlogo/
# data(ggseqlogo_sample)
# seqs_dna: sets of binding sites for 12 transcription factors obtained from FASTA files in JASPAR. This is represented as a named list of character vectors, where the names represent the JASPAR ID.
# pfms_dna: a list of position frequency matrices for four transcription factors obtained from JASPAR. This is represented as a list of matrices, where the names represent the JASPAR ID.
# seqs_aa: sets of kinase-substrate phosphorylation sites obtained from Wagih et al. This is represented as a named list of character vectors where the names represent the names of the kinases associated with the phosphosites.

# seqs_numeric = chartr('ATGC','1234', seqs_dna$MA0001.1)
# ggseqlogo(seqs_numeric, method='p', namespace=1:4) 

# ggseqlogo(seqs_dna$MA0001.1, method='p') # this is what I want
# ggseqlogo(pfms_dna$MA0018.2, method='p') # this is what I want

# (mf_names = paste0("Hsapiens-jaspar2022-", df1$TF, "-", df1$TFmotif_id)) # "Hsapiens-jaspar2022-SOX10-MA0442.2" "Hsapiens-jaspar2022-SOX4-MA0867.2" 
# mf = MotifDb[mf_names]
# pfm1 = mf@listData[[mf_names[1]]] # mf[[1]]
# ggseqlogo(pfm1, method='p') # great
# ----

#### 2) Get TF pfms
df = rbind(df1, df2)
(mf_names = paste0("Hsapiens-jaspar2022-", df$TF, "-", df$TFmotif_id)) # "Hsapiens-jaspar2022-SOX10-MA0442.2" "Hsapiens-jaspar2022-SOX4-MA0867.2" "Hsapiens-jaspar2022-ZNF135-MA1587.1"
mf = MotifDb[mf_names]

pfm1 = mf@listData[[mf_names[1]]] # mf[[1]]
ggseqlogo(pfm1, method='p') # great

pfm2 = mf@listData[[mf_names[2]]] # mf[[2]]
ggseqlogo(pfm2, method='p') # great

pfm3 = mf@listData[[mf_names[3]]] # mf[[3]]
ggseqlogo(pfm3, method='p') # great

#### 3) Load motifbreakR results to identify SNV location relative to the TF motifs
lnames = load("STRONGmotifbreakRresults_jaspar2022hsmotifs_Proband_SNVsWithin500bpOfDownregDistATAC.rda") # "results"        "df_results"     "strong_results"
Proband_StrongRes = strong_results
df = left_join(df, Proband_StrongRes)
df$motifPos = sapply(df$motifPos, function(x) paste(x, collapse = ","))
df$SNV_relPos = sapply(df$motifPos, function(x) 1-as.integer(unlist(str_split(x, ","))[1]))

# ----
### How to annotate TF binding sequence logo?
# Refer to https://omarwagih.github.io/ggseqlogo/
# data(ggseqlogo_sample)
# ggplot() + 
#   annotate('rect', xmin = 0.5, xmax = 3.5, ymin = -0.05, ymax = 1.9, alpha = .1, col='black', fill='yellow') +
#   geom_logo(seqs_dna$MA0001.1, stack_width = 0.90) + 
#   annotate('segment', x = 4, xend=8, y=1.2, yend=1.2, size=2) + 
#   annotate('text', x=6, y=1.3, label='Text annotation') + 
#   theme_logo()
# 
# ggplot() + 
#   annotate('rect', xmin = df$SNV_relPos[1] - 0.5, xmax = df$SNV_relPos[1] + 0.5, ymin = -0.05, ymax = 2, alpha = .1, col='black', fill='red') +
#   geom_logo(pfm1, stack_width = 0.90) + 
#   #geom_logo(pfm1, method='p', stack_width = 0.90) + 
#   #annotate('segment', x = 4, xend=8, y=1.2, yend=1.2, size=2) + 
#   annotate('text', x = df$SNV_relPos[1], y=2.5, label = paste0(df$SNV_id[1])) + 
#   theme_logo() +
#   xlab(paste(df$TF[1], "binding motif", df$TFmotif_id[1]))
# #xlab(gsub("Hsapiens-jaspar2022-","", mf_names[1]))
# Great
# ----

#### 4) Plot for SNV that disrupts TF binding within 500bp of ADSL/MAP2K2 distal ATAC peak
for (i in 1:nrow(df)) {
  RelPos = df$SNV_relPos[i]
  PFM = get(paste0("pfm", i))
  MotifRef = names(which.max(PFM[,RelPos]))
  
  if (MotifRef != df$REF[i]) {
    rownames(PFM) = c("T","G","C","A") # originally "A" "C" "G" "T"
  }
  
  p = ggplot() + 
    annotate('rect', xmin = RelPos - 0.5, xmax = RelPos + 0.5, ymin = -0.05, ymax = 2, alpha = .1, col='black', fill='red') +
    geom_logo(PFM, stack_width = 0.90) + 
    #geom_logo(pfm1, method='p', stack_width = 0.90) + 
    #annotate('segment', x = 4, xend=8, y=1.2, yend=1.2, size=2) + 
    annotate('text', x = RelPos, y=2.2, label = paste0(df$SNV_id[i]), size = 6) + 
    theme_logo() +
    xlab(paste(df$TF[i], "binding motif", df$TFmotif_id[i]))
  #xlab(gsub("Hsapiens-jaspar2022-","", mf_names[i]))
  
  assign(paste0("p", i), p)
}

pdf("13_6_exampleSNV_disrupt_TFmotif_within500bpofADSLorMAP2K2distaldownregATAC.pdf", height = 6, width = 10)
plot_grid(p1, p2, p3, align = "h") # looks great
dev.off()

pdf("13_6_exampleSNV_disrupt_TFmotif_within500bpofADSLorMAP2K2distaldownregATAC_1col.pdf", height = 10, width = 5)
plot_grid(p1, p2, p3, align = "h", ncol = 1)
dev.off()

save(df, pfm1, pfm2, pfm3, mf_names, p1, p2, p3, file = "13_6_exampleSNV_disrupt_TFmotif_within500bpofADSLorMAP2K2distaldownregATAC.rda")

## Li Ma suggests:
# Reverse the binding motif if SNV on the reverse strand. (edited in the loop that creates the plots)
# Add ATAC peak down-reg plot on the side? 
# Experimental validation?

#### 5) Check whether the two SNVs are tested in Matilde/Alex's MPRA

### A. Overlap Matilde's MPRA sigVars
rare = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/33_Matilde_ASD_rare_variants/SSC_ASD_proband_putative_functional_variants_MPRA_forJenny_MC_2022-09-01.txt", header = T)

df$variant.position = paste0(df$seqnames, ":", df$start)
df$variant.position %in% rare$variant.position # F

### B. Overlap Alex's MPRA sigVars
common = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/34_Alex_ASD_common_variants/mpra.phNP2merge.common.exp.SigVar.p=0_05.20211006.ann2_HighlyExpressedVariants.txt", header = T)
common$variant.position = paste0(common$chr,":",common$pos)

df$variant.position %in% common$variant.position # F

common = read.table("~/Documents/Documents/Geschwind_lab/LAB/Database_download/34_Alex_ASD_common_variants/mpra.phNP2merge.common.frVar.20220414_NoFilter.txt", header = T)
common$variant.position = sapply(common$var, function(x) unlist(str_split(x, "_"))[1])

df$variant.position %in% common$variant.position # F

### Conclusion: The SNVs of ADSL and MAP2K2 distal ATAC peaks are not sigVars in Matilde/Alex's MPRA.




