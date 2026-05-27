# Modified from 20_1_HighConfidenceLoops_EachCellType_VaryingPercentSamples.R

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)

## 1) Get gene promoter loops
# setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/20_HighConfidenceLoops_varyingPercentSamples_diffPromoterLoops/")

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/q001_AllFithic2SettingsSame")

load("Bulk_consensus_Autosome_loops_inAtLeast8outof36samples.rda") # consensus_Bulk
load("NeuNp_consensus_Autosome_loops_inAtLeast9outof15samples.rda") # consensus_NeuNp
load("NeuNn_consensus_Autosome_loops_inAtLeast7outof15samples.rda") # consensus_NeuNn
load("UpsetTbl_8Bulk7NeuNn9NeuNp_LoopFilter_20kbto5MbRmSegdup.rda") # Upset_tbl

consensus_Bulk = consensus_Bulk[consensus_Bulk$loop %in% Upset_tbl$loop,] # 71957 loops
consensus_NeuNp = consensus_NeuNp[consensus_NeuNp$loop %in% Upset_tbl$loop,] # 79504 loops
consensus_NeuNn = consensus_NeuNn[consensus_NeuNn$loop %in% Upset_tbl$loop,] # 80876 loops

# Limit to genes expressed in Jill's dataset
load("../txs_df.rda")
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

txs_df_BrExp = txs_df[txs_df$gene_id %in% Jill_DEG$ensembl_gene_id,]
promoter_bins = unique(txs_df_BrExp$chr_promoterbin) # 28292 promoter bins

getPL <- function(Loop_df) {
  Loop_df$anchor1 = paste0(Loop_df$V1, "_", Loop_df$V2)
  Loop_df$anchor2 = paste0(Loop_df$V1, "_", Loop_df$V3)
  PL_df = Loop_df[which(Loop_df$anchor1 %in% promoter_bins | Loop_df$anchor2 %in% promoter_bins),]
  return(PL_df)
}

PL_Bulk = getPL(consensus_Bulk) # 34657 consensus promoter loops
PL_NeuNp = getPL(consensus_NeuNp) # 38142 consensus promoter loops
PL_NeuNn = getPL(consensus_NeuNn) # 46180 consensus promoter loops

save(PL_Bulk, PL_NeuNp, PL_NeuNn, file = "Consensus_promoter_loops_8Bulk7NeuNn9NeuNp.rda")
write.table(PL_Bulk[,c(1,2,1,3)], file = "Consensus_promoter_loops_Bulk_location_8Bulk7NeuNn9NeuNp.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(PL_NeuNp[,c(1,2,1,3)], file = "Consensus_promoter_loops_NeuNp_location_8Bulk7NeuNn9NeuNp.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(PL_NeuNn[,c(1,2,1,3)], file = "Consensus_promoter_loops_NeuNn_location_8Bulk7NeuNn9NeuNp.txt", quote = F, row.names = F, col.names = F, sep = "\t")
# Transfer the three Consensus_promoter_loops_xxx.txt files to /u/project/geschwind/jennybea/ASD_project_2ndbatch/HiC/FitHiC/python_output_ASDL1subsetAndL2_Dup15q/consensusPL_AllFihic2SettingsSame/

## 2) Get loop by sample CC data frames
# See whether Bulk promoter loops all have CC in previous runs
rm(list = ls())
lnames = load("Consensus_promoter_loops_8Bulk7NeuNn9NeuNp.rda")
lnames = load("../../11_Loops_ASDL1subsetAndL2_Dup15q_sepCelltype/00_Dataset_Complete_BulkSamples.rda") # datCC, datMeta, datSeq, datSeq_numeric
rownames(PL_Bulk) = paste0(PL_Bulk$V1, "_", PL_Bulk$V2, "_", PL_Bulk$V1, "_", PL_Bulk$V3)
all(rownames(PL_Bulk) %in% rownames(datCC)) # True. Hooray, no need to re-get Bulk CCs

datCC$loop_alt = rownames(datCC)
PL_Bulk_location = PL_Bulk[,c(1,2,1,3)]; colnames(PL_Bulk_location) = c("V1", "V2", "V3", "V4")
PL_Bulk_location$loop_alt = rownames(PL_Bulk)
PL_Bulk_location = left_join(PL_Bulk_location, datCC)
idx = which(colnames(PL_Bulk_location) %in% colnames(PL_Bulk))
CCdf = PL_Bulk_location[,c(1:4, idx[4:length(idx)])]
save(CCdf, file = "../PL_by_sample_datasets/Consensus_promoter_loops_Bulk_8Bulk7NeuNn9NeuNp_LoopBySample_CC.rda")

# For Dup15q/CTL NeuN samples:
# Run 22_9_2_ContactCount_FilteredConsensusPromoterLoops_forEachSample_Dup15q_fithic2SettingsSameAsASD.sh
# Run 22_9_3_LoopContactCount_and_qval_by_Sample_Dup15q_fithic2SettingsSameAsASD.R to Combine CC of each NeuN sample on each chr into a dataframe

## 3) Calculate logCPM
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/PL_by_sample_datasets")

load("Consensus_promoter_loops_Bulk_8Bulk7NeuNn9NeuNp_LoopBySample_CC.rda")
Bulk_CC = CCdf
load("Consensus_promoter_loops_NeuNp_8Bulk7NeuNn9NeuNp_LoopBySample_CC.rda")
NeuNp_CC = CCdf
load("Consensus_promoter_loops_NeuNn_8Bulk7NeuNn9NeuNp_LoopBySample_CC.rda")
NeuNn_CC = CCdf
rm(CCdf)

# Remove Outlier and TechRep removal
#copied from 17_3_DiffLoopAnalysis_LinearModelBuild.R
Sample_rm = c("B4498", "B7575","B5242_Rep", "B5342_Rep", "B7436_NeuNn_Rep", "B5163_NeuNp")

rmOutlierTechrep <- function(df) {
  colnames(df)[which(colnames(df) == "B5242A")] = "B5242"
  colnames(df)[which(colnames(df) == "B5242B")] = "B5242_Rep"
  colnames(df)[which(colnames(df) == "B5342A")] = "B5342"
  colnames(df)[which(colnames(df) == "B5342B")] = "B5342_Rep"
  colnames(df)[which(colnames(df) == "B2987")] = "B7575"
  #colnames(df)[match(c("B5242A", "B5242B", "B5342A", "B5342B", "B2987"), colnames(df))] = c("B5242", "B5242_Rep", "B5342", "B5342_Rep", "B7575")
  
  idx_rm = which(colnames(df) %in% Sample_rm) # 6
  if (length(idx_rm) > 0) {df = df[,-idx_rm]}
  return(df)
}

Bulk_CC = rmOutlierTechrep(Bulk_CC)
NeuNp_CC = rmOutlierTechrep(NeuNp_CC)
NeuNn_CC = rmOutlierTechrep(NeuNn_CC)

save(Bulk_CC, NeuNp_CC, NeuNn_CC, file = "Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_CC_rmOutlierTechRep.rda")

# Get logCPM
Bulk_logCPM = as.data.frame(apply(Bulk_CC[,5:ncol(Bulk_CC)], 2, function(x) log2((x+0.5)/(sum(x+0.5)/1e6))))
NeuNp_logCPM = as.data.frame(apply(NeuNp_CC[,5:ncol(NeuNp_CC)], 2, function(x) log2((x+0.5)/(sum(x+0.5)/1e6))))
NeuNn_logCPM = as.data.frame(apply(NeuNn_CC[,5:ncol(NeuNn_CC)], 2, function(x) log2((x+0.5)/(sum(x+0.5)/1e6))))
rownames(Bulk_logCPM) = paste0(Bulk_CC$V1, "_", Bulk_CC$V2, "_", Bulk_CC$V4)
rownames(NeuNp_logCPM) = paste0(NeuNp_CC$V1, "_", NeuNp_CC$V2, "_", NeuNp_CC$V4)
rownames(NeuNn_logCPM) = paste0(NeuNn_CC$V1, "_", NeuNn_CC$V2, "_", NeuNn_CC$V4)

save(Bulk_logCPM, NeuNp_logCPM, NeuNn_logCPM, file = "Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_rmOutlierTechRep.rda")

pdf(paste0("Histogram_logCPM_ConsensusPLs_8Bulk7NeuNn9NeuNp.pdf"), height = 4, width = 8)
par(mfrow = c(1,3))
hist(as.matrix(Bulk_logCPM), breaks = 100, main = "Bulk", ylim = c(0, 7e4), xlim = c(-2,8), yaxp = c(0, 7e4, 7), xlab = "logCPM") # Histogram of logCPM - Bulk
hist(as.matrix(NeuNp_logCPM), breaks = 100, main = "NeuNp", ylim = c(0, 4e4), xlim = c(-2,8), yaxp = c(0, 4e4, 4), xlab = "logCPM")
hist(as.matrix(NeuNn_logCPM), breaks = 100, main = "NeuNn", ylim = c(0, 4e4), xlim = c(-2,8), yaxp = c(0, 4e4, 4), xlab = "logCPM")
dev.off()
# Observation:
# Bulk is about twice the height of NeuNp or NeuNn. Makes sense, as Bulk has twice the number of samples

png("QQplot_Bulk_logCPM_ConsensusPLs_8Bulk7NeuNn9NeuNp.png", res = 150)
qqnorm(as.matrix(Bulk_logCPM))
qqline(as.matrix(Bulk_logCPM), col = "red", lwd = 2)
# The lower part is off diagonal.
dev.off()

png("QQplot_NeuNp_logCPM_ConsensusPLs_8Bulk7NeuNn9NeuNp.png", res = 150)
qqnorm(as.matrix(NeuNp_logCPM))
qqline(as.matrix(NeuNp_logCPM), col = "red", lwd = 2)
# The lower part is off diagonal.
dev.off()

png("QQplot_NeuNn_logCPM_ConsensusPLs_8Bulk7NeuNn9NeuNp.png", res = 150)
qqnorm(as.matrix(NeuNn_logCPM))
qqline(as.matrix(NeuNn_logCPM), col = "red", lwd = 2)
# The lower part is off diagonal.
dev.off()
