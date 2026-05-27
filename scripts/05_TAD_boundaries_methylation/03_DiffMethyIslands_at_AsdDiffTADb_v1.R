### Goal: 1) To ask whether CTCF binding sites are differentially methylated at differential/weakened TAD boundaries. 2) To ask whether differentially methylated regions (DMRs) are enriched at differential/weakened TAD boundaries.

### Datasets:
# Differential TAD boundaries in ASD samples were determined by 6_5_differential_ISnonsegdup_analysis_v4_afterGRC_TADbOnly_includeTransCisInLm.R
# DMRs: re-analyzed Wong et al 2019 using 02_1_redo_DMR_analysis.R

#### ----- Procedures --------
### (1) Identify enriched TF binding sites in up- or down-reg DMRs - Expect CTCF in hypermethylated DMRs, because we see reduced CTCF binding and hypermethylation antagonize CTCF binding (Wang 2012 Genome Res).
### (2) Identify CTCF-binding sites by motif scan (fimo p<1e5 used by Wang 2012 Genome Res) at CTCF-bound ATAC peaks.
### (3) Hypothesis: Significant increase in %DMRup at weakened TAD boundaries
### (4) Hypothesis: Significant increase in %DMRup CTCF sites at weakened TAD boundaries
### (5) Hypothesis: CpG islands with CTCF binding sites at weakened TAD boundaries are significantly more likely to be hyper-methylated
### (6) Example weakened TADb with hyper-methylated CpG island
# ---------------------------

rm(list = ls())
options(stringsAsFactors = F, scipen = 999)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/03_use_reanalyzed_DiffMethyIslands")

library(tidyverse)    
library(ggplot2)    

### (1) Identify enriched TF binding sites in up- or down-reg DMRs

## Load DMRs
lnames = load("../02_redo_DMG_analysis/08_DMR_fdr005.rda") # DMR, HM450k_final, Probe_Island_grouped
all(rownames(DMR) == Probe_Island_grouped$UCSC_CpG_Islands_Name) # T

length(which(DMR$magnitude < 0 & DMR$fdr < 0.05)) # 2261
length(which(DMR$magnitude > 0 & DMR$fdr < 0.05)) # 6638 - hypermetylation!

# Refer to 9_3_HOMERofDiffATACupORdown_Jaspar2022.R
## Create bed files of the up- and down-regulated DMRs for HOMER

DMRs_up = DMR[DMR$magnitude > 0 & DMR$fdr < 0.05,] # 6646
DMRs_down = DMR[DMR$magnitude < 0 & DMR$fdr < 0.05,] # 2269
# Most DMRs show hyper-methylation. If they locate at CTCF binding sites, this would be consistent with reduced CTCF binding and weakened TAD boundaries.
DMRs_up_location = as.data.frame(t(as.data.frame(str_split(rownames(DMRs_up), ":|-"))))
DMRs_down_location = as.data.frame(t(as.data.frame(str_split(rownames(DMRs_down), ":|-"))))
DMRs_up = cbind(DMRs_up_location, DMRs_up)
DMRs_down = cbind(DMRs_down_location, DMRs_down)

write.table(DMRs_up, file = "01_DMRup.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(DMRs_down, file = "01_DMRdown.bed", quote = F, sep = "\t", col.names = F, row.names = F)
save(DMRs_up, DMRs_down, file = "01_DMRs.rda")

## Run HOMER: 03_1_HOMER_MotifEnrichment_in_DMRsUpDown_Jaspar2022.sh

### TF binding motif enrichment plot
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/03_use_reanalyzed_DiffMethyIslands/")
(motif_files = list.files(pattern = "*_knownResults.txt")) # 2

for (motif_file in motif_files) {
  new_name = substr(motif_file,1,nchar(motif_file) - 4)
  (new_name = gsub("_Jaspar2022", "", new_name))
  (new_name = gsub("01_", "", new_name))
  
  new_motif = read.table(paste0(motif_file)) # "HOMER_TFenrichment_results/", 
  colnames(new_motif) = c("Consensus", "P_value", "logP_value", "q_value_Benjamini", 
                          "Number_of_Target_Sequences_with_Motif","Percent_of_Target_Sequences_with_Motif",
                          "Number_of_Background_Sequences_with_Motif","Percent_of_Background_Sequences_with_Motif")
  new_motif$Motif_name = rownames(new_motif)
  assign(new_name, new_motif)
}

range(DMRdown_knownResults$q_value_Benjamini) # 0.58 1
range(DMRup_knownResults$q_value_Benjamini) # 0 1. CTCF ranks the thrid!
DMRdown_knownResults$P_value[startsWith(DMRdown_knownResults$Motif_name, "CTCF")] # 0.1, 1, 1
DMRup_knownResults$P_value[startsWith(DMRup_knownResults$Motif_name, "CTCF")] # MA0139.1 got 1e-8; the other two got 1.

## Generate a matrix comparing the q-value of each motif.
ref = DMRup_knownResults
DMRdown_knownResults = DMRdown_knownResults[match(rownames(ref), rownames(DMRdown_knownResults)),]
#upATAC_threshold2_knownResults = upATAC_threshold2_knownResults[match(rownames(ref), rownames(upATAC_threshold2_knownResults)),]

Motif_comparison_main = 
  data_frame(Motif_name = ref$Motif_name,
             Consensus = ref$Consensus,
             
             pval_up = DMRup_knownResults$P_value,
             pval_down = DMRdown_knownResults$P_value,
             
             qval_up = DMRup_knownResults$q_value_Benjamini,
             qval_down = DMRdown_knownResults$q_value_Benjamini,
  )

idx = which(apply(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "qval"))], 1, function(x) any(x < 0.05)))
Motif_comparison_main = 
  Motif_comparison_main[idx,] # only FDR<0.05 are shown

Motif_comparison_main$minus_log_qval_up = -log10(Motif_comparison_main$qval_up)
Motif_comparison_main$minus_log_qval_down = -log10(Motif_comparison_main$qval_down)

Motif_comparison_main$Motif_abbr = sapply(Motif_comparison_main$Motif_name, function(x) substr(x, 1, which(strsplit(x, "")[[1]] == "/")[1] -1))
any(duplicated(Motif_comparison_main$Motif_abbr)) # F, all motifs are unique

library(reshape2)
Motif_comparison_main_melt = melt(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "minus_log_qval") | colnames(Motif_comparison_main) == "Motif_abbr")])
colnames(Motif_comparison_main_melt)[3] = c("significance")
Motif_comparison_main_melt$variable = as.character(Motif_comparison_main_melt$variable)
Motif_comparison_main_melt$variable = gsub("minus_log_qval_", "", Motif_comparison_main_melt$variable)

Motif_comparison_main_melt$Motif_abbr =
  factor(Motif_comparison_main_melt$Motif_abbr,
         levels = Motif_comparison_main$Motif_abbr) 
Motif_comparison_main_melt$variable = ifelse(Motif_comparison_main_melt$variable == "down", "hypo-methylated", "hyper-methylated")
Motif_comparison_main_melt$variable = factor(Motif_comparison_main_melt$variable, levels = c("hyper-methylated", "hypo-methylated"))
#Motif_comparison_main_melt = Motif_comparison_main_melt[!is.na(Motif_comparison_main_melt$Motif_abbr),]

save(Motif_comparison_main, Motif_comparison_main_melt, file = "01_Motif_Enrichment_in_DMRupORdown_HOMERq005_Jaspar2022.rda")

# plot
rm(list = ls())
lnames = load("01_Motif_Enrichment_in_DMRupORdown_HOMERq005_Jaspar2022.rda")

TextSize = 5
LegendSize = 0.5

Motif_plot = Motif_comparison_main_melt %>%
  ggplot(aes(x = Motif_abbr, y = significance, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Enrichment of TFBS in ASD differentially methylated CpG islands") + # (Only TF motifs with FDR<0.05 are shown)
  geom_hline(aes(yintercept=-log10(0.05), linetype="FDR = 0.05"), color = "black", size=0.3) +
  ylab("-log10(FDR) of TF enrichment") +
  xlab("Transcription factors") +
  labs(fill = "DMRs") +
  theme_bw() +
  theme(text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize),
        axis.text.x = element_text(size = TextSize, angle = 45, hjust = 1),
        legend.text = element_text(size = TextSize), # required to have legends at 5 pt
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        legend.margin = margin(0,0,-0.3,0,unit = "lines"),
        legend.key.size = unit(LegendSize, 'lines'),
        legend.position = "top"
  ) +
  scale_linetype_manual(name = "Threshold", values = c("FDR = 0.05" = 2)) + #, "FDR = 0.2" = 3
  scale_fill_manual(values = c("hypo-methylated" = "deepskyblue", "hyper-methylated" = "#F8766D")) # purple #C77CFF

pdf("01_Motif_Enrichment_in_DMRupORdown_HOMERq005_Jaspar2022.pdf", width = 3.5, height = 2.5)
Motif_plot
dev.off()

# Conclusion: CTCF is enriched at hyper-methylated CpG islands.

### (2) Identify CTCF-binding sites and Overlap with ChIP-seq from brain samples 
## Done in 01_DiffMethylProbes_at_AsdDiffTADb_v1.R
# Conclusion: There is only a 40-60% overlap between the ChIP-seq peaks and CTCF BS within ATAC peaks. Could be a difference in brain region and age. But it is not ideal.

### (3) Hypothesis: Significant increase in %DMRup at weakened TAD boundaries
# We have n weakened TAD boundaries (FDR<0.05), and x% DMP up-reg out of all the probed CG pairs in this region. We will draw a background distribution by randomly select n TAD boundaries and inquire the % DMP up-reg in those n TADbs. Finally assign a p-val.
# May limit to CTCF within 50bp of probed CG pairs later
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/03_use_reanalyzed_DiffMethyIslands/")

library(Repitools)
library(GenomicRanges)

## Load Differential TAD boundaries
ASD_TAD_dir = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/6_TADs/IS_GENOVA/differential_IS/v4_afterGRC_TADbOnly_includeTransCisInLm/"

lnames = load(paste0(ASD_TAD_dir, "p_magnitude_Bulk_lmeDiagnosisAgeSexBatch2.TemporalABNNICHDvalid.interactioncis.longRange_includingTechRepsAndOutliers.rda")) # p, magnitude, diffIS
length(which(diffIS$FDR < 0.05)) # 1051 diff TADb
(N = length(which(diffIS$FDR < 0.05 & diffIS$Magnitude > 0))) # 754 weakened TADb 
length(which(diffIS$FDR < 0.05 & diffIS$Magnitude < 0)) # 297 strengthened TADb
nrow(diffIS) # Total 11167 TAD boundaries

## Load DMRs (re-analyzed)
#lnames = load("01_DMRs.rda")
lnames = load("../02_redo_DMG_analysis/08_DMR_fdr005.rda") # DMR, HM450k_final, Probe_Island_grouped
# lnames = load("../02_redo_DMG_analysis/06_final_lm.rda") # expression_lm, mod_mat, Covariates, datMeth_CTX_island, HM450k_final, Probe_Island_grouped

## % DMR up-reg at weakened TADb
TADb_location = as.data.frame(str_split_fixed(rownames(diffIS), "_", 3))
TADb_location$V2 = as.integer(TADb_location$V2); TADb_location$V3 = as.integer(TADb_location$V3)
colnames(TADb_location) = c("chr", "start", "end")

weakened_TADb = TADb_location[which(diffIS$FDR < 0.05 & diffIS$Magnitude > 0), ]
weakened_TADb_gr = annoDF2GR(weakened_TADb)

all(Probe_Island_grouped$UCSC_CpG_Islands_Name == rownames(DMR)) # T
idx_up = which(DMR$fdr < 0.05 & DMR$magnitude > 0) # 6638
idx_down = which(DMR$fdr < 0.05 & DMR$magnitude < 0) # 2261
DMRs_up = DMR[idx_up,]
DMRs_down = DMR[idx_down,]
DMRs_up = cbind(DMRs_up, Probe_Island_grouped[idx_up,])
DMRs_down = cbind(DMRs_down, Probe_Island_grouped[idx_down,])
DMRup_gr = annoDF2GR(DMRs_up)
DMRdown_gr = annoDF2GR(DMRs_down)

Islands_gr = annoDF2GR(Probe_Island_grouped)

save(Islands_gr, DMRup_gr, DMRdown_gr, weakened_TADb_gr, Probe_Island_grouped, DMRs_up, DMRs_down, weakened_TADb, TADb_location, file = "03_Probe_weakenedTAD_gr.rda")

rm(list = ls())
lnames = load("03_Probe_weakenedTAD_gr.rda")

df_DMRup_in_weakenedTADb = as.data.frame(findOverlaps(weakened_TADb_gr, DMRup_gr))
df_DMRup_in_weakenedTADb$Island = DMRs_up$UCSC_CpG_Islands_Name[df_DMRup_in_weakenedTADb$subjectHits]
DMRup_in_weakenedTADb = unique(df_DMRup_in_weakenedTADb$Island) # 216 CpG Islands hyper-methylated in weakened TADb

df_Islands_in_weakenedTADb = as.data.frame(findOverlaps(weakened_TADb_gr, Islands_gr))
df_Islands_in_weakenedTADb$Island = Probe_Island_grouped$UCSC_CpG_Islands_Name[df_Islands_in_weakenedTADb$subjectHits]
Islands_in_weakenedTADb = unique(df_Islands_in_weakenedTADb$Island) # 614 CpG Islands in total

p_DMRup_obs = length(DMRup_in_weakenedTADb)/length(Islands_in_weakenedTADb) # 35.2% CpG islands up-reg in weakened TADb. 

## Draw background distribution of %DMP up-reg at TADb

p_DMRup_null = c()
set.seed(1002)
N = 754
B = 500
for (i in 1:B) {
  idx = sample(1:nrow(TADb_location), N)
  TAD_idx = TADb_location[idx,]
  TADb_idx_gr = annoDF2GR(TAD_idx)
  
  df_DMRup_in_TADbidx = as.data.frame(findOverlaps(TADb_idx_gr, DMRup_gr))
  df_DMRup_in_TADbidx$Island = DMRs_up$UCSC_CpG_Islands_Name[df_DMRup_in_TADbidx$subjectHits]
  Island_DMRup_in_TADbidx = unique(df_DMRup_in_TADbidx$Island)
  
  df_Island_in_TADbidx = as.data.frame(findOverlaps(TADb_idx_gr, Islands_gr))
  df_Island_in_TADbidx$Island = Probe_Island_grouped$UCSC_CpG_Islands_Name[df_Island_in_TADbidx$subjectHits]
  Islands_All_in_TADbidx = unique(df_Island_in_TADbidx$Island)
  
  p_DMRup_idx = length(Island_DMRup_in_TADbidx)/length(Islands_All_in_TADbidx)
  
  p_DMRup_null = c(p_DMRup_null, p_DMRup_idx)
}

## Draw histogram and assign p-val
(DMRup_pval = length(which(p_DMRup_null >= p_DMRup_obs))/B) # 0.02! It is significant!
pdf("03_Histogram_percent_DMRup_in_weakened_TAD_boundaries.pdf", height = 4, width = 6)
hist(p_DMRup_null, breaks = 50, main = paste0("p = ", DMRup_pval),
     xlab = "% hyper-methylated CpG Islands at weakened TAD boundaries",
     xlim = c(0.2,0.4))
abline(v = p_DMRup_obs, col = "red", lwd = 2)
dev.off()

save(p_DMRup_null, p_DMRup_obs, DMRup_pval, file = "03_Histogram_percent_DMRup_in_weakened_TAD_boundaries.rda")

### (5) Hypothesis: CpG islands with CTCF binding sites at weakened TAD boundaries are significantly more likely to be hyper-methylated
# We have n weakened TAD boundaries (FDR<0.05), and x% CpG islands with CTCF_BS show hypermethylation. We will draw a background distribution by randomly select n TAD boundaries and inquire the % CpG islands with CTCF_BS that show hypermethylation n TADbs. Finally assign a p-val.

## Load DMRs, CTCF-bound ATAC peaks, weakened TADs, 
rm(list = ls())
lnames = load("../02_redo_DMG_analysis/08_DMR_fdr005.rda") # DMR, HM450k_final Probe_Island_grouped 
lnames = load("03_Probe_weakenedTAD_gr.rda") # Islands_gr, DMRup_gr, DMRdown_gr, weakened_TADb_gr, Probe_Island_grouped, DMRs_up, DMRs_down, weakened_TADb, TADb_location
lnames = load("04_CTCFBS_inIslands.rda") # CTCF_inIslands_gr, CTCF_inIslands, CTCF_gr, CTCF_inATAC

## Restrict to probed CpG islands with CTCF_BS
Islands_wCTCF = as.data.frame(findOverlaps(Islands_gr, CTCF_gr))
Islands_wCTCF = Probe_Island_grouped[unique(Islands_wCTCF$queryHits),] # 13696 out of 24717 CpG Islands have CTCF binding
Islands_wCTCF_gr = annoDF2GR(Islands_wCTCF)

## % CpG Islands with CTCF_BS showing hyper-methylation at weakened TADb
df_upDMRwCTCF = DMRs_up[DMRs_up$UCSC_CpG_Islands_Name %in% Islands_wCTCF$UCSC_CpG_Islands_Name,] # 5177 DMR up-reg with CTCF_BS
gr_upDMRwCTCF = annoDF2GR(df_upDMRwCTCF)

df_upDMRwCTCF_in_weakenedTADb = as.data.frame(findOverlaps(gr_upDMRwCTCF, weakened_TADb_gr))
df_upDMRwCTCF_in_weakenedTADb$Island = df_upDMRwCTCF$UCSC_CpG_Islands_Name[df_upDMRwCTCF_in_weakenedTADb$queryHits]
upDMRwCTCF = unique(df_upDMRwCTCF_in_weakenedTADb$Island) # 173 CTCF_BS containing up-regulated DMRs overlap with weakened TAD boundaries!

df_IslandwCTCF_in_weakenedTADb = as.data.frame(findOverlaps(Islands_wCTCF_gr, weakened_TADb_gr))
df_IslandwCTCF_in_weakenedTADb$Island = Islands_wCTCF$UCSC_CpG_Islands_Name[df_IslandwCTCF_in_weakenedTADb$queryHits]
IslandwCTCF = unique(df_IslandwCTCF_in_weakenedTADb$Island) # 396 CTCF_BS containing CpG islands overlap with weakened TAD boundaries!

p_upDMRwCTCF_obs = length(upDMRwCTCF)/length(IslandwCTCF) # 43.7% CpG islands with CTCF binding sites in weakened TADb were hyper-methylated. 

## Draw background distribution of %CTCF_BS with up-reg DMR at TADb

p_upDMRwCTCF_null = c()
set.seed(1002)
N = 754
B = 500
for (i in 1:B) {
  idx = sample(1:nrow(TADb_location), N)
  TAD_idx = TADb_location[idx,]
  TADb_idx_gr = annoDF2GR(TAD_idx)
  
  df_upDMRwCTCF_in_TADbidx = as.data.frame(findOverlaps(gr_upDMRwCTCF, TADb_idx_gr))
  df_upDMRwCTCF_in_TADbidx$Island = df_upDMRwCTCF$UCSC_CpG_Islands_Name[df_upDMRwCTCF_in_TADbidx$queryHits]
  upDMRwCTCF = unique(df_upDMRwCTCF_in_TADbidx$Island)
  
  df_IslandwCTCF_in_TADbidx = as.data.frame(findOverlaps(Islands_wCTCF_gr, TADb_idx_gr))
  df_IslandwCTCF_in_TADbidx$Island = Islands_wCTCF$UCSC_CpG_Islands_Name[df_IslandwCTCF_in_TADbidx$queryHits]
  IslandwCTCF = unique(df_IslandwCTCF_in_TADbidx$Island)
  
  p_upDMRwCTCF_idx = length(upDMRwCTCF)/length(IslandwCTCF)
  p_upDMRwCTCF_null = c(p_upDMRwCTCF_null, p_upDMRwCTCF_idx)
}

## Draw histogram and assign p-val
(upDMRwCTCF_pval = length(which(p_upDMRwCTCF_null >= p_upDMRwCTCF_obs))/B) # 0.098!

pdf("05_Histogram_percent_upDMRwithCTCFBS_in_weakened_TAD_boundaries.pdf", height = 4, width = 6)
hist(p_upDMRwCTCF_null, breaks = 50, main = paste0("p = ", upDMRwCTCF_pval),
     xlab = "% hyper-methylated CTCF BS-containing CpG Islands\nat weakened TAD boundaries",
     xlim = c(0.25,0.55))
abline(v = p_upDMRwCTCF_obs, col = "red", lwd = 2)
dev.off()

## Conclusion: CTCF BS-containing CpG islands at weakened TAD boundaries are more likely to be hyper-methylated (p < 0.1).

### (6) Example weakened TADb with hyper-methylated CpG island
# Find a sample region, plot CpG island methylation profile, Hi-C heatmap showing TAD boundary weakening. Ideally BDNF gene. 

## Is BDNF gene located at weakened TAD boundaries?
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")
colnames(Jill_DEG)[c(4,7:8)] = c("chr", "start", "end")
Jill_DEG$chr = paste0("chr", Jill_DEG$chr)
DEG_gr = annoDF2GR(Jill_DEG[,-5])

rownames(weakened_TADb) = paste0(weakened_TADb$chr, "_", weakened_TADb$start, "_", weakened_TADb$end)
df_upDMRwCTCF_in_weakenedTADb$TADb = rownames(weakened_TADb)[df_upDMRwCTCF_in_weakenedTADb$subjectHits]
df_upDMRwCTCF_in_weakenedTADb = cbind(df_upDMRwCTCF_in_weakenedTADb, weakened_TADb[df_upDMRwCTCF_in_weakenedTADb$subjectHits,])

weakenedTADb_w_upDMRwCTCF = unique(df_upDMRwCTCF_in_weakenedTADb[,-c(1:3)]) # 157 TADb encompass hypermethylated CpG islands with CTCF binding sites.
weakenedTADb_w_upDMRwCTCF_gr = annoDF2GR(weakenedTADb_w_upDMRwCTCF)

Genes_at_weakenedTADb_w_upDMRwCTCF = as.data.frame(findOverlaps(weakenedTADb_w_upDMRwCTCF_gr, DEG_gr))
Genes_at_weakenedTADb_w_upDMRwCTCF = cbind(Genes_at_weakenedTADb_w_upDMRwCTCF, Jill_DEG[Genes_at_weakenedTADb_w_upDMRwCTCF$subjectHits,c(1:3, 10:13)])
Genes_at_weakenedTADb_w_upDMRwCTCF = cbind(Genes_at_weakenedTADb_w_upDMRwCTCF, weakenedTADb_w_upDMRwCTCF[Genes_at_weakenedTADb_w_upDMRwCTCF$queryHits,])
Genes_at_weakenedTADb_w_upDMRwCTCF = unique(Genes_at_weakenedTADb_w_upDMRwCTCF[,-c(1:2)])

# BDNF is not here. 

## Any other SFARI genes?
SFARI = read.csv("~/Documents/Documents/Geschwind_lab/LAB/Database_download/SFARI_human_genes/SFARI-Gene_genes_01-11-2022release_03-03-2022export.csv")
SFARI_genes = SFARI[SFARI$gene.score %in% c("S", "1", "2"),] # 554 genes

Genes_at_weakenedTADb_w_upDMRwCTCF$SFARI = SFARI_genes$gene.score[match(Genes_at_weakenedTADb_w_upDMRwCTCF$ensembl_gene_id, SFARI_genes$ensembl.id)]
table(Genes_at_weakenedTADb_w_upDMRwCTCF$SFARI) 

ASD_TAD_dir = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/6_TADs/IS_GENOVA/differential_IS/v4_afterGRC_TADbOnly_includeTransCisInLm/"
lnames = load(paste0(ASD_TAD_dir, "p_magnitude_Bulk_lmeDiagnosisAgeSexBatch2.TemporalABNNICHDvalid.interactioncis.longRange_includingTechRepsAndOutliers.rda")) # p, magnitude, diffIS
Genes_at_weakenedTADb_w_upDMRwCTCF$DiffIS_Magnitude = diffIS$Magnitude[match(Genes_at_weakenedTADb_w_upDMRwCTCF$TADb, rownames(diffIS))]
Genes_at_weakenedTADb_w_upDMRwCTCF$DiffIS_FDR = diffIS$FDR[match(Genes_at_weakenedTADb_w_upDMRwCTCF$TADb, rownames(diffIS))]

Sfari_genes_at_weakenedTADb_w_upDMRwCTCF = Genes_at_weakenedTADb_w_upDMRwCTCF[complete.cases(Genes_at_weakenedTADb_w_upDMRwCTCF),]
Sfari_genes_at_weakenedTADb_w_upDMRwCTCF$external_gene_name
#  [1] *"ZMIZ1"    "APBB1"    "DLG2"     *"C12orf57" "PRICKLE1" "INTS6"    *"CCNK"     *"USP7"  "QRICH1"   "UBN2"     "UBR5"    
Sfari_genes_at_weakenedTADb_w_upDMRwCTCF[Sfari_genes_at_weakenedTADb_w_upDMRwCTCF$WholeCortex_ASD_FDR < 0.05, c(2,4,5,12)]
#                            external_gene_name WholeCortex_ASD_logFC WholeCortex_ASD_FDR SFARI
# chr11_6440000_6480000.2                 APBB1               -0.1680            0.000158     2
# chr12_42840000_42880000.1            PRICKLE1               -0.2120            0.010400     2
# chr13_52000000_52040000                 INTS6                0.0690            0.006370     2
# chr16_9040000_9080000                    USP7               -0.0952            0.007920     S
# chr8_103240000_103280000.1               UBR5                0.0577            0.029300     2

save(Genes_at_weakenedTADb_w_upDMRwCTCF, weakenedTADb_w_upDMRwCTCF, Sfari_genes_at_weakenedTADb_w_upDMRwCTCF, file = "06_Genes_at_weakenedTADb_w_upDMRwCTCF_labelSFARI.rda")

### Plot 
library(plotgardener)
# a. ASD-CTL Hi-C heatmap
# b. CpG islands at the region (red/blue/grey as hyper/hypo/ns-methylated) (DMR below)
# c. CTCF_BS bound by CTCF in ATAC peaks (CTCF_inATAC below)
# d. genes
# e. Bulk ATAC track

rm(list = ls())

## Load DMRs, weakened TADb, CTCF BS
lnames = load("../02_redo_DMG_analysis/08_DMR_fdr005.rda") # DMR, HM450k_final Probe_Island_grouped 
lnames = load("03_Probe_weakenedTAD_gr.rda") # Islands_gr, DMRup_gr, DMRdown_gr, weakened_TADb_gr, Probe_Island_grouped, DMRs_up, DMRs_down, weakened_TADb, TADb_location
lnames = load("04_CTCFBS_inIslands.rda") # CTCF_inIslands_gr, CTCF_inIslands, CTCF_gr, CTCF_inATAC
lnames = load("06_Genes_at_weakenedTADb_w_upDMRwCTCF_labelSFARI.rda") # Genes_at_weakenedTADb_w_upDMRwCTCF, weakenedTADb_w_upDMRwCTCF, Sfari_genes_at_weakenedTADb_w_upDMRwCTCF

# Try SFARI genes first
Gene = "UBR5" # USP7, PRICKLE1, UBR5, C12orf57, ZMIZ1, DLG2
Gene_id = Sfari_genes_at_weakenedTADb_w_upDMRwCTCF$ensembl_gene_id[Sfari_genes_at_weakenedTADb_w_upDMRwCTCF$external_gene_name == Gene]
(TADb_weakened = rownames(Sfari_genes_at_weakenedTADb_w_upDMRwCTCF)[Sfari_genes_at_weakenedTADb_w_upDMRwCTCF$external_gene_name == Gene]) # chr16_9040000_9080000 for USP7
TADb_weakened = unlist(str_split(TADb_weakened, "\\."))[1]

## Identify DMRs
DMR$DMR = 
  ifelse(DMR$fdr < 0.05 & DMR$magnitude > 0, "hyper_methylated", 
         ifelse(DMR$fdr < 0.05 & DMR$magnitude < 0, "hypo_methylated", "no_signif_change"))
#DMR_location = Probe_Island_grouped[match(rownames(DMR), Probe_Island_grouped$UCSC_CpG_Islands_Name),3:5]
all(rownames(DMR) == Probe_Island_grouped$UCSC_CpG_Islands_Name) # T
DMR = cbind(Probe_Island_grouped[,3:5], DMR)

rm(list = setdiff(ls(), c("DMR", "CTCF_inATAC", "Gene", "Gene_id", "TADb_weakened")))

## For now, use ATAC bigwig of sample CQ56-2_unique.bw
CTL_ATAC_bw = readBigwig(file = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/6_logFC_cor_DAR_DEG/IGV/bw/CQ56-2_unique.bw") 
colnames(CTL_ATAC_bw)[1] = "chrom"

## Hi-C
ASD_GENOVA_dir = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/6_TADs/IS_GENOVA/6_9_aggregated_matrix_Batch1and2/W16_MinstrengthMinus1/"
ASD_TopDom_dir = "~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/6_TADs/TopDom/6_9_aggregated_matrix_Batch1and2/"

## Plotgardener 
# refer to 16_2_v2_caQTLBryois_plot_EffectOnATAC_vs_ASDgwasZ_ColorByPromDistUpDown.R and Rscripts_v2/Figure6d_plotgardener_NACC1.R

TextSize = 5
PointSize = 1
pre_name = "06_"
i=1; DOMAIN_COL = "red"; ws = 5
# UBR5:
HALF_window = 1.25e6; CUTOFF_CTL = 800
# Looks the best, but there is another TADb that has hyper-methylated CpG island with CTCF_BS that does not show differential IS.
# looks good, just too busy (too many CpG islands and ATAC peaks in the region)
# Conclusion: UBR5 looks the best.

candidate_clear = data.frame(
  chr = unlist(str_split(TADb_weakened, "_"))[1],
  TADb_start = as.integer(unlist(str_split(TADb_weakened, "_"))[2]),
  TADb_end = as.integer(unlist(str_split(TADb_weakened, "_"))[3]),
  external_gene_name = Gene,
  ensembl_gene_id = Gene_id
)

i=1
CHR = candidate_clear$chr[i]
View_center = (candidate_clear$TADb_start + candidate_clear$TADb_end)/2
START = round((View_center - HALF_window)/40e3) * 40e3  
END = round((View_center + HALF_window)/40e3) * 40e3
CTCF_HALFwidth = (END-START)/500 # 300bp for window size 600kb

# Bulk CTL Hi-C matrix
SampleName = "Batch1and2_CTL"
load(paste0(ASD_GENOVA_dir, SampleName, "_loadcontacts.rda"))
load(paste0(ASD_TopDom_dir, SampleName, "_TADcallsDomains_TopDom2GENOVA_windowsize",ws,".rda"))

ABS = SampleName_40kb$IDX
chr_range = range(ABS$V4[ABS$V1 == CHR])
ABS_subset = ABS %>%
  filter(V1 == CHR, V3 >= START, V2 <= END)

SampleName_40kb_MAT = SampleName_40kb$MAT
SampleName_40kb_chr = SampleName_40kb_MAT %>%
  filter(V1 >= chr_range[1], V2 <= chr_range[2], V1 != V2)
# Clean up Diagonal
SampleName_40kb_chr$V3[which(SampleName_40kb_chr$V1 == SampleName_40kb_chr$V2)] = 0

SampleName_40kb_subset = SampleName_40kb_chr %>%
  filter(V1 >= min(ABS_subset$V4), V2 <= max(ABS_subset$V4))
colnames(SampleName_40kb_subset)[3] = "counts"
# Set diagonal as max
Diagonal_subset = data.frame(V1 = seq(min(ABS_subset$V4), max(ABS_subset$V4), 1), V2 = seq(min(ABS_subset$V4), max(ABS_subset$V4), 1), counts = max(SampleName_40kb_subset$counts))
SampleName_40kb_subset = rbind(SampleName_40kb_subset, Diagonal_subset)
# Add genomic coordinates
SampleName_40kb_subset$V4 = ABS_subset$V2[match(SampleName_40kb_subset$V1, ABS_subset$V4)] # Note that I should use V2 instead of mid
SampleName_40kb_subset$V5 = ABS_subset$V2[match(SampleName_40kb_subset$V2, ABS_subset$V4)]
SampleName_40kb_subset = SampleName_40kb_subset[,c(4,5,3,1,2)]

Domains_relevant = Domains[Domains$chrom == CHR & Domains$start > START & Domains$end < END,]

tmp_list = list(MAT = SampleName_40kb_subset, Domains_relevant = Domains_relevant, SUM_CC = sum(SampleName_40kb_chr$V3))
assign(paste0("list_", SampleName), tmp_list)

# Bulk ASD-CTL subtraction Hi-C matrix
SampleName2 = "Batch1and2_ASD"
load(paste0(ASD_GENOVA_dir, SampleName2, "_loadcontacts.rda"))
load(paste0(ASD_TopDom_dir, SampleName2, "_TADcallsDomains_TopDom2GENOVA_windowsize",ws,".rda"))

SampleName_40kb_MAT = SampleName_40kb$MAT
SampleName_40kb_chr = SampleName_40kb_MAT %>%
  filter(V1 >= chr_range[1], V2 <= chr_range[2], V1 != V2)
# Clean up Diagonal
SampleName_40kb_chr$V3[which(SampleName_40kb_chr$V1 == SampleName_40kb_chr$V2)] = 0

SampleName_40kb_subset = SampleName_40kb_chr %>%
  filter(V1 >= min(ABS_subset$V4), V2 <= max(ABS_subset$V4))
colnames(SampleName_40kb_subset)[3] = "counts"
# Set diagonal as max
Diagonal_subset = data.frame(V1 = seq(min(ABS_subset$V4), max(ABS_subset$V4), 1), V2 = seq(min(ABS_subset$V4), max(ABS_subset$V4), 1), counts = max(SampleName_40kb_subset$counts))
SampleName_40kb_subset = rbind(SampleName_40kb_subset, Diagonal_subset)
# Add genomic coordinates
SampleName_40kb_subset$V4 = ABS_subset$V2[match(SampleName_40kb_subset$V1, ABS_subset$V4)] # Note that I should use V2 instead of mid
SampleName_40kb_subset$V5 = ABS_subset$V2[match(SampleName_40kb_subset$V2, ABS_subset$V4)]
SampleName_40kb_subset = SampleName_40kb_subset[,c(4,5,3,1,2)]

Domains_relevant2 = Domains[Domains$chrom == CHR & Domains$start > START & Domains$end < END,] # same as Domains_relevant

tmp_list = list(MAT = SampleName_40kb_subset, Domains_relevant = Domains_relevant2, SUM_CC = sum(SampleName_40kb_chr$V3))
assign(paste0("list_", SampleName2), tmp_list)

# Normalize for total contact count on this chr
nf_ASD = list_Batch1and2_ASD$SUM_CC/list_Batch1and2_CTL$SUM_CC
list_Batch1and2_ASD$MAT$counts = list_Batch1and2_ASD$MAT$counts/nf_ASD

SampleName_40kb_subtract = list_Batch1and2_CTL$MAT
SampleName_40kb_subtract$counts = list_Batch1and2_ASD$MAT$counts - list_Batch1and2_CTL$MAT$counts
tmp_list = list(MAT = SampleName_40kb_subtract, Domains_relevant = Domains_relevant2, SUM_CC = 0)
assign("list_Batch1and2_subtract", tmp_list)

# ATAC prep
CTL_ATAC_bw_subset = CTL_ATAC_bw %>%
  filter(chrom == CHR, end > START, start < END)
range(CTL_ATAC_bw_subset$score) # 0-5.12 for USP7 gene, 0-2.86 for PRICKLE1

# CpG island prep
DMR_subset = DMR %>%
  filter(chr == CHR, end > START, start < END)
colnames(DMR_subset)[1] = "chrom"

# CTCF prep
CTCF_subset = CTCF_inATAC %>%
  filter(chr == CHR, end > START, start < END)
colnames(CTCF_subset)[1] = "chrom"

# Hi-C heatmap color max
Gb_scale = "Mb"
Gene_highlight = candidate_clear$external_gene_name[i]
  
save(candidate_clear, i, CHR, START, END, list_Batch1and2_CTL, list_Batch1and2_ASD, list_Batch1and2_subtract, ws, TADb_weakened, CTCF_HALFwidth, CTL_ATAC_bw_subset, Gene, Gene_highlight, CTCF_subset, DMR_subset, CUTOFF_CTL, Gb_scale, TextSize, PointSize, DOMAIN_COL, pre_name, file = paste0(pre_name, "plotgardner_",Gene,"_Window_",CHR, "_",START/1000,"to",END/1000,"kb_CUTOFFctl",CUTOFF_CTL,"_DomainTopDom2GENOVAws",ws,".rda")) # _DomainTopDom2GENOVAws5
  
rm(list = ls())
lnames = load("06_plotgardner_UBR5_Window_chr8_102000to104520kb_CUTOFFctl800_DomainTopDom2GENOVAws5.rda")

left_blank = 0.5
pdf(paste0(pre_name, "plotgardner_",Gene,"_Window_",CHR, "_",START/1000,"to",END/1000,"kb_CUTOFFctl",CUTOFF_CTL,"_DomainTopDom2GENOVAws",ws,".pdf"), height = 3.2, width = 2.3 + left_blank)
pageCreate(width = 2.3 + left_blank, height = 3.2, showGuides = F)
  
region <- pgParams(
  chrom = CHR,
  chromstart = START, chromend = END,
  assembly = "hg19",
  width = 2
)

# CTL_Bulk
hicPlot_Bulk = plotHicTriangle(
  data = list_Batch1and2_CTL$MAT[,1:3], resolution = 40000,
  zrange = c(0, 800),
  #zrange = c(0, CUTOFF_CTL),
  params = region,
  x = left_blank + 0.15, y = 0.05, height = 1,
  just = c("left", "top"),
  default.units = "inches"
)
plotText(
  label = "Bulk_CTL", fontsize = TextSize, fontface = "bold",
  x = left_blank + 0.38, y = 0.35
)
annoHeatmapLegend(
  plot = hicPlot_Bulk, x = left_blank + 2, y = 0.05,
  width = 0.05, height = 0.6,
  just = c("left", "top"),
  fontcolor = "black",
  fontsize = TextSize
)
annoDomains(
  plot = hicPlot_Bulk, 
  data = list_Batch1and2_CTL$Domains_relevant, # list_Batch1and2_CTL$Domains_relevant is off by 1 40kb bin
  #data = manual_Domains,
  linecolor = DOMAIN_COL,  # darkgrey
  lwd = PointSize
)

# Subtraction matrix
hicPlot_Subtraction = plotHicTriangle(
  data = list_Batch1and2_subtract$MAT[,1:3], resolution = 40000,
  zrange = c(0, 35),
  params = region,
  x = left_blank + 0.15, y = 1.1, height = 1, # adjust y
)
plotText(
  label = "ASD - CTL", fontsize = TextSize, fontface = "bold",
  x = left_blank + 0.38, y = 1.35
)
annoHeatmapLegend(
  plot = hicPlot_Subtraction, x = left_blank + 2, y = 1.1,
  width = 0.05, height = 0.6,
  just = c("left", "top"),
  fontcolor = "black",
  fontsize = TextSize
)
annoDomains(
  plot = hicPlot_Subtraction,
  data = list_Batch1and2_subtract$Domains_relevant, # list_Batch1and2_CTL$Domains_relevant is off by 1 40kb bin
  #data = manual_Domains, # list_Batch1and2_CTL$Domains_relevant is off by 1 40kb bin
  linecolor = DOMAIN_COL,
  lwd = PointSize,
  lty = "dashed"
)

plotText(
  label = "Hi-C", fontsize = TextSize, fontface = "bold",
  x = 0.3, y = 1
)

# Genome labels
annoGenomeLabel(
  plot = hicPlot_Bulk, x = left_blank + 0.15, y = 2.1,
  scale = Gb_scale, # If not pretty, remove this line to allow full expression of xxx,xxx,xxx bp
  fontsize = TextSize, 
  just = c("left", "top")
)

nf = 2.3/6 # 0.38

# CpG islands (DMRs)
DMR_subset$col = 
  ifelse(DMR_subset$DMR == "hyper_methylated", "light salmon",
         ifelse(DMR_subset$DMR == "hyper_methylated", "steel blue","darkgrey"))
DMR_subset_expanded = DMR_subset
DMR_subset_expanded$start = DMR_subset$start - 2*CTCF_HALFwidth
DMR_subset_expanded$end = DMR_subset$end + 2*CTCF_HALFwidth

idx_DMR = which(DMR_subset_expanded$end > as.integer(unlist(str_split_fixed(TADb_weakened, "_", 3))[2]) & DMR_subset_expanded$start < as.integer(unlist(str_split_fixed(TADb_weakened, "_", 3))[3]) & DMR_subset_expanded$DMR == "hyper_methylated")
idx_DMR_rm = which(DMR_subset_expanded$end > DMR_subset_expanded$start[idx_DMR] & DMR_subset_expanded$start < DMR_subset_expanded$end[idx_DMR] & DMR_subset_expanded$DMR == "no_signif_change")
if(length(idx_DMR_rm) > 0) {DMR_subset_expanded = DMR_subset_expanded[-idx_DMR_rm,]}

pileupPlot <- plotRanges(
  data = DMR_subset_expanded,
  params = region,
  #fill = "grey",
  fill = DMR_subset_expanded$col,
  #fill = colorby("DMR", palette =
  #                 colorRampPalette(c("steel blue", "light salmon"))),
  x = left_blank + 0.15, y = 2.22, height = 0.17*nf,
  boxHeight = unit(3*nf, "mm"), spaceHeight = 0.4*nf,
  spaceWidth = 0.0005
)
plotText(
  label = "CpG islands", fontsize = TextSize, fontface = "bold",
  x = 0.3, y = 2.22 + 0.17*nf/2
)

# CTCF BS
idx_CTCF = which(CTCF_subset$end > as.integer(unlist(str_split_fixed(TADb_weakened, "_", 3))[2]) & CTCF_subset$start < as.integer(unlist(str_split_fixed(TADb_weakened, "_", 3))[3]))
CTCF_subset$col = "grey"
CTCF_subset$col[idx_CTCF] = "black"

CTCF_subset_expanded = CTCF_subset
CTCF_subset_expanded$start = CTCF_subset$start - CTCF_HALFwidth
CTCF_subset_expanded$end = CTCF_subset$end + CTCF_HALFwidth

pileupPlot <- plotRanges(
  data = CTCF_subset_expanded,
  params = region,
  #fill = "black",
  fill = CTCF_subset_expanded$col,
  x = left_blank + 0.15, y = 2.32, height = 0.17*nf,
  boxHeight = unit(3*nf, "mm"), spaceHeight = 0.4*nf,
  #spaceWidth = 0.00001
  spaceWidth = 0.005 # for UBR5
)
plotText(
  label = "CTCF BS", fontsize = TextSize, fontface = "bold",
  x = 0.3, y = 2.32 + 0.17*nf/2
)

# ATAC peaks
signalPlot <- plotSignal(
  data = CTL_ATAC_bw_subset, params = region,
  range = c(0, max(CTL_ATAC_bw_subset$score)),  # this sets the range of the signal
  x = left_blank + 0.15, y = 2.42 + 0.17*nf/2, width = 2, height = 0.3,
  linecolor = "black", baseline.color = "lightgrey"
)
plotText(
  label = "ATAC", fontsize = TextSize, fontface = "bold",
  x = 0.3, y = 2.57 + 0.17*nf/2
)

# Genes
geneOrder = Gene
genesPlot <- plotGenes(
  params = region,
  fill = c("black", "black"),
  fontcolor = c("black", "black"),
  fontsize = TextSize, 
  x = left_blank + 0.15, y = 2.75 + 0.17*nf/2, width = 2, height = 0.35, 
  # adjust y = 3.25 when Triangle baseline at 3
  # remember to adjust width to the same as Triangle Hi-C map
  geneHighlights = data.frame(gene = Gene_highlight, color = "red"), geneBackground = "black",
  geneOrder = geneOrder
) 
plotText(
  label = "Genes", fontsize = TextSize, fontface = "bold",
  x = 0.3, y = 2.75 + 0.17*(1+nf/2)
)

dev.off()



