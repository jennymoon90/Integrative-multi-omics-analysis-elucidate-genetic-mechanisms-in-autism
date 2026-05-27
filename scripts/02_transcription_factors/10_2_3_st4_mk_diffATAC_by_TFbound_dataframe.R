qrsh -l h_data=10G,h_rt=8:00:00 
module load R # Loading R/4.0.2
R

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/BINDetect/")
overview_dir = "TFoi_overview/"
outdir = "DiffATAC_intersectTFbound/"
DiffATAC_DIR="/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/5_DiffATAC"
ATACup = read.table(paste0(DiffATAC_DIR, "/promoter_proximal_HiCdistal_ATACup_ATACfdrthreshold0.1_logfcthreshold0.2.bed"))
ATACdown = read.table(paste0(DiffATAC_DIR, "/promoter_proximal_HiCdistal_ATACdown_ATACfdrthreshold0.1_logfcthreshold0.2.bed"))
str(ATACup) # 4841 obs of 11 variables
str(ATACdown) # 4165 obs of 11 variables
ATACup = unique(ATACup[,1:3]); ATACdown = unique(ATACdown[,1:3])
colnames(ATACup) = colnames(ATACdown) = c("chr", "start", "end")
str(ATACup) # 4758 obs of 3 variables
str(ATACdown) # 3997 obs of 3 variables

all_files = list.files(path = overview_dir, pattern = "Hsapiens.*_ASCorCTLbound.txt")
TFoi = read.table(paste0(overview_dir, "TFoi_topTFsFootprintConfirmed_total39.txt"))
TFbound_files = paste0(overview_dir, TFoi$V1, "_overview_ASCorCTLbound.txt")
all(TFbound_files %in% paste0(overview_dir, all_files)) # TRUE, note that ASDorCTLbound was named ASCorCTLbound by mistake in the previous script.

length(TFbound_files) # 39 TFs of interest - top 5% by differential binding score, human motif and footprint confirmed normal shape

TFn = 1
for (f in TFbound_files) {
  TFmotif = unlist(str_split(f, "-", 3))[3] # the first 3 is to ensure getting multiple motifs of the same TF.
  TFmotif = sub("_overview_ASDorCTLbound.txt", "", TFmotif)
  print(paste(TFn, TFmotif))
  TFn = TFn + 1
  
  df = read.table(f)
  df = unique(cbind(df[,7:9], 1)) # Doing the unique() is important, otherwise the resulting ATACup will be very big with many duplicated rows.
  colnames(df) = c("chr", "start", "end", TFmotif)

  ATACup = left_join(ATACup, df)
  ATACdown = left_join(ATACdown, df)
}
# 2022/11/14 9:45PM, takes < 1min

save(ATACup, ATACdown, file = paste0(outdir, "ATACupORdown_fdr01logfc02_overlap_TFoiBound.rda"))
# Transfer the rda file to MAC results/10_TOBIAS/BINDetect/

#####################
## Now work on MAC ##
#####################

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/BINDetect")

load("ATACupORdown_fdr01logfc02_overlap_TFoiBound.rda") # ATACup ATACdown
ATACup[is.na(ATACup)] = 0
ATACdown[is.na(ATACdown)] = 0

## 1) How many TFs bind in each differential ATAC peak?
ATACup$sum = apply(ATACup[,4:ncol(ATACup)], 1, sum)
ATACdown$sum = apply(ATACdown[,4:ncol(ATACdown)], 1, sum)

hist(ATACup$sum) 
hist(ATACdown$sum) 
# both mostly 0, if not show a gradient decrease probability of multiple TF binding.

diffATAC_sums = data_frame(DiffATAC = c(rep("Up", nrow(ATACup)), rep("Down", nrow(ATACdown))), N_TFbound = c(ATACup$sum, ATACdown$sum))
range(diffATAC_sums$N_TFbound) # 0-22
diffATAC_sums$DiffATAC = factor(diffATAC_sums$DiffATAC, levels = c("Up", "Down"))

pdf("Number_TFs_bound_inEaDiffATAC.pdf", height = 4, width = 6)
diffATAC_sums %>%
  ggplot(aes(x = N_TFbound, fill = DiffATAC)) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 1) +
  theme_bw() +
  xlab("Number of TFs bound in each differential ATAC peak") +
  ylab("Frequency (Number of ATAC peaks)") +
  annotate("text", x = 12, y = 1000, label = paste0("Range of the number of TFs \n bound in each diff. ATAC peak: 0-", range(diffATAC_sums$N_TFbound)[2]))
dev.off()
# More up-reg ATAC peaks than down-reg ones have 0 TF of interest bound. More down-reg ATAC peaks than up-reg ones have >= 1 TF of interest bound.

## 2) Which factors bind to a big proportion of the diffATAC peaks?
ATACconsensus = read.table("../../4_ConsensusPeaks/1_DiffBind/Peaksets_minOverlap5.bed")
N_ATACup_bound_byEaTF = read.table("ATACup_ATACfdrthreshold0.1_logfcthreshold0.2_overlapTF.summary")
N_ATACdown_bound_byEaTF = read.table("ATACdown_ATACfdrthreshold0.1_logfcthreshold0.2_overlapTF.summary")
N_ATACconsensus_bound_byEaTF = read.table("ATACconsensus_overlapTF.summary")

N_ATACup_bound_byEaTF$V1 = gsub("Hsapiens-jaspar2016-", "", N_ATACup_bound_byEaTF$V1)
N_ATACdown_bound_byEaTF$V1 = gsub("Hsapiens-jaspar2016-", "", N_ATACdown_bound_byEaTF$V1)
N_ATACconsensus_bound_byEaTF$V1 = gsub("Hsapiens-jaspar2016-", "", N_ATACconsensus_bound_byEaTF$V1)

N_ATACup_bound_byEaTF$V1 = gsub("Hsapiens-HOCOMOCOv10-", "", N_ATACup_bound_byEaTF$V1)
N_ATACdown_bound_byEaTF$V1 = gsub("Hsapiens-HOCOMOCOv10-", "", N_ATACdown_bound_byEaTF$V1)
N_ATACconsensus_bound_byEaTF$V1 = gsub("Hsapiens-HOCOMOCOv10-", "", N_ATACconsensus_bound_byEaTF$V1)

# Perform Fisher's Exact Test to assess association between TF binding and ATAC dysregulation (up and down respectively)
N_up_down = cbind(N_ATACdown_bound_byEaTF, N_ATACup_bound_byEaTF$V2, N_ATACconsensus_bound_byEaTF$V2)
colnames(N_up_down) = c("TFmotif", "N_ATACdown", "N_ATACup", "N_ATACconsensus")
TFdown = c("CTCF", "SOX9", "ETV5", "ELK4", "NR2F1") # "CPEB1" is associated with down-reg ATAC. However, its differential binding score is up.
N_up_down$type = "Up"
N_up_down$type[which(grepl(paste(TFdown, collapse = "|"), N_up_down$TFmotif))] = "Down"

for (i in 1:nrow(N_up_down)) {
  if (N_up_down$type[i] == "Down") {
    bound_ATACdown = N_up_down$N_ATACdown[i]
    unbound_ATACdown = nrow(ATACdown) - bound_ATACdown
    bound_ATACnotd = N_up_down$N_ATACconsensus[i] - bound_ATACdown
    unbound_ATACnotd = nrow(ATACconsensus) - nrow(ATACdown) - bound_ATACnotd
    
    fisher_matrix = matrix(c(bound_ATACdown, unbound_ATACdown, bound_ATACnotd, unbound_ATACnotd), nrow = 2)
    fisher_res = fisher.test(fisher_matrix)
    N_up_down$Fisher_p[i] = fisher_res$p.value
  } else {
    bound_ATACup = N_up_down$N_ATACup[i]
    unbound_ATACup = nrow(ATACup) - bound_ATACup
    bound_ATACnotu = N_up_down$N_ATACconsensus[i] - bound_ATACup
    unbound_ATACnotu = nrow(ATACconsensus) - nrow(ATACup) - bound_ATACnotu
    
    fisher_matrix = matrix(c(bound_ATACup, unbound_ATACup, bound_ATACnotu, unbound_ATACnotu), nrow = 2)
    fisher_res = fisher.test(fisher_matrix)
    N_up_down$Fisher_p[i] = fisher_res$p.value
  }
}

N_up_down$Fisher_fdr = p.adjust(N_up_down$Fisher_p, method = "fdr")
N_up_down$Fisher_pbonf = p.adjust(N_up_down$Fisher_p, method = "bonferroni")
N_up_down$TFmotif[N_up_down$Fisher_pbonf < 0.05] # 11 
N_up_down$TFmotif[N_up_down$Fisher_fdr < 0.05] # 15 including CTCF, (ELK4), SOX9, NR2F1, JDP2, JUNB, FOS, JUN, FOSL1, (BATF, NFE2, JUNB), JUND, CPEB1. ETV5 is just above 0.05. No RFX2/3/4/5. The ones in brackets are only significant using fdr but not bonferroni. Since they are not independent (some belong to the same TF cluster), makes more sense to use FDR.

# Calculate percentage
colnames(N_ATACup_bound_byEaTF) = colnames(N_ATACdown_bound_byEaTF) = c("TFmotif", "N_diffATAC")
N_ATACup_bound_byEaTF$Perc_diffATAC = N_ATACup_bound_byEaTF$N_diffATAC/nrow(ATACup)
N_ATACdown_bound_byEaTF$Perc_diffATAC = N_ATACdown_bound_byEaTF$N_diffATAC/nrow(ATACdown)

N_ATACup_bound_byEaTF = N_ATACup_bound_byEaTF %>%
  arrange(desc(N_diffATAC))
N_ATACdown_bound_byEaTF = N_ATACdown_bound_byEaTF %>%
  arrange(desc(N_diffATAC))

N_ATACup_bound_byEaTF$DiffATAC = "Up"; N_ATACdown_bound_byEaTF$DiffATAC = "Down"
diffATAC_eaTF = rbind(N_ATACup_bound_byEaTF, N_ATACdown_bound_byEaTF)
diffATAC_eaTF$DiffATAC = factor(diffATAC_eaTF$DiffATAC, levels = c("Up", "Down"))
 
tmp = data_frame(TFmotif = N_ATACup_bound_byEaTF$TFmotif, N_upATAC = N_ATACup_bound_byEaTF$N_diffATAC, N_downATAC = N_ATACdown_bound_byEaTF$N_diffATAC[match(N_ATACup_bound_byEaTF$TFmotif, N_ATACdown_bound_byEaTF$TFmotif)],
                 Perc_upATAC = N_ATACup_bound_byEaTF$Perc_diffATAC, Perc_downATAC = N_ATACdown_bound_byEaTF$Perc_diffATAC[match(N_ATACup_bound_byEaTF$TFmotif, N_ATACdown_bound_byEaTF$TFmotif)])
tmp$ratio_UpOverDown = tmp$Perc_upATAC/tmp$Perc_downATAC

diffATAC_eaTF2 = tmp; rm(tmp)
diffATAC_eaTF2 = diffATAC_eaTF2 %>%
  arrange(diffATAC_eaTF2$Perc_upATAC) # ratio_UpOverDown does not look good, N_upATAC looks more ordered
  # arrange(diffATAC_eaTF2$Perc_downATAC) # not as clear as ranking by Perc_upATAC
diffATAC_eaTF$TFmotif = factor(diffATAC_eaTF$TFmotif, levels = diffATAC_eaTF2$TFmotif)

# Add Fisher's Exact p.adjust (fdr)
diffATAC_eaTF$Fisher_fdr = 1
for (i in 1:nrow(diffATAC_eaTF)) {
  idx = which(N_up_down$TFmotif == diffATAC_eaTF$TFmotif[i] & N_up_down$type == diffATAC_eaTF$DiffATAC[i])
  if (length(idx) > 0) {
    diffATAC_eaTF$Fisher_fdr[i] = N_up_down$Fisher_fdr[idx]
  }
}
diffATAC_eaTF$Fisher_label = ifelse(diffATAC_eaTF$Fisher_fdr < 0.05, "*", ifelse(diffATAC_eaTF$Fisher_fdr < 0.1, "#", ""))
diffATAC_eaTF$Fisher_nudge = ifelse(diffATAC_eaTF$DiffATAC == "Up", -0.25, 0.25)
diffATAC_eaTF$Fisher_size = case_when(
  diffATAC_eaTF$Fisher_label == "*" ~ 20,
  diffATAC_eaTF$Fisher_label == "#" ~ 10,
  is.na(diffATAC_eaTF$Fisher_label) ~ 0
  )


pdf("Proportion_DiffATAC_bound_byEaTF.pdf", height = 5, width = 8)
diffATAC_eaTF %>%
  ggplot(aes(x = TFmotif, y = Perc_diffATAC, fill = DiffATAC)) + # y = N_diffATAC
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 0.25) + # it does help to view the little ones more clearly
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 25, unit = "pt")) +
  xlab("TF motifs") +
  ylab("Proportion of diff. ATAC peaks bound by the TF") + # Number of
  geom_text(aes(label = Fisher_label, y = Perc_diffATAC + 0.005), nudge_x = diffATAC_eaTF$Fisher_nudge, size = diffATAC_eaTF$Fisher_size / .pt) + # label = ifelse(Fisher_fdr < 0.05, "*", "")
  ylim(0, max(diffATAC_eaTF$Perc_diffATAC) + 0.02) +
  annotate("text", x = 2, y = 0.15, label = "Association b/w TF binding and ATAC dysreg.\n* FDR < 0.05\n# FDR < 0.1", hjust = 0) # hjust = 0 to left-align
dev.off()
# Observation: 
# a. CTCF bind to a large proportion (nearly 7.5-20%) of both up and down-reg ATAC peaks.
# b. CTCF, SOX9, NR2F1, ELK4 and CPEB1 show evidently greater binding to down-reg ATAC peaks. ETV5 is just below the FDR<0.05 threshold
# c. The TFs with up-reg binding score do not show higher proportion of up-reg ATAC peaks bound by them compared to down-reg ATAC peaks.
# Note the asterisks is for significance level from Fisher's exact test on association between TF binding and peak dysregulation (bound vs unbound, down-reg vs not, contingency table).

save(diffATAC_eaTF, file = "Proportion_DiffATAC_bound_byEaTF_and_FisherExactTestOnDependencyBwTFbindingAndATACdysreg.rda")

## 3) Which factors co-bind?
# Use 10_2_3_st5_VennDiagraom_downregATAC_coBoundbyTFoi.R

