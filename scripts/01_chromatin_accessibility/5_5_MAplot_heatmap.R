rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/")

## Load differential ATAC results
load("p_magnitude_lmeDiagnosisASD_ASD_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.rda")
DiffATAC = as.data.frame(cbind(magnitude[,1], p[,1]))
colnames(DiffATAC) = c("ATAC_logFC", "ATAC_p")
DiffATAC$ATAC_FDR = p.adjust(DiffATAC$ATAC_p, method = "fdr")

DiffATAC_location = as.data.frame(str_split_fixed(rownames(DiffATAC), "_", 3))
colnames(DiffATAC_location) = c("chr", "start", "end")
DiffATAC_location$start = as.integer(DiffATAC_location$start)
DiffATAC_location$end = as.integer(DiffATAC_location$end)
DiffATAC = cbind(DiffATAC_location, DiffATAC)

save(DiffATAC, file = "DiffATAC.rda")

## MA-plot
load("CQN.rda")
stopifnot(rownames(RPKM.cqn) == rownames(DiffATAC))
Mean = apply(RPKM.cqn, 1, mean)
MA_df = cbind(DiffATAC, Mean)

## Define differential ATAC peaks
for (fdr_t in c(0.05, 0.1)) {
  for (logfc_t in seq(0, 0.3, 0.1)) {
    MA_df$col = case_when(
      MA_df$ATAC_FDR >= fdr_t ~ "lightgrey",
      MA_df$ATAC_FDR < fdr_t & MA_df$ATAC_logFC > logfc_t ~ "red",
      MA_df$ATAC_FDR < fdr_t & MA_df$ATAC_logFC < -logfc_t ~ "blue"
    )
    
    tmp = as.data.frame(table(MA_df$col)) # blue (down) 2248, red (up) 2785, lightgrey 122938
    
    save(MA_df, file = paste0("MAdf_FDRthreshold",fdr_t, "_logFCthreshold",logfc_t,".rda"))
    MA_df = MA_df %>%
      arrange(desc(ATAC_FDR))
    
    pdf(paste0("MAplot_FDRthreshold",fdr_t, "_logFCthreshold",logfc_t,".pdf"), height = 4, width = 6)
    plot(MA_df$Mean, MA_df$ATAC_logFC, pch = 19, col = MA_df$col,
         xlab = "Log2 mean accessibility", ylab = "Log2 fold change",
         main = paste0("Total ", nrow(MA_df), " ATAC peaks"))
    text(x = 2, y = 1.2, paste0("Up-reg: ", tmp$Freq[tmp$Var1 == "red"], " peaks\nDown-reg: ", tmp$Freq[tmp$Var1 == "blue"], " peaks"))
    dev.off()
  }
}

## Used FDR < 0.05 

