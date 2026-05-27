rm(list = ls())
setwd("/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_3_14_GeneticVariants_at_CTCFboundATAC_withGEdownreg/")
library(motifbreakR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

args = commandArgs(trailingOnly = T)
i = args[1]

lnames = load("../10_3_13_GeneticVariants_at_CTCFpromoter/mf_jaspar2022hs.rda") # mf
lnames = load(paste0("snps_mb_CTCFboundDEGdownSNPs_dbSNP155GRCh37_",i,".rda")) # snps.mb
print(length(snps.mb))

## st3. motifbreakR
results = motifbreakR(snpList = snps.mb, filterp = TRUE, 
                      pwmList = mf, # mf is from jaspar 2022
                      threshold = 1e-4)
# R/4.2.2 Error in d$value$value : $ operator is invalid for atomic vectors

save(results, file = paste0("motifbreakRresults_jaspar2022hsmotifs_CTCFboundDEGdownSNPs_dbSNP155GRCh37_",i,".rda"))

## st4. Filter for "strong" predicted effects
results_names = names(results) # 2678 rsids
uniq_results_names = unique(names(results)) # 568 rsids

df_results = as.data.frame(unname(results)) # 2678 obs x 24 columns
strong_results = df_results[df_results$effect == "strong", ] # 1831 rows
print(length(unique(strong_results$SNP_id))) # 449 SNPs

save(results, df_results, strong_results, file = paste0("STRONGmotifbreakRresults_jaspar2022hsmotifs_CTCFboundDEGdownSNPs_dbSNP155GRCh37_",i,".rda"))

