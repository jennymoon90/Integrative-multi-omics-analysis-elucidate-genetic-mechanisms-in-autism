### Get top PCs that explain total 80% variance in the multiHiCcompare normalized logCPM, identify biological and technical covariates that significantly correlate with these top PCs.

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/23_differential_loop_using_package/mimickATAC_lmBuild/")

### Load data ### 
## Load un-normalized loop logCPM data
lnames = load("../../25_nlme_DiffLoopAnalysis/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_includingOutlier.rda") # "Bulk_logCPM" "NeuNp_logCPM" "NeuNn_logCPM" from 25_3_Loop_logCPM_includingOutliers.R # old: 25_1_Loop_logCPM_includingTechReps.R

## Load Covariates
#lnames = load("../multiHiCcompare_real_ASD_eliminateSegDup_distanceOver0/hicexp_normalized_chr1_noTechRep_rmSegDup_distOver0.rda"); rm(hicexp)
lnames = load("../../25_nlme_DiffLoopAnalysis/datMetaSeq.rda") # "datMetaSeq" "Sample_rm" from 25_1_Loop_logCPM_includingTechReps.R

colnames(Bulk_logCPM)[which(colnames(Bulk_logCPM) == "B5242A")] = "B5242"
colnames(Bulk_logCPM)[which(colnames(Bulk_logCPM) == "B5242B")] = "B5242_Rep"
colnames(Bulk_logCPM)[which(colnames(Bulk_logCPM) == "B5342A")] = "B5342"
colnames(Bulk_logCPM)[which(colnames(Bulk_logCPM) == "B5342B")] = "B5342_Rep"
Covariates = datMetaSeq[match(colnames(Bulk_logCPM), rownames(datMetaSeq)),] # 40 obs. of 19 variables

colnames(Covariates)
apply(Covariates[,13:15], 1, sum) # all sum up to 1.

# Format biological and technical Covariates
Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "ASD"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal", "Temporal"))
Covariates$PMI[is.na(Covariates$PMI)] = mean(Covariates$PMI, na.rm = T)
Covariates = Covariates[,c("Subject", "Diagnosis", "Age", "PMI", "Sex", "BrainBank", "Batch", "Cortex", "Region", "Valid", "Duplicate", "Trans", "cis_shortRange", "cis_longRange", "ReadDepth")]
colnames(Covariates) = gsub("_", ".", colnames(Covariates))

## Identify top PCs
norm <- t(scale(t(Bulk_logCPM),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:27)]) ## Top 27 PCs explain >80% of total variance.
topPC <- PC$rotation[,1:27]

### Correlation ###
## Correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates)[-1]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
colnames(Cor)[ncol(topPC) + c(1, 4:14)]
idx_spearman = ncol(topPC) + c(1, 4:14)
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

## Find out which Covariates significantly correlate with the topPCs (like what Yuyan described in Feng 2022 BioRxiv)
Cor_sig = matrix(nrow = nrow(Cor), ncol = ncol(Cor))
rownames(Cor_sig) = colnames(Cor_sig) = colnames(Cor)
for (i in 1:ncol(mod_mat_withPC)) {
  for (j in 1:ncol(mod_mat_withPC)) {
    tmp = cor.test(mod_mat_withPC[,i], mod_mat_withPC[,j])
    Cor_sig[i,j] = tmp$p.value
    # if (tmp$p.value < 0.05) {print(paste(i, colnames(Cor)[i], j, colnames(Cor)[j]))}
  }
}

n_tests = (ncol(Cor_sig)-1)^2 - (ncol(topPC)-1)^2
(p_cor_threshold = 0.05/n_tests)
Cov_pool_idx = which(apply(Cor_sig[1:ncol(topPC),], 2, function(x) any(x < p_cor_threshold)))
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # Duplicate, Trans, cis.longRange

# Bonferroni correction of the p-val is too stringent, use fdr and see how many covariates are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisASD! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.1)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Still no DiagnosisASD! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.2)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Finally there is DiagnosisASD, similar to ATAC

## Corrplot (the 2nd plot of each denotes fdr<0.2 by *)
library(corrplot)
sig_level = 0.2
# change the old corrplot name as xxx_FrontalAsBase
pdf(paste0("Corrplot_topPCofBulklogCPM_Covariates_includingTechRepsAndOutliers.pdf"), height = 35, width = 45) 
corrplot(Cor,method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor, p.mat = Cor_fdr, insig = "label_sig", sig.level = sig_level, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = sig_level, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = sig_level, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
dev.off()

### Select covariates ###
## I know that Batch/Cortex/Region/BrainBank are correlated
unique(Covariates[,c("Cortex", "Region", "Batch")])
# Frontal       BA9         2
# Temporal      BA38        2
# Parietal      BA7         1
# Parietal      BA3_1_2_5   1
# Frontal       BA44_45     2
# Parietal      BA39_40     1
unique(Covariates[,c("BrainBank", "Batch")])
# ABN     2
# NICHD-BTB     2
# NICHD-BTB     1
# Harvard-ATP     1
# Harvard-ATP     2
unique(Covariates[,c("Cortex", "Region", "Batch", "BrainBank")])
# Frontal        BA9     2         ABN
# Temporal      BA38     2   NICHD-BTB
# Parietal       BA7     1   NICHD-BTB
# Parietal BA3_1_2_5     1   NICHD-BTB
# Parietal       BA7     1 Harvard-ATP
# Frontal    BA44_45     2   NICHD-BTB
# Frontal        BA9     2   NICHD-BTB
# Temporal      BA38     2 Harvard-ATP
# Parietal BA3_1_2_5     1 Harvard-ATP
# Parietal   BA39_40     1 Harvard-ATP

Cov_pool
# [1] "DiagnosisASD"         "SexF"                 "BrainBankHarvard-ATP"
# [4] "BrainBankNICHD-BTB"   "Batch"                "CortexFrontal"       
# [7] "RegionBA39_40"        "RegionBA44_45"        "RegionBA7"           
# [10] "RegionBA9"            "Valid"                "Duplicate"           
# [13] "Trans"                "cis.shortRange"       "cis.longRange"  

## First include all candidate covariates that cor > 0.3 with PC1-PC4 and that are not significantly correlated with each other.
Cov_pool[sapply(Cov_pool, function(x) any(abs(Cor[1:4,x]) > 0.3))]
# [1] "DiagnosisASD"         "BrainBankHarvard-ATP" "BrainBankNICHD-BTB"  
# [4] "Batch"                "CortexFrontal"        "RegionBA7"           
# [7] "RegionBA9"            "Valid"                "Duplicate"           
# [10] "Trans"                "cis.shortRange"       "cis.longRange"                      
Cov_pool[!sapply(Cov_pool, function(x) any(abs(Cor[1:4,x]) > 0.3))]
# These ones later: "SexF" "RegionBA39_40" "RegionBA44_45"     

# Prioritize biological before technical covariate.
test_row = which(rownames(Cor_sig) == "DiagnosisASD")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include DiagnosisASD

test_row = which(rownames(Cor_sig) == "BrainBankHarvard-ATP")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
test_row = which(rownames(Cor_sig) == "BrainBankNICHD-BTB")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include BrainBank

test_row = which(rownames(Cor_sig) == "Batch")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # CortexFrontal, RegionBA7, RegionBA9, Duplicate
test_row = which(rownames(Cor_sig) == "CortexFrontal")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Batch, RegionBA9, Duplicate
test_row = which(rownames(Cor_sig) == "RegionBA7")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Batch, Duplicate
test_row = which(rownames(Cor_sig) == "RegionBA9")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) #  Batch, CortexFrontal, Duplicate

sum(abs(Cor[1:4,which(colnames(Cor) == "Batch")])) # 1.862537
sum(abs(Cor[1:4,which(colnames(Cor) == "CortexFrontal")])) # 1.339502
sum(abs(Cor[1:4,which(colnames(Cor) == "RegionBA7")])) # 1.058884
sum(abs(Cor[1:4,which(colnames(Cor) == "RegionBA9")])) # 1.44846
sum(abs(Cor[1:4,which(colnames(Cor) == "Duplicate")])) # 1.387015
# Include Batch among the 4, and this would exclude Duplicate.

# Now technical covariates:
test_row = which(rownames(Cor_sig) == "Valid")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Valid

test_row = which(rownames(Cor_sig) == "Trans")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PC4, cis.shortRange, cis.longRange
test_row = which(rownames(Cor_sig) == "cis.shortRange")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Trans, cis.longRange

sum(abs(Cor[1:4,which(colnames(Cor) == "Trans")])) # 1.364515
sum(abs(Cor[1:4,which(colnames(Cor) == "cis.shortRange")])) # 1.23334
sum(abs(Cor[1:4,which(colnames(Cor) == "cis.longRange")])) # 1.314903
# Include Trans

# Now come to the less important but significant ones
test_row = which(rownames(Cor_sig) == "SexF")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include SexF

test_row = which(rownames(Cor_sig) == "RegionBA39_40")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # None
# Include RegionBA39_40

test_row = which(rownames(Cor_sig) == "RegionBA44_45")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include RegionBA44_45

## Check VIF
library(olsrr)
# Base model
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans, data = cur_data)"
# All VIF < 2

# Will BrainBank cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard, data = cur_data)"
# All VIF < 2
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankNICHD, data = cur_data)"
# All VIF < 2
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + BrainBankNICHD, data = cur_data)"
# VIF of Batch and BrainBankHarvard > 2.5.
sum(abs(Cor[1:4,which(colnames(Cor) == "BrainBankHarvard-ATP")])) # 0.7398315
sum(abs(Cor[1:4,which(colnames(Cor) == "BrainBankNICHD-BTB")])) # 0.7128782
# Include just BrainBankHarvard

# Will Sex cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + SexF, data = cur_data)"
# All VIF < 2. Fine to include Sex

# Will RegionBA39_40 or RegionBA44_45 cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + RegionBA39_40, data = cur_data)"
# All VIF < 2.
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + RegionBA44_45, data = cur_data)"
# All VIF < 2.
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + RegionBA39_40 + RegionBA44_45, data = cur_data)"
# All VIF < 2. Fine to include both BA regions

# Will combining Sex and BA regions cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + RegionBA39_40 + RegionBA44_45 + SexF, data = cur_data)"
# All VIF < 2. Fine to include both Sex and the 2 BA regions

i=1
cur_data = Bulk_logCPM[i,]
cur_data = as.data.frame(cbind(t(cur_data), mod_mat))
colnames(cur_data)[1] = c("logCPM")
colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
colnames(cur_data)[which(colnames(cur_data) == "BrainBankHarvard-ATP")] = "BrainBankHarvard"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign

## Final model: 
expression_lm = "lm(logCPM ~ DiagnosisASD + Batch + Valid + Trans + BrainBankHarvard + SexF + RegionBA39_40 + RegionBA44_45, data = cur_data)"

## Run linear regression first (later lmm). 
runlm <- function(cur_data,expression_lm) {
  lm1 <- eval(parse(text = expression_lm));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients
  return(tabOut)
}

Covariates$RegionBA39_40 = ifelse(Covariates$Region == "BA39_40", 1, 0)
Covariates$RegionBA44_45 = ifelse(Covariates$Region == "BA44_45", 1, 0)
Covariates$BrainBankHarvard = ifelse(Covariates$BrainBank == "Harvard-ATP", 1, 0)
Covariates$RegionBA38 = ifelse(Covariates$Region == "BA38", 1, 0)

expression_lm = gsub("DiagnosisASD", "Diagnosis", expression_lm)
expression_lm = gsub("SexF", "Sex", expression_lm)
expression_lm

#expression_lm = "lm(logCPM ~ Diagnosis + Age + Batch + Valid + Trans + BrainBankNICHD, data = cur_data)"

n = length(unlist(str_split(expression_lm, "\\+")))
p = magnitude = matrix(nrow = nrow(Bulk_logCPM), ncol = n) 

for (i in 1:nrow(Bulk_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = Bulk_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  tabOut <- try(runlm(cur_data,expression_lm),silent=F)
  magnitude[i,] <- tabOut[-1,1]
  p[i,] <- tabOut[-1,4]
}
# 2022/09/08 4:42 PM

tabOut
colnames(p) = colnames(magnitude) = rownames(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(Bulk_logCPM)

e1 = unlist(str_split(expression_lm, "~ "))[2]
(e2 = unlist(str_split(e1, ","))[1])
e3 = gsub("_", ".", e2)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_Bulk_lm",e3,"_includingTechRepsAndOutliers.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {pos = nchar(e2)}

pdf(paste0("HistPval_Bulk_lm",e3,"_includingTechRepsAndOutliers.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the mis-chosen lm (DiagnosisBatchValidcis.shortRangeNICHDBA39.40BA44.45): Looks very nice! -> now run lmm
# Use the correctly-chosen lm (DiagnosisBatchValidTransRangeATPBA39.40BA44.45): Only a little enriched at p<0.05. Run lmm or Try substitute BrainBankHarvard with BrainBankNICHD.

## Run linear regression - using linear mixed model, including random effect of subjects
library(nlme)
runlme <- function(thisdat,expression) {
  lm1 <- eval(parse(text=expression));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients$fixed
  lm1.anova = anova(lm1)
  return(list(tabOut, lm1.anova))
}

# Variable twisting: 
#expression_lm = "lm(logCPM ~ Diagnosis + Batch + Valid + Trans + BrainBankNICHD + RegionBA39_40 + RegionBA44_45 + Sex, data = cur_data)"
#expression_lm = "lm(logCPM ~ Diagnosis + Batch + Valid + Trans + BrainBankNICHD + Sex + RegionBA39_40 + RegionBA44_45 + RegionBA38, data = cur_data)"

expression_model = gsub("lm","lme",expression_lm)
expression_model = gsub(",",", random=~1|Subject,", expression_model)
expression_model

n = length(unlist(str_split(expression_model, "\\+")))
p = magnitude = matrix(nrow = nrow(Bulk_logCPM), ncol = n) 

for (i in 1:nrow(Bulk_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = Bulk_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  lm1.out <- try(runlme(cur_data,expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of ATAC peak', i, rownames(Bulk_logCPM)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}
# 2022/09/08 6:05-6:10 PM (Use the original lm with 8 variables)
## 2022/09/08 5:19-5:24 PM (Use Trans instead of cis.shortRange, still 8 variables)
## 2022/09/08 5:34-5:39 PM (Use Trans instead of cis.shortRange and add BA38, 9 variables)
# Many Errors.
length(which(is.na(p[,1]))) 
# With the original lme using cis.shortRange: 83 Bulk loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Age + Batch + Valid + Trans +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# With Trans instead of cis.shortRange: 77 Bulk loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Age + Batch + Valid + Trans +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# With Trans instead of cis.shortRange, and add BA38: 81 Bulk loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Age + Batch + Valid + Trans +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(Bulk_logCPM)

e1 = unlist(str_split(expression_model, "~ "))[2]
(e2 = unlist(str_split(e1, ", data"))[1])
e3 = unlist(str_split(e2, ","))[1]
e3 = gsub("_", ".", e3)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_Bulk_lme",e3,"_includingTechRepsAndOutliers.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {pos = nchar(e2)}

pdf(paste0("HistPval_Bulk_lme",e3,"_includingTechRepsAndOutliers.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear mixed model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the original lm: Looks very nice! Better than lm(). This is similar to ATAC results.

range(p[,1], na.rm = T) # 7.3e-6 to 1
length(which(p[,1] < 0.05)) # 2369 nominally significant loops, more than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
fdr = p.adjust(p[,1])
range(fdr, na.rm = T) # 0.25 to 1. Better than not including the 2 ASD outliers. And better than using BrainBankNICHD.

## Check correlation between topPCs of regressed logCPM with covariates
colnames(magnitude)
colnames(mod_mat)
colnames(magnitude)[which(colnames(magnitude) == "BrainBankNICHD")] = "BrainBankNICHD-BTB"
colnames(magnitude)[which(colnames(magnitude) == "BrainBankHarvard")] = "BrainBankHarvard-ATP"
mod_used = mod_mat[,colnames(magnitude)[-1]]

logCPM.Regressed = Bulk_logCPM - as.matrix(magnitude[,-1]) %*% t(mod_used)
idx = which(apply(logCPM.Regressed, 1, function(x) all(is.na(x))))
if (length(idx) > 0) {logCPM.Regressed = logCPM.Regressed[-idx,]}

norm <- t(scale(t(logCPM.Regressed),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:23)]) # Top 23 PCs explain >80% of variance.
topPC <- PC$rotation[,1:23] ## these first 15 explain ~80% of the variance

mod_mat_withPC = cbind(topPC, mod_mat)
Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
idx_spearman = ncol(topPC) + c(1, 4:14)
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

pdf(paste0("Corrplot_TopPCofBulkRegressedLogCPM_Covariates_lme",e3,"_includingTechRepsAndOutliers.pdf"), height = 30, width = 40)
corrplot(Cor[1:ncol(topPC),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
dev.off()
# Original lme: Very nice and clean! Not sure whether adding Age would help. 

## Final model for Bulk:
expression_model
# "lme(logCPM ~ Diagnosis + Batch + Valid + Trans + BrainBankHarvard + Sex + RegionBA39_40 + RegionBA44_45, random=~1|Subject, data = cur_data)"
