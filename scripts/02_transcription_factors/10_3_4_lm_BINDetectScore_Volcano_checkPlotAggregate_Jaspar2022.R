rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(olsrr)

#setwd("/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect/all_bindetect_results")
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/ind_BINDetect_Jaspar2022/")
lnames = load("Jaspar2022hs_TFmotifs_mean_score_by_sample.rda") # df
colnames(df)[2:5] = gsub("X", "", colnames(df)[2:5])
colnames(df) = gsub("_mean_score", "", colnames(df))
colnames(df) = gsub("\\.", "-", colnames(df))

df$motif_id = gsub("Hsapiens-", "", df$motif_id)
rownames(df) = df$motif_id
df_Scores = df[,-1]; rm(df)

## Load Covariates
# lnames = load("/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/3_QC/BiolTechCov_and_QC.rda") # df_wAllCov
lnames = load("../../3_QC/BiolTechCov_and_QC.rda") # df_wAllCov
Covariates = df_wAllCov; rm(df_wAllCov)

all(colnames(df_Scores) %in% Covariates$sample_name) # TRUE
Covariates = Covariates[match(colnames(df_Scores), Covariates$sample_name),] # 38 samples of 48 columns
colnames(Covariates) = gsub("_", ".", colnames(Covariates))
colnames(Covariates)[which(colnames(Covariates) == "subject")] = "Subject"

Covariates$PMI = as.numeric(Covariates$PMI) 
Covariates$PMI[is.na(Covariates$PMI)] = mean(Covariates$PMI, na.rm = T)
Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "ASD"))
Covariates$Sex = factor(Covariates$Sex, levels = c("M", "F"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal", "Temporal"))

# Limit biol and tech cov to the relevant ones selected based on the "unorganized" script
# copied from 5_4_v7_PCA_build_linear_model_SignificantCorWithPC_VIF2p5_lmeSubject_organized.R
Covariates = Covariates[,c("Subject","Diagnosis", "Age", "PMI", "Sex", "BrainBank", "Batch", "Cortex", "Region", "tssenrich.score", "FRiP", "picard.alignment.PCT.READS.ALIGNED.IN.PAIRS", "picard.duplication.PERCENT.DUPLICATION")]
colnames(Covariates)[(ncol(Covariates)-1):ncol(Covariates)] = c("PCT.READS.ALIGNED.IN.PAIRS", "PERCENT.DUPLICATION")

## Remove technical replicates for building linear model
# But when I run lme, I shall include Subject as a random effect
TechRep_rm = c("B4334-1", "B4721A", "B5000B", "B5718D", "B5813B", "CQ56-2")
df_Scores2 = df_Scores[,which(! colnames(df_Scores) %in% TechRep_rm)]
Covariates2 = Covariates[which(! rownames(Covariates) %in% TechRep_rm),]

stopifnot(rownames(Covariates) == colnames(df_Scores))
stopifnot(rownames(Covariates2) == colnames(df_Scores2))

## Identify top PCs
norm <- t(scale(t(df_Scores2),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:2)]) ## Top 2 PCs explain 97% of total variance.
topPC <- PC$rotation[,1:2]

### Correlation ###
## Correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates2)[-1]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates2)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
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
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # tssenrich.score FRiP
# Diagnosis does not significantly correlate with any topPCs. Fine.

# Bonferroni correction of the p-val is too stringent, use fdr and see how many covariates are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisASD! Batch, CortexFrontal, RegionBA9, tssenrich.score, FRiP, PCT.READS.ALIGNED.IN.PAIRS, PERCENT.DUPLICATION.
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.1)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Still no DiagnosisASD! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.35)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Finally there is DiagnosisASD! Essentially all covs are included
SIG_T = 0.35

## Corrplot (the 2nd plot of each denotes fdr<0.2 by *)
library(corrplot)
# change the old corrplot name as xxx_FrontalAsBase
pdf(paste0("Corrplot_topPCofBINDetectMeanScore_Covariates.pdf"), height = 35, width = 45) 
corrplot(Cor,method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor, p.mat = Cor_fdr, insig = "label_sig", sig.level = SIG_T, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = SIG_T, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = SIG_T, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
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

Cov_pool # all covariates except for SexF

## First include all candidate covariates that cor > 0.5 with PC1 or PC2 and that are not significantly correlated with each other.
# Prioritize biological factors to technical ones

test_row = which(rownames(Cor_sig) == "DiagnosisASD")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Diagnosis

test_row = which(rownames(Cor_sig) == "Age")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Age

test_row = which(rownames(Cor_sig) == "PMI")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include PMI

test_row = which(rownames(Cor_sig) == "Batch")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # "CortexFrontal" "RegionBA7"  "RegionBA9" "PCT.READS.ALIGNED.IN.PAIRS" "PERCENT.DUPLICATION"
test_row = which(rownames(Cor_sig) == "CortexFrontal")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # "Batch" "RegionBA9"  "PERCENT.DUPLICATION"
test_row = which(rownames(Cor_sig) == "RegionBA9")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # "Batch" "CortexFrontal"  "PERCENT.DUPLICATION"
test_row = which(rownames(Cor_sig) == "PCT.READS.ALIGNED.IN.PAIRS")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # "Batch" "PERCENT.DUPLICATION"
test_row = which(rownames(Cor_sig) == "PERCENT.DUPLICATION")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # "Batch" "CortexFrontal"              "RegionBA9" "PCT.READS.ALIGNED.IN.PAIRS"

sum(abs(Cor[1:2,which(colnames(Cor) == "Batch")]) * varexp[1:2]) # 0.384493
sum(abs(Cor[1:2,which(colnames(Cor) == "CortexFrontal")]) * varexp[1:2]) # 0.3884356
sum(abs(Cor[1:2,which(colnames(Cor) == "RegionBA9")]) * varexp[1:2]) # 0.4542132
sum(abs(Cor[1:2,which(colnames(Cor) == "PCT.READS.ALIGNED.IN.PAIRS")]) * varexp[1:2]) # 0.242287
sum(abs(Cor[1:2,which(colnames(Cor) == "PERCENT.DUPLICATION")]) * varexp[1:2]) # 0.3814582
# Initial decision: Include RegionBA9, exempt "CortexFrontal" "Batch" "PERCENT.DUPLICATION"
# May try: Include Batch, exempt "CortexFrontal" "RegionBA7"  "RegionBA9" "PCT.READS.ALIGNED.IN.PAIRS" "PERCENT.DUPLICATION"

Cov_pool

test_row = which(rownames(Cor_sig) == "RegionBA38")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # CortexTemporal, essentially the same
# Include RegionBA38
test_row = which(rownames(Cor_sig) == "RegionBA39_40")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include RegionBA39_40
test_row = which(rownames(Cor_sig) == "RegionBA44_45")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include RegionBA44_45
test_row = which(rownames(Cor_sig) == "RegionBA7")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Batch
# May Include RegionBA7

test_row = which(rownames(Cor_sig) == "tssenrich.score")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PC1, PC2, FRiP
test_row = which(rownames(Cor_sig) == "FRiP")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PC1, tssenrich.score

sum(abs(Cor[1:2,which(colnames(Cor) == "tssenrich.score")]) * varexp[1:2]) # 0.6219323
sum(abs(Cor[1:2,which(colnames(Cor) == "FRiP")]) * varexp[1:2]) # 0.6073503
# Initial decision: Take in tssenrich.score and may spare FRiP if VIF > 2.5
# Choose between tssenrich.score and FRiP. Take in FRiP and may spare tssenrich.score if VIF > 2.5

test_row = which(rownames(Cor_sig) == "BrainBankHarvard-ATP")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include BrainBankHarvard-ATP
test_row = which(rownames(Cor_sig) == "BrainBankNICHD-BTB")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # None
# Include BrainBankNICHD-BTB

## Check VIF
library(olsrr)
# Base model
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score, data = cur_data)"
# all VIF < 2.5 (note that BA9 has VIF > 2)
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + RegionBA7, data = cur_data)"
# VIF RegionBA9 > 3, can't include BA7
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + BrainBankNICHD, data = cur_data)"
# all VIF < 2.5, fine to add BrainBankNICHD
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + BrainBankHarvard, data = cur_data)"
# VIF RegionBA9 > 3, can't include BrainBankHarvard
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + PCT.READS.ALIGNED.IN.PAIRS, data = cur_data)"
# VIF RegionBA9 & PCT.READS.ALIGNED.IN.PAIRS > 4, can't include PCT.READS.ALIGNED.IN.PAIRS


i=1
cur_data = df_Scores2[i,]
cur_data = as.data.frame(cbind(t(cur_data), mod_mat))
colnames(cur_data)[1] = c("score")
colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
colnames(cur_data)[which(colnames(cur_data) == "BrainBankHarvard-ATP")] = "BrainBankHarvard"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign

## Final model: 
expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + BrainBankNICHD, data = cur_data)"
# expression_lm = "lm(score ~ DiagnosisASD + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + FRiP, data = cur_data)"

## Run linear regression - using linear mixed model, including random effect of subjects
library(nlme)

runlme <- function(thisdat,expression) {
  lm1 <- eval(parse(text=expression));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients$fixed
  lm1.anova = anova(lm1)
  return(list(tabOut, lm1.anova))
}

Covariates$Subject = factor(Covariates$Subject)
Covariates$RegionBA9 = ifelse(Covariates$Region == "BA9", 1, 0)
Covariates$RegionBA38 = ifelse(Covariates$Region == "BA38", 1, 0)
Covariates$RegionBA39_40 = ifelse(Covariates$Region == "BA39_40", 1, 0)
Covariates$RegionBA44_45 = ifelse(Covariates$Region == "BA44_45", 1, 0)
Covariates$BrainBankNICHD = ifelse(Covariates$BrainBank == "BrainBankNICHD-BTB", 1, 0)

expression_lm
expression_model = "lme(score ~ Diagnosis + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score + BrainBankNICHD, random=~1|Subject, data = cur_data)" # BrainBankNICHD bothers with nlme
expression_model = "lme(score ~ Diagnosis + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + tssenrich.score, random=~1|Subject, data = cur_data)"
expression_model = "lme(score ~ Diagnosis + Age + PMI + RegionBA9 + RegionBA38 + RegionBA39_40 + RegionBA44_45 + FRiP, random=~1|Subject, data = cur_data)"

n = length(unlist(str_split(expression_model, "\\+")))
p = magnitude = matrix(nrow = nrow(df_Scores), ncol = n) 

for (i in 1:nrow(df_Scores)) {
  cur_data = df_Scores[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("score")
  #colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
  lm1.out <- try(runlme(cur_data,expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of ATAC peak', i, rownames(df_Scores)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}
# All motifs show error (Singularity in backsolve at level 0, block 1)
# Removing BB solved the Singularity Issue.
length(which(is.na(p[,1]))) # 0

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(df_Scores)

expression_model
name1 = gsub("lme\\(score ~ ", "", expression_model)
name2 = gsub(", random=~1\\|Subject, data = cur_data\\)", "", name1)
(name3 = gsub(" \\+ ", "_", name2))

save(p, magnitude, df_Scores, Covariates, file = paste0("p_magnitude_lme",name3,".rda"))

pdf(paste0("HistPval_lme",name3,".pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear model: TF binding score ~ ", name2, " + 1|Subject"))
dev.off()
# base: Diagnosis_Age_PMI_RegionBA9_RegionBA38_RegionBA39_40_RegionBA44_45_tssenrich
# base: sharp p<0.05 but abnormal historgram: low between p = 0.02-0.18, high at p=0.2
# new: Diagnosis_Age_PMI_RegionBA9_RegionBA38_RegionBA39_40_RegionBA44_45_FRiP
# new: great, very sharp p<0.05 and normal-looking histogram except for a little peak at p=0.4

## Check if there is any significant hits (fdr < 0.05)
fdr = p.adjust(p[,1], method = "fdr")
range(fdr, na.rm = T)
length(which(fdr < 0.05))
length(which(fdr < 0.1))
length(which(fdr < 0.2))
length(which(p[,1] < 0.05))

# new: Diagnosis_Age_PMI_RegionBA9_RegionBA38_RegionBA39_40_RegionBA44_45_FRiP
# new: FDR range from 0.13 to 1. 0 fdr < 0.1, 308 fdr < 0.2, 227 p < 0.05

## Plot volcano plot
rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/ind_BINDetect_Jaspar2022/")
lnames = load("p_magnitude_lmeDiagnosis_Age_PMI_RegionBA9_RegionBA38_RegionBA39_40_RegionBA44_45_FRiP.rda")

colnames(p); colnames(magnitude)
df = data.frame(TF = rownames(p), magnitude = magnitude[,1], p = p[,1], sig =  -log10(p[,1]))

plot(df$magnitude, df$sig, pch = 19, col = alpha("black", 0.5))
abline(v = 0.0068, h = 1.7)
abline(v = -0.0044, h = 1) # to include CTCF

df$col = "NS"
#df$col[which(df$magnitude > 0 & df$sig > 1.9)] = "Up_inASD"
df$col[which(df$magnitude < -0.0044 & df$sig > 1)] = "Down_inASD"
df$col[which(df$magnitude > 0.0068 & df$sig > 1)] = "Up_inASD"
#df$col[which(df$magnitude < -0.0054 & df$sig > 1)] = "Down_inASD"

df$TF_shortname = sapply(df$TF, function(x) unlist(str_split(x, "-"))[2])
df$TFlabel = ifelse(df$col != "NS", df$TF_shortname, "")

library(ggplot2)
library(ggrepel)
pdf("Volcano_BindingScore_lme_xxBA9xxFRiP.pdf", height = 6, width = 10)
df %>%
  ggplot(aes(x = magnitude, y = sig, col = col, label = TFlabel)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, size = 3) +
  theme_bw() +
  xlab("Differential TF binding score (ASD vs. CTL)") +
  ylab("-log10(p-value)") +
  xlim(-0.017, 0.017) +
  scale_color_manual(values = c("NS" = alpha("black", 0.5), "Up_inASD" = alpha("#F8766D", 0.8), "Down_inASD" = alpha("dodgerblue", 0.8)))
dev.off()
# Observation:
# Many overlap with the old Volcano plot without linear model.

unique(df$TFlabel[df$col == "Up_inASD"]) 
# [1] "FOS"          "FOSL2"        "MEF2C"        "DBP"          "JDP2"         "MEF2B"   
# [7] "OLIG2"        "MEF2D"        "BHLHE23"      "OLIG1"        "OLIG3"        "NFE2"     
# [13] "TEF"          "FOS::JUN"     "FOSL1::JUN"   "FOSL2::JUN"   "JUN::JUNB"    "FOS::JUNB"
# [19] "FOSB::JUNB"   "FOSL1::JUNB"  "FOSL2::JUNB"  "FOS::JUND"    "FOSL1::JUND"  "FOSL2::JUND" 
# [25] "BACH2"        "BATF::JUN"    "NEUROG1"      "BATF"         "BATF3"        "FOSL1"   
# [31] "HLF"          "JUNB"         "JUND"         "MEF2A"        "MAFG::NFE2L1" "NFIL3"   
# [37] "BNC2"         "ZNF211"       "BHLHA15"      "BHLHE22"      "BACH1" 
unique(df$TFlabel[df$col == "Down_inASD"]) 
# [1] "CTCF"   "PLAG1"  "ZBTB33" "HINFP"  "ZIC1"   "GLIS2"  "ZIC4"   "TFAP2A" "TFAP2B" "TFAP2C"
# [11] "TFDP1"  "KLF15"  "KLF3"   "MAZ"    "PLAGL2" "TFAP2E" "ZFP57"  "ZIC5"   "ZNF682" "ZNF454"
# [21] "ZNF610" "ZNF93" (KZNF)  "ZBTB14" "ZNF417" "CTCFL"  "YY2"    "PATZ1"  "ZNF320" "ZNF530" "ZNF574"
# [31] "ZNF669" "ZNF692" "THAP1" 
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1457042/ for a comprehensive catalog of human KZNFs. (Huntley 2006 Genome Res) - can't open the catalog database link
# See https://www.biorxiv.org/content/biorxiv/early/2022/08/15/2022.02.03.479047.full.pdf for TOBIAS findings with ZNF in plants
# No ZNF here is DEG in Jill's

## Now check the footprints to filter for the valid ones. 
TFmotif_for_PlotAggregate = paste0("Hsapiens-", rownames(df)[which(df$col != "NS")])
TFmotif_for_PlotAggregate = gsub("\\(", "", TFmotif_for_PlotAggregate)
TFmotif_for_PlotAggregate = gsub("\\)", "", TFmotif_for_PlotAggregate)
TFmotif_for_PlotAggregate = gsub("::", "", TFmotif_for_PlotAggregate)

write_lines(TFmotif_for_PlotAggregate, file = "../PlotAggregate/DiffTFbindingscoreFromlmeVolcano_Jaspar2022TFmotif_for_PlotAggregate.txt")
# Run script 10_3_5_PlotAggregate_filterTFoi_Jaspar2022.sh

# Invalid ones below:
invalid = c("BHLHE22", "CTCFL", "GLIS2", "HINFP", "KLF3", "KLF15", "MAZ", "PATZ1", "PLAG1", "PLAGL2", "TFAP2*", "TFDP1", "THAP1", "ZBTB14", "ZFP57", "ZIC*", "ZNF*") # But ZNF211 is fine
# refer to Bentson 2020 Nat Commun Fig2

df$TF_footprint = ifelse(grepl(paste0(invalid, collapse = "|"), df$TF_shortname), "bad", "good")
df$TF_footprint[df$TF_shortname == "ZNF211"] = "good"
table(df$TF_footprint) # good 601, bad 90

# Re-plot
df$TFlabel = ifelse(df$col != "NS" & df$TF_footprint != "bad", df$TF_shortname, "")

pdf("Volcano_BindingScore_lme_xxBA9xxFRiP_LabelFootprintPass.pdf", height = 6, width = 10)
df %>%
  ggplot(aes(x = magnitude, y = sig, col = col, label = TFlabel)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, size = 3) +
  theme_bw() +
  xlab("Difference in mean TF binding score (ASD vs. CTL)") +
  ylab("-log10(p-value)") +
  xlim(-0.017, 0.017) +
  scale_color_manual(values = c("NS" = alpha("black", 0.5), "Up_inASD" = alpha("#F8766D", 0.8), "Down_inASD" = alpha("dodgerblue", 0.8)))
dev.off()

DiffBS_lm_Volcano = df
rownames(DiffBS_lm_Volcano) = paste0("Hsapiens-", rownames(DiffBS_lm_Volcano))
tmp = str_split_fixed(rownames(DiffBS_lm_Volcano), "-",3)
DiffBS_lm_Volcano$TFmotif = tmp[,3]

# Attach TF family
BINDetect_res = read.table("../BINDetect_Jaspar2022/bindetect_results.txt", header = T)
DiffBS_lm_Volcano$cluster = BINDetect_res$cluster[match(rownames(DiffBS_lm_Volcano),BINDetect_res$motif_id)]

save(DiffBS_lm_Volcano, file = "DiffBS_lm_Volcano.rda")

