### Get top PCs that explain total 80% variance in the un-normalized logCPM, identify biological and technical covariates that significantly correlate with these top PCs.

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/23_differential_loop_using_package/mimickATAC_lmBuild/")

### Load data ### 
## Load un-normalized loop logCPM data
lnames = load("../../25_nlme_DiffLoopAnalysis/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_includingOutlier.rda") # "Bulk_logCPM" "NeuNp_logCPM" "NeuNn_logCPM" from 25_1_Loop_logCPM_includingTechReps.R

## Load Covariates
lnames = load("../../25_nlme_DiffLoopAnalysis/datMetaSeq.rda") # "datMetaSeq" "Sample_rm" from 25_1_Loop_logCPM_includingTechReps.R

Covariates = datMetaSeq[match(colnames(NeuNp_logCPM), rownames(datMetaSeq)),] # 16 obs. of 19 variables

colnames(Covariates)
apply(Covariates[,13:15], 1, sum) # all sum up to 1.

# Format biological and technical Covariates
Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "Dup15q"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal"))
Covariates$PMI[is.na(Covariates$PMI)] = mean(Covariates$PMI, na.rm = T)
Covariates = Covariates[,c("Subject", "Diagnosis", "Age", "PMI", "Sex", "BrainBank", "Cortex", "Region", "Valid", "Duplicate", "Trans", "cis_shortRange", "cis_longRange", "ReadDepth")]
colnames(Covariates) = gsub("_", ".", colnames(Covariates))

## Identify top PCs
norm <- t(scale(t(NeuNp_logCPM),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:11)]) ## Top 11 PCs explain >80% of total variance.
topPC <- PC$rotation[,1:11]

### Correlation ###
# At first, I include cis.interaction%, but then using fdr<0.2 I was not able to include as many covariates as possible. So now I remove cis.interaction% to reduce the multiple testing burden.
# Covariates = Covariates[,- which(startsWith(colnames(Covariates), "cis"))]
# Deleting the cis.interactions% from Covariates get the same covariates pool. Do not delete cis.interactions%

# What about combining BA9 and BA10?
Covariates$Region[which(Covariates$Region %in% c("BA9", "BA10"))] = "BA9_10"
Covariates$Region = factor(Covariates$Region, levels = c("BA7", "BA9_10", "BA4_6", "BA24"))

## Correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates)[-1]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
colnames(Cor)[ncol(topPC) + c(1, 4:9)] # 4:10 original before combining BA9_10
idx_spearman = ncol(topPC) + c(1, 4:9) # 4:10 original before combining BA9_10
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
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # Trans, cis.longRange

# Bonferroni correction of the p-val is too stringent, use fdr and see how many covariates are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisDup15q! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.1)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # There is DiagnosisDup15q! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.2)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Originally same variables as fdr<0.1. To be consistent with ASD, just use fdr<0.2 as cutoff. After using BA7 as the base, now adding RegionBA9_10. Great.
# Deleting the cis.interactions% from Covariates get the same covariates pool. Do not delete cis.interactions%

## Corrplot (the 2nd plot of each denotes fdr<0.2 by *)
library(corrplot)
sig_level = 0.2
# change the old corrplot name as xxx_FrontalAsBase
pdf(paste0("Corrplot_topPCofNeuNplogCPM_Covariates_includingTechRepsAndOutliers_UseBA7asBase.pdf"), height = 35, width = 45) 
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
## I know that Cortex/Region/BrainBank are correlated
unique(Covariates[,c("Cortex", "Region", "BrainBank")])
# Frontal    BA9   NICHD-BTB
# Frontal  BA4_6   NICHD-BTB
# Frontal   BA24   NICHD-BTB
# Frontal   BA24 Harvard-ATP
# Frontal  BA4_6 Harvard-ATP
# Frontal    BA9 Harvard-ATP
# Parietal   BA7 Harvard-ATP
# Frontal   BA10 Harvard-ATP

Cov_pool
# [1] "DiagnosisDup15q" "SexF"            "RegionBA4_6"     "Valid"           "Trans" 
# [6] "cis.shortRange"  "cis.longRange" 
# Now adding "RegionBA9_10"

## First include all candidate covariates that cor > 0.4 with PC1-PC3 and that are not significantly correlated with each other.
Cov_pool[sapply(Cov_pool, function(x) any(abs(Cor[1:3,x]) > 0.4))]
# [1] "DiagnosisDup15q" "SexF" "Trans" "cis.shortRange" "cis.longRange"  

Cov_pool[!sapply(Cov_pool, function(x) any(abs(Cor[1:3,x]) > 0.4))]
# These ones later: "RegionBA4_6" "Valid". Now adding "RegionBA9_10"    

# Prioritize biological before technical covariate.
test_row = which(rownames(Cor_sig) == "DiagnosisDup15q")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include DiagnosisDup15q

test_row = which(rownames(Cor_sig) == "SexF")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include SexF

# Now technical covariates:
test_row = which(rownames(Cor_sig) == "Trans")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PC2, cis.shortRange, cis.longRange
test_row = which(rownames(Cor_sig) == "cis.shortRange")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Trans, cis.longRange
test_row = which(rownames(Cor_sig) == "cis.shortRange")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Trans, cis.shortRange

sum(abs(Cor[1:3,which(colnames(Cor) == "Trans")])) # 1.182531
sum(abs(Cor[1:3,which(colnames(Cor) == "cis.shortRange")])) # 1.139213
sum(abs(Cor[1:3,which(colnames(Cor) == "cis.longRange")])) # 1.20132
# Include cis.longRange
sum(abs(Cor[1:2,which(colnames(Cor) == "Trans")])) # 1.045653
sum(abs(Cor[1:2,which(colnames(Cor) == "cis.shortRange")])) # 1.015922
sum(abs(Cor[1:2,which(colnames(Cor) == "cis.longRange")])) # 0.9199614
# Include Trans

# Now come to the less important but significant ones
test_row = which(rownames(Cor_sig) == "RegionBA4_6")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include RegionBA4_6
test_row = which(rownames(Cor_sig) == "RegionBA9_10")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Age
# Include RegionBA9_10

test_row = which(rownames(Cor_sig) == "Valid")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Valid

## Check VIF
library(olsrr)
# Base model
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + cis.longRange, data = cur_data)"
# cis.longRange VIF = 2.1, VIF of the other two < 2. Fine.

# Will BA4_6 and BA9_10 cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + cis.longRange + RegionBA4_6, data = cur_data)"
# cis.longRange VIF = 2.3, VIF of the other two < 2. Fine.
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + cis.longRange + RegionBA9_10, data = cur_data)"
# cis.longRange VIF > 3. Cannot include RegionBA9_10
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA9_10, data = cur_data)"
# Trans VIF = 2.2, VIF of the other two < 2. Fine.
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA4_6 + RegionBA9_10, data = cur_data)"
# Trans VIF = 2.2, VIF of the other two < 2. Fine to include both BA regions.

# Will Valid cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + Valid, data = cur_data)"
# Trans VIF = 2.4, VIF of the other two < 2. Fine.

# Will combining BA regions and Valid cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA9_10 + RegionBA4_6 + Valid, data = cur_data)"
# Trans VIF > 2.8. Cannot include all.
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + cis.longRange + RegionBA9_10 + RegionBA4_6 + Valid, data = cur_data)"
# cis.longRange VIF > 3.5. Cannot include all. Spare Valid and just include the two BA regions.
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA4_6 + Valid, data = cur_data)"
# Trans VIF = 2.46. Fine.
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA9_10 + Valid, data = cur_data)"
# Trans VIF > 2.8. Nope.

i=1
cur_data = NeuNp_logCPM[i,]
cur_data = as.data.frame(cbind(t(cur_data), mod_mat))
colnames(cur_data)[1] = c("logCPM")
colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
colnames(cur_data)[which(colnames(cur_data) == "BrainBankHarvard-ATP")] = "BrainBankHarvard"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign

## Final model: 
# [old without using BA7 as base or combining BA9_10]
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + cis.longRange + RegionBA4_6 + Valid, data = cur_data)"
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA4_6 + Valid, data = cur_data)"
# New
expression_lm = "lm(logCPM ~ DiagnosisDup15q + SexF + Trans + RegionBA4_6 + RegionBA9_10, data = cur_data)"


## Run linear regression first (later lmm). 
runlm <- function(cur_data,expression_lm) {
  lm1 <- eval(parse(text = expression_lm));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients
  return(tabOut)
}

Covariates$RegionBA4_6 = ifelse(Covariates$Region == "BA4_6", 1, 0)
Covariates$RegionBA9_10 = ifelse(Covariates$Region == "BA9_10", 1, 0)

expression_lm = gsub("DiagnosisDup15q", "Diagnosis", expression_lm)
expression_lm = gsub("SexF", "Sex", expression_lm)
expression_lm

n = length(unlist(str_split(expression_lm, "\\+")))
p = magnitude = matrix(nrow = nrow(NeuNp_logCPM), ncol = n) 

for (i in 1:nrow(NeuNp_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = NeuNp_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  tabOut <- try(runlm(cur_data,expression_lm),silent=F)
  magnitude[i,] <- tabOut[-1,1]
  p[i,] <- tabOut[-1,4]
}
# 2022/09/08 9:56 PM (use cis.longRange)
# 2022/09/08 10:20 PM (use Trans)
# New setting BA7 as base and combine BA9_10
# 2022/09/08 11:01 PM (use Trans and BA9_10)

tabOut
colnames(p) = colnames(magnitude) = rownames(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(NeuNp_logCPM)

e1 = unlist(str_split(expression_lm, "~ "))[2]
(e2 = unlist(str_split(e1, ","))[1])
e3 = gsub("_", ".", e2)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_NeuNp_lm",e3,"_includingTechRepsAndOutliers.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {npos = nchar(e2)}

pdf(paste0("HistPval_NeuNp_lm",e3,"_includingTechRepsAndOutliers.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the lm with cis.longRange: Only a little enriched at p<0.05 -> Run lme.
# Use the lm with Trans: Similar. Only a little enriched at p<0.05 -> Run lme.
# Use the lm with Trans and BA9_10: Only a little enriched at p<0.05 -> Run lme.

## Run linear regression - using linear mixed model, including random effect of subjects
library(nlme)
runlme <- function(thisdat,expression) {
  lm1 <- eval(parse(text=expression));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients$fixed
  lm1.anova = anova(lm1)
  return(list(tabOut, lm1.anova))
}

# Modifications:
expression_lm = "lm(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 + RegionBA9_10 + Valid, data = cur_data)"
expression_lm = "lm(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 + Valid, data = cur_data)"

expression_model = gsub("lm","lme",expression_lm)
expression_model = gsub(",",", random=~1|Subject,", expression_model)
expression_model

n = length(unlist(str_split(expression_model, "\\+")))
p = magnitude = matrix(nrow = nrow(NeuNp_logCPM), ncol = n) 

for (i in 1:nrow(NeuNp_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = NeuNp_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  lm1.out <- try(runlme(cur_data,expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of ATAC peak', i, rownames(NeuNp_logCPM)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}
# 2022/09/08 9:58-10:02 PM (Use cis.longRange)
# 2022/09/08 10:21-10:25 PM (Use Trans)
# New setting BA7 as base and combining BA9_10
# 2022/09/08 11:02-11:06 PM (Use Trans and BA4_6 and BA9_10)
# 2022/09/08 11:28-11:32 PM (Use Trans and BA4_6 and BA9_10 and Valid)
# 2022/09/08 11:39-11:43 PM (Use Trans and BA4_6 and Valid)
# 2022/09/08 11:?-11:? PM (Use Trans and BA9_10 and Valid)
# Some Error messages.
length(which(is.na(p[,1]))) 
# With the lme using cis.shortRange: 135 NeuNp loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Sex + cis.longRange + RegionBA4_6 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# With the lme using Trans: 134 NeuNp loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Sex + Trans + RegionBA4_6 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# New setting BA7 as base and combine BA9_10
# With the lme using Trans and BA4_6 and BA9_10: 79 NeuNp loops show error (Error in lme.formula(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# With the lme using Trans and BA4_6 and BA9_10 and Valid: 155 NeuNp loops show error (Error in lme.formula(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).
# With the lme using Trans and BA4_6 and Valid: 134 NeuNp loops show error (Error in lme.formula(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(NeuNp_logCPM)

e1 = unlist(str_split(expression_model, "~ "))[2]
(e2 = unlist(str_split(e1, ", data"))[1])
e3 = unlist(str_split(e2, ","))[1]
e3 = gsub("_", ".", e3)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_NeuNp_lme",e3,"_includingTechRepsAndOutliers_UseBA7asBase.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
pos2 = str_locate_all(pattern = ",", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {npos = pos2}

pdf(paste0("HistPval_NeuNp_lme",e3,"_includingTechRepsAndOutliers_UseBA7asBase.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear mixed model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the lm with cis.longRange: Looks nice! Better than lm(). This is similar to ATAC results.
# Use the lm with Trans: Similar to using cis.longRange. Looks nice! Better than lm(). This is similar to ATAC results.
# New setting BA7 as base
# Use the lm with Trans and BA4_6 and BA9_10: Similar to the model wo BA9_10. Looks nice! Better than lm(). This is similar to ATAC results.
# Use the lm with Trans and BA4_6 and BA9_10 and Valid: Similar to the model wo BA9_10. Looks nice! Better than lm(). This is similar to ATAC results.
# Use the lm with Trans and BA4_6 and Valid: No change from using BA10 as base.

range(p[,1], na.rm = T) # 6.6e-7 (cis.longRange) or 7.9e-7 (Trans) or 1.8e-6 (Trans and BA9_10) or 1.8e-6 (Trans and BA9_10 and Valid) to 1
length(which(p[,1] < 0.05)) # 2815 (cis.longRange) or 2835 (Trans) or 2734 (Trans and BA9_10) or 2851 (Trans and BA9_10 and Valid) nominally significant loops, more than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
fdr = p.adjust(p[,1])
range(fdr, na.rm = T) # 0.025 (cis.longRange) or 0.03 (Trans) or 0.06 (Trans and BA9_10) or 0.07 (Trans and BA9_10 and Valid) to 1. Good
length(which(fdr < 0.05)) # 3 (cis.longRange) or 2 (Trans) or 0 (Trans and BA9_10 and or wo Valid) significant loops, fewer than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
length(which(fdr < 0.1)) # 6 or 2 (Trans and BA9_10 and or wo Valid) significant loops, fewer than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
length(which(fdr < 0.2)) # 10 (cis.longRange) or 9 (Trans) or 3 (Trans and BA9_10) or 7 (Trans and BA9_10 and Valid) significant loops, fewer than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
# Consider using limma-voom, see if I would get better results.
# Including BA9_10 gets even worse results. Better to use Valid wo BA9_10

## Check correlation between topPCs of regressed logCPM with covariates
colnames(magnitude)
colnames(mod_mat)
colnames(magnitude)[which(colnames(magnitude) == "BrainBankNICHD")] = "BrainBankNICHD-BTB"
colnames(magnitude)[which(colnames(magnitude) == "BrainBankHarvard")] = "BrainBankHarvard-ATP"
mod_used = mod_mat[,colnames(magnitude)[-1]]

logCPM.Regressed = NeuNp_logCPM - as.matrix(magnitude[,-1]) %*% t(mod_used)
idx = which(apply(logCPM.Regressed, 1, function(x) all(is.na(x))))
if (length(idx) > 0) {logCPM.Regressed = logCPM.Regressed[-idx,]}

norm <- t(scale(t(logCPM.Regressed),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:8)]) # Top 8 PCs explain >80% of variance.
topPC <- PC$rotation[,1:8] ## these first 15 explain ~80% of the variance

mod_mat_withPC = cbind(topPC, mod_mat)
Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
idx_spearman = ncol(topPC) + c(1, 4:9) # 4:10 before combining BA9_10
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

pdf(paste0("Corrplot_TopPCofNeuNpRegressedLogCPM_Covariates_lme",e3,"_includingTechRepsAndOutliers.pdf"), height = 30, width = 40)
corrplot(Cor[1:ncol(topPC),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
dev.off()
# Use lme with cis.longRange: Age, PMI, NICHD, BA9, Duplicate, Trans, ReadDepth all remain correlated with the 2 topPCs.
# Use lme with Trans: Age, PMI, NICHD, BA24, BA9, Duplicate, Trans, ReadDepth all remain correlated with the 2 topPCs.
# Need to refine the fdr cutoff when selecting covariates. Or maybe remove cis.interaction to decrease multiple comparisons.
# Use lme with Trans and BA9_10: NICHD, Valid, ReadDepth all remain correlated with the 3 topPCs. Try to include Valid, although VIF > 2.5
table(Covariates[,c("Diagnosis", "BrainBank")])
# Diagnosis Harvard-ATP NICHD-BTB
# CTL              3         6
# Dup15q           7         0
# BrainBank is associated with Diagnosis
# I know that Trans is different b/w CTL and Dup15q as well.

## Any correlation with Bulk ASD vs CTL?
#lnames = load("p_magnitude_NeuNp_lmeDiagnosisSexTransRegionBA4.6RegionBA9.10_includingTechRepsAndOutliers.rda") # p magnitude
#lnames = load("p_magnitude_NeuNp_lmeDiagnosisSexTransRegionBA4.6RegionBA9.10Valid_includingTechRepsAndOutliers.rda") # p magnitude
lnames = load("p_magnitude_NeuNp_lmeDiagnosisSexTransRegionBA4.6Valid_includingTechRepsAndOutliers_UseBA7asBase.rda") # p magnitude
magnitude_NeuNp = magnitude
lnames = load("p_magnitude_Bulk_lmeDiagnosisBatchValidTransBrainBankHarvardSexRegionBA39.40RegionBA44.45_includingTechRepsAndOutliers.rda") # p magnitude

tmp = intersect(rownames(magnitude_NeuNp), rownames(magnitude)) # 24185 loops shared
df = as.data.frame(cbind(magnitude_NeuNp[tmp,1], magnitude[tmp,1]))
colnames(df) = c("logFC_NeuNp", "logFC_ASD")

plot(df$logFC_ASD, df$logFC_NeuNp)
# Little correlation, using Trans and BA4_6 and BA9_10
# Little correlation, using Trans and BA4_6 and BA9_10 and Valid
# Little correlation, using Trans and BA4_6 and Valid

## Final model for NeuNp:
expression_model
# "lme(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 + Valid, random=~1|Subject, data = cur_data)"

