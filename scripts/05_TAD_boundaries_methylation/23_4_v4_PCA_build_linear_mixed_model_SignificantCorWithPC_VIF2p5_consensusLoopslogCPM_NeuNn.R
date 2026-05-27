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

Covariates = datMetaSeq[match(colnames(NeuNn_logCPM), rownames(datMetaSeq)),] # 16 obs. of 19 variables

colnames(Covariates)
apply(Covariates[,13:15], 1, sum) # all sum up to 1.

# Format biological and technical Covariates
Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "Dup15q"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal"))
Covariates$PMI[is.na(Covariates$PMI)] = mean(Covariates$PMI, na.rm = T)
# Combining BA9 and BA10 to BA9_10
Covariates$Region[which(Covariates$Region %in% c("BA9", "BA10"))] = "BA9_10"
Covariates$Region = factor(Covariates$Region, levels = c("BA7", "BA9_10", "BA4_6", "BA24"))
Covariates = Covariates[,c("Subject", "Diagnosis", "Age", "PMI", "Sex", "BrainBank", "Cortex", "Region", "Valid", "Duplicate", "Trans", "cis_shortRange", "cis_longRange", "ReadDepth")]
colnames(Covariates) = gsub("_", ".", colnames(Covariates))

## Identify top PCs
norm <- t(scale(t(NeuNn_logCPM),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:11)]) ## Top 11 PCs explain >80% of total variance.
topPC <- PC$rotation[,1:11]

### Correlation ###
## Correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates)[-1]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
colnames(Cor)[ncol(topPC) + c(1, 4:9)] # after combining BA9_10
idx_spearman = ncol(topPC) + c(1, 4:9) # after combining BA9_10
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
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # none

# Bonferroni correction of the p-val is too stringent, use fdr and see how many covariates are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisDup15q! 
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.2)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisDup15q! Just PMI.
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.5)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Finally there is DiagnosisDup15q! 

## Corrplot (the 2nd plot of each denotes fdr<0.2 by *)
library(corrplot)
sig_level = 0.2
# change the old corrplot name as xxx_FrontalAsBase
pdf(paste0("Corrplot_topPCofNeuNnlogCPM_Covariates_includingTechRepsAndOutliers_UseBA7asBase.pdf"), height = 35, width = 45) 
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
# [1] "DiagnosisDup15q"    "Age"                "PMI"               "BrainBankNICHD-BTB"
# [5] "RegionBA9_10"       "RegionBA4_6"        "Valid"             "Duplicate"      
# [9] "Trans"              "cis.shortRange"     "cis.longRange"     

## First include all candidate covariates that cor > 0.4 with PC1-PC3 and that are not significantly correlated with each other.
Cov_pool[sapply(Cov_pool, function(x) any(abs(Cor[1:3,x]) > 0.4))]
# [1] "DiagnosisDup15q" "Age"             "PMI"             "RegionBA9_10"   
# [5] "Duplicate"       "Trans"           "cis.shortRange"  "cis.longRange" 

Cov_pool[!sapply(Cov_pool, function(x) any(abs(Cor[1:3,x]) > 0.4))]
# These ones later: "BrainBankNICHD-BTB" "RegionBA4_6"        "Valid"    

# Prioritize biological before technical covariate.
test_row = which(rownames(Cor_sig) == "DiagnosisDup15q")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include DiagnosisDup15q

test_row = which(rownames(Cor_sig) == "Age")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # RegionBA9_10
test_row = which(rownames(Cor_sig) == "RegionBA9_10")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Age

sum(abs(Cor[1:3,which(colnames(Cor) == "Age")])) # 1.052687
sum(abs(Cor[1:3,which(colnames(Cor) == "RegionBA9_10")])) # 0.8336506
# Include Age and exempt RegionBA9_10

test_row = which(rownames(Cor_sig) == "PMI")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include PMI (remember that some values are imputed with the mean of other samples)

# Now technical covariates:
test_row = which(rownames(Cor_sig) == "Duplicate")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Duplicate

test_row = which(rownames(Cor_sig) == "Trans")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # cis.longRange
test_row = which(rownames(Cor_sig) == "cis.longRange")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Trans

sum(abs(Cor[1:3,which(colnames(Cor) == "Trans")])) # 1.410719
sum(abs(Cor[1:3,which(colnames(Cor) == "cis.longRange")])) # 1.47002
sum(abs(Cor[1:2,which(colnames(Cor) == "Trans")])) # 0.848734
sum(abs(Cor[1:2,which(colnames(Cor) == "cis.longRange")])) # 0.9449162
# Include cis.longRange

test_row = which(rownames(Cor_sig) == "cis.shortRange")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include cis.longRange

# Now come to the less important but significant ones
test_row = which(rownames(Cor_sig) == "BrainBankNICHD-BTB")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include BrainBankNICHD-BTB

test_row = which(rownames(Cor_sig) == "RegionBA4_6")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include RegionBA4_6

test_row = which(rownames(Cor_sig) == "Valid")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include Valid

## Check VIF
library(olsrr)
# Base model
expression_lm = "lm(logCPM ~ DiagnosisDup15q + Age + PMI + Duplicate + cis.longRange + cis.shortRange, data = cur_data)"
# Age and PMI have VIF = 2.2, VIF of the other two < 2. Fine.

# Will BrainBankNICHD-BTB cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + Age + PMI + Duplicate + cis.longRange + cis.shortRange + BrainBankNICHD, data = cur_data)"
# VIF Dup15q and PMI > 2.5, VIF NICHD > 4.
table(Covariates[,c("Diagnosis", "BrainBank")])
# Diagnosis Harvard-ATP NICHD-BTB
# CTL              3         5
# Dup15q           8         0
# Do not include NICHD.

# Will BA4_6 cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + Age + PMI + Duplicate + cis.longRange + cis.shortRange + RegionBA4_6, data = cur_data)"
# PMI and cis.longRange VIF > 2.9, do not include BA4_6

# Will Valid cause a problem?
expression_lm = "lm(logCPM ~ DiagnosisDup15q + Age + PMI + Duplicate + cis.longRange + cis.shortRange + Valid, data = cur_data)"
# VIF of PMI, Duplicate, Valid > 2.5, Age VIF > 3.5. Do not include Valid.

i=1
cur_data = NeuNn_logCPM[i,]
cur_data = as.data.frame(cbind(t(cur_data), mod_mat))
colnames(cur_data)[1] = c("logCPM")
colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
colnames(cur_data)[which(colnames(cur_data) == "BrainBankHarvard-ATP")] = "BrainBankHarvard"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign

## Final model: 
# [Using BA7 as base, combined BA9_10]
expression_lm = "lm(logCPM ~ DiagnosisDup15q + Age + PMI + Duplicate + cis.longRange + cis.shortRange, data = cur_data)"

## Run linear regression first (later lmm). 
runlm <- function(cur_data,expression_lm) {
  lm1 <- eval(parse(text = expression_lm));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients
  return(tabOut)
}

#Covariates$RegionBA4_6 = ifelse(Covariates$Region == "BA4_6", 1, 0)
#Covariates$RegionBA9_10 = ifelse(Covariates$Region == "BA9_10", 1, 0)

expression_lm = gsub("DiagnosisDup15q", "Diagnosis", expression_lm)
expression_lm

n = length(unlist(str_split(expression_lm, "\\+")))
p = magnitude = matrix(nrow = nrow(NeuNn_logCPM), ncol = n) 

for (i in 1:nrow(NeuNn_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = NeuNn_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  tabOut <- try(runlm(cur_data,expression_lm),silent=F)
  magnitude[i,] <- tabOut[-1,1]
  p[i,] <- tabOut[-1,4]
}
# Use BA7 as base and combine BA9_10
# 2022/09/09 12:31 PM (use cis.longRange + cis.shortRange)

tabOut
colnames(p) = colnames(magnitude) = rownames(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(NeuNn_logCPM)

e1 = unlist(str_split(expression_lm, "~ "))[2]
(e2 = unlist(str_split(e1, ","))[1])
e3 = gsub("_", ".", e2)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_NeuNn_lm",e3,"_includingTechRepsAndOutliers.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {npos = nchar(e2)}

pdf(paste0("HistPval_NeuNn_lm",e3,"_includingTechRepsAndOutliers.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the lm with cis.longRange + cis.shortRange: Looks nice. -> Run lme.

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
#expression_lm = "lm(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 + RegionBA9_10 + Valid, data = cur_data)"
#expression_lm = "lm(logCPM ~ Diagnosis + Sex + Trans + RegionBA4_6 + Valid, data = cur_data)"

expression_model = gsub("lm","lme",expression_lm)
expression_model = gsub(",",", random=~1|Subject,", expression_model)
expression_model

n = length(unlist(str_split(expression_model, "\\+")))
p = magnitude = matrix(nrow = nrow(NeuNn_logCPM), ncol = n) 

for (i in 1:nrow(NeuNn_logCPM)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = NeuNn_logCPM[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("logCPM")
  lm1.out <- try(runlme(cur_data,expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of ATAC peak', i, rownames(NeuNn_logCPM)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}

# Use BA7 as base and combining BA9_10
# 2022/09/09 12:33-12:37 PM (Use cis.longRange + cis.shortRange)
# Some Error messages.
length(which(is.na(p[,1]))) 
# With the lme using cis.shortRange + cis.shortRange: 91 NeuNn loops show error (Error in lme.formula(logRPKM ~ Diagnosis + Age + PMI + Duplicate + cis.longRange +  : nlminb problem, convergence error code = 1; message = singular convergence (7)).

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(NeuNn_logCPM)

e1 = unlist(str_split(expression_model, "~ "))[2]
(e2 = unlist(str_split(e1, ", data"))[1])
e3 = unlist(str_split(e2, ","))[1]
e3 = gsub("_", ".", e3)
(e3 = gsub(" \\+ ", "", e3))
save(p, magnitude, file = paste0("p_magnitude_NeuNn_lme",e3,"_includingTechRepsAndOutliers_UseBA7asBase.rda"))

pos = str_locate_all(pattern = "\\+", e2)[[1]][,1]
pos2 = str_locate_all(pattern = ",", e2)[[1]][,1]
if (length(pos) >= 5) {npos = pos[5]} else {npos = pos2}

pdf(paste0("HistPval_NeuNn_lme",e3,"_includingTechRepsAndOutliers_UseBA7asBase.pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = paste0("Linear mixed model: ", substr(e2,1,npos), "\n", substr(e2,npos+1,nchar(e2))))
dev.off()
# Use the lm with cis.longRange + cis.shortRange: Looks nice! Even better than lm(). This is similar to ATAC results.

range(p[,1], na.rm = T) # 2e-7 to 1
length(which(p[,1] < 0.05)) # 4663 nominally significant loops, more than the 2877 loops from limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
fdr = p.adjust(p[,1])
range(fdr, na.rm = T) # 0.09 to 1. Good
length(which(fdr < 0.1)) # 7 significant loops, fewer than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
length(which(fdr < 0.2)) # 7 significant loops, fewer than the limma results with previous linear model (using highest adjR2) from 22_11_Limma_onEachCelltype.R
# Consider using limma-voom, see if I would get better results.

## Check correlation between topPCs of regressed logCPM with covariates
colnames(magnitude)
colnames(mod_mat)
colnames(magnitude)[which(colnames(magnitude) == "BrainBankNICHD")] = "BrainBankNICHD-BTB"
colnames(magnitude)[which(colnames(magnitude) == "BrainBankHarvard")] = "BrainBankHarvard-ATP"
mod_used = mod_mat[,colnames(magnitude)[-1]]

logCPM.Regressed = NeuNn_logCPM - as.matrix(magnitude[,-1]) %*% t(mod_used)
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
idx_spearman = ncol(topPC) + c(1, 4:9) # after combining BA9_10
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

pdf(paste0("Corrplot_TopPCofNeuNnRegressedLogCPM_Covariates_lme",e3,"_includingTechRepsAndOutliers.pdf"), height = 30, width = 40)
corrplot(Cor[1:ncol(topPC),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
dev.off()
# Use lme with cis.longRange + cis.shortRange: NICHD and Trans inevitably correlate with PC1 as they are associated with Dup15q. ReadDepth and Valid also show association with the topPCs, but Valid violates the VIF<2.5 criteria.

## Final model for NeuNn:
expression_model
# "lme(logCPM ~ Diagnosis + Age + PMI + Duplicate + cis.longRange + cis.shortRange, random=~1|Subject, data = cur_data)"

## Any correlation with Bulk ASD vs CTL?
lnames = load("p_magnitude_NeuNn_lmeDiagnosisAgePMIDuplicatecis.longRangecis.shortRange_includingTechRepsAndOutliers_UseBA7asBase.rda") # p magnitude
magnitude_NeuNn = magnitude
lnames = load("p_magnitude_Bulk_lmeDiagnosisBatchValidTransBrainBankHarvardSexRegionBA39.40RegionBA44.45_includingTechRepsAndOutliers.rda") # p magnitude

tmp = intersect(rownames(magnitude_NeuNn), rownames(magnitude)) # 22699 loops shared
df = as.data.frame(cbind(magnitude_NeuNn[tmp,1], magnitude[tmp,1]))
colnames(df) = c("logFC_NeuNn", "logFC_ASD")

plot(df$logFC_ASD, df$logFC_NeuNn)
# There is correlation! Great!
library(ggplot2)
model = lm(logFC_NeuNn ~ 0 + logFC_ASD, df)
model_res = summary(model)
lm_coef = round(model_res$coefficients[1,1], digits = 2) # 0.3
r_coef = round(sqrt(model_res$r.squared), digits = 2) # 0.17
p_coef = round(model_res$coefficients[1,4], digits = 2) # 0 p<2.2e-16

df_location = as.data.frame(str_split_fixed(rownames(df), "_",3))
df_location$V2 = as.integer(df_location$V2)
df_location$V3 = as.integer(df_location$V3)
df$PLsWithinDup15qRegion = ifelse(df_location$V1 == "chr15" & df_location$V3 < 33e6, "red", "black")

pdf("Cor_logFC_ASDbulkvsDup15qNeuNn_lme.pdf", height = 5, width = 8)
plot(df$logFC_ASD, df$logFC_NeuNn, 
     pch = 19, 
     col = alpha(df$PLsWithinDup15qRegion, 0.5),
     xlab = "logFC bulk ASD vs. CTL", ylab = "logFC NeuNn Dup15q vs. CTL", main = "Promoter loops")
abline(coef = c(0, lm_coef), col = "red")
#text(0.5, 1, paste("R =", r_coef, "\np =", p_coef, "\nslope =", lm_coef), col = "red")
text(0.3, 1, paste("R =", r_coef, "\np < 2e-16\nslope =", lm_coef), col = "red")
legend("topleft", legend=c("Loops within the Dup15q region",
                           "Loops outside of the Dup15q region"),
       fill=alpha(c("red", "black"), 0.7))
dev.off()
