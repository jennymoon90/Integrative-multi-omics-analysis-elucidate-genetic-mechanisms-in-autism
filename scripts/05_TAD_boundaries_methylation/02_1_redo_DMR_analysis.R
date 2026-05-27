rm(list = ls())
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/DNA_methylation/results/02_redo_DMG_analysis")

# ---- Procedures ----
### (1) Load beta-values of the filtered HM450k data
### (2) Check whether data has been normalized using the wateRmelon package
### (3) Decide how to group the probes
### (4) Group probes within the same CpG island together, and the rest as themselves
### (5) Get the top PCs
### (6) Determine the linear mixed model
### (7) Run lme
### (8) Set FDR < 0.05 as threshold for DMGs
# --------------------

### (1) Load beta-values of the filtered HM450k data

lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/32_Gokul_H3K27Ac_datasets/GoogleDrive/Differential_DNA_methylation/ASD_Methylation.RData") # datMeth_anno, datMeta_FC/TC/CB, datMeth_FC/TC/CB
nrow(datMeth_CB) # 417460, the same as Wong et al 2019 Methods section stated: The final dataset concists of a total of 417460 probes from 76 PFC samples (n=36 iASD patients, 7 dup15q patients, 33 controls), 77 TC samples (33 iASD patients, 6 dup15q patients, 38 controls)

range(datMeta_FC$Age) # 2 to 71
hist(datMeta_FC$Age, breaks = 50)

### (2) Check whether data has been normalized using the wateRmelon package

# BiocManager::install('wateRmelon')
# library(wateRmelon)
# packageVersion("wateRmelon") # ‘2.4.0’
# Refer to https://bioconductor.org/packages/release/bioc/vignettes/wateRmelon/inst/doc/wateRmelon.html#installation

## Test whether the data is already normalized.
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/32_Gokul_H3K27Ac_datasets/GoogleDrive/Differential_DNA_methylation/Methylation_wholeGene_forRedoDMG.RData") # datMeth_CTX_wholeGene, CpG_list, datMeta_CTX
all(rownames(datMeta_CTX) %in% c(rownames(datMeta_FC), rownames(datMeta_TC))) # T
Gokul_value = datMeth_CTX_wholeGene[1,1] # 0.2634332
Gokul_sample = colnames(datMeth_CTX_wholeGene)[1] # AN17515_ba9
Gokul_gene = rownames(datMeth_CTX_wholeGene)[1] # ENSG00000103479
Gokul_probes = unlist(CpG_list[Gokul_gene])
Data_value = mean(datMeth_FC[which(rownames(datMeth_FC) %in% Gokul_probes) , which(colnames(datMeth_FC) == Gokul_sample)]) # 0.2634332
all.equal(Data_value, Gokul_value) # T. 

## Conclusion: the dataset has already been normalized by wateRmelon, and can be used directly for differential analysis.

rm(list = ls())
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/32_Gokul_H3K27Ac_datasets/GoogleDrive/Differential_DNA_methylation/ASD_Methylation.RData") # datMeth_anno, datMeta_FC/TC/CB, datMeth_FC/TC/CB

all(colnames(datMeth_FC) == rownames(datMeta_FC)) # T
all(colnames(datMeth_TC) == rownames(datMeta_TC)) # T
all(rownames(datMeth_FC) == rownames(datMeth_TC)) # T

datMeth_CTX = cbind(datMeth_FC, datMeth_TC)
datMeta_CTX = rbind(datMeta_FC, datMeta_TC)
rm(list = setdiff(ls(), c("datMeth_CTX", "datMeta_CTX", "datMeth_anno")))
save(list = ls(), file = "02_datMeta_datMeth_anno_CTX.rda")

### (3) Decide how to group the probes
rm(list = ls())
lnames = load("02_datMeta_datMeth_anno_CTX.rda")

# The number of probes are too huge. May need to cluster them to regions based on distance. -> Read Loyfer 2023 Nature "A DNA methylation atlas of normal human cell types"

## Load HM450k probe annotations
HM450k = read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/39_WongCC_DNAmethylationASDsupplementary/HM450k_Illumina_csv/humanmethylation450_15017482_v1-2.xlsx", skip = 7) # the xlsx file was download from Illumina https://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
colnames(HM450k) # IlmnID, CHR, MAPINFO store the cg_ID, chr, pos.
HM450k = HM450k[,c(1,12:13,25:26, 2:11,14:24,27:ncol(HM450k))]
colnames(HM450k)[1:3] = c("Probe", "chr", "pos")
HM450k$chr = paste0("chr", HM450k$chr)

unique(HM450k$chr[match(rownames(datMeth_CTX), HM450k$Probe)]) # There is chr1-Y.
# I will ignore chrX and chrY for now.

# Shall we group by CpG island?
length(unique(HM450k$UCSC_CpG_Islands_Name[!is.na(HM450k$UCSC_CpG_Islands_Name)])) # 27176, this is an ideal scale for differential analysis
length(HM450k$UCSC_CpG_Islands_Name[!is.na(HM450k$UCSC_CpG_Islands_Name)]) # 309465
# Many CpG probes belong to the same CpG island!
length(HM450k$UCSC_CpG_Islands_Name[is.na(HM450k$UCSC_CpG_Islands_Name)]) # 176963 CpGs not related to CpG islands

length(unique(HM450k$UCSC_CpG_Islands_Name[!is.na(HM450k$UCSC_CpG_Islands_Name) & ! HM450k$chr %in% c("chrX", "chrY", "chrMULTI")])) # 26328, this is an ideal scale for differential analysis. Increased by 2 compared to using the wrong genome build
length(HM450k$UCSC_CpG_Islands_Name[!is.na(HM450k$UCSC_CpG_Islands_Name) & ! HM450k$chr %in% c("chrX", "chrY", "chrMULTI")]) # 300885
# Many CpG probes belong to the same CpG island!
length(HM450k$UCSC_CpG_Islands_Name[is.na(HM450k$UCSC_CpG_Islands_Name) & ! HM450k$chr %in% c("chrX", "chrY", "chrMULTI")]) # 173895 CpGs not related to CpG islands. Increased by 32 compared to using the wrong genome build
# will remove a few thousand CpG probes, but does not change the scale.

library(stringr)
HM450k_CpGisland = as.data.frame(str_split_fixed(HM450k$UCSC_CpG_Islands_Name, ":|-", 3))
HM450k_CpGisland$V2 = as.integer(HM450k_CpGisland$V2); HM450k_CpGisland$V3 = as.integer(HM450k_CpGisland$V3)
colnames(HM450k_CpGisland) = c("chr", "start", "end")

HM450k_CpGisland$width = HM450k_CpGisland$end - HM450k_CpGisland$start
hist(HM450k_CpGisland$width) # some can be > 40kb, but most < 2kb
mean(HM450k_CpGisland$width, na.rm = T) # 935 bp
median(HM450k_CpGisland$width, na.rm = T) # 683 bp
# Refer to https://www.youtube.com/watch?v=2964vuECi-A. CpG islands are defined as stretches of > 200bp (Illumina suggest > 500bp) with GC% > 50% (Illumina suggest > 55%). 40% of gene promoters contain CpG islands! CpG shelves are ~ 4kb from islands and CpG shores are ~ 2kb from islands

# Any DMP in UCSC_CpG_Islands?
lnames = load("../01_use_DMRresultsFromWong/01_1_DMRs.rda") # DMRs_expanded
DMRs_expanded$UCSC_CpG_Islands_Name = HM450k$UCSC_CpG_Islands_Name[match(DMRs_expanded$Probe, HM450k$Probe)]
table(is.na(DMRs_expanded$UCSC_CpG_Islands_Name)) # 90 F, 67 T.
length(unique(DMRs_expanded$UCSC_CpG_Islands_Name)) # 89.
DMRs_expanded$UCSC_CpG_Islands_Name[which(duplicated(DMRs_expanded$UCSC_CpG_Islands_Name) & ! is.na(DMRs_expanded$UCSC_CpG_Islands_Name))] # "chr4:100870377-100871994" "chr6:30523775-30525189".
# The 2 DMPs in chr4:100870377-100871994, one is S_Shelf, the other is S_Shore, and their methylation levels change in different direction.
# The 2 DMPs in chr6:30523775-30525189, both are in Islands and both show hypermethylation in ASD.
DMRs_expanded$Relation_to_UCSC_CpG_Island = HM450k$Relation_to_UCSC_CpG_Island[match(DMRs_expanded$Probe, HM450k$Probe)]
# Also check chr20:10652573-10655611, which has only one DMP in the island

#tmp = HM450k[which(HM450k$UCSC_CpG_Islands_Name == "chr2:85765695-85766983"),]
tmp = HM450k[which(HM450k$UCSC_CpG_Islands_Name == "chr4:100870377-100871994"),]
tmp = HM450k[which(HM450k$UCSC_CpG_Islands_Name == "chr6:30523775-30525189"),]
tmp = HM450k[which(HM450k$UCSC_CpG_Islands_Name == "chr20:10652573-10655611"),]
cin = gsub(":", "_", unique(tmp$UCSC_CpG_Islands_Name))
cin = gsub("-", "_", cin)

library(tidyverse) 
tmp = tmp %>%
  arrange(pos)
tmp$Relation_to_UCSC_CpG_Island # N_Shelf - N_Share - Island - S_Shore - S_Shelf

# Generate heatmap of normalized beta-values to decide how to group CpGs into regions
library(pheatmap)

datMeta_CTX = datMeta_CTX %>%
  arrange(Phenotype, Age, BrainRegion)
datMeth_CTX = datMeth_CTX[, match(rownames(datMeta_CTX), colnames(datMeth_CTX))]

tmp_data = datMeth_CTX[match(tmp$Probe, rownames(datMeth_CTX)),] # 2 NA probes
#tmp_data = tmp_data[complete.cases(tmp_data),]
annotation_row = as.data.frame(tmp[match(rownames(tmp_data), tmp$Probe), "Relation_to_UCSC_CpG_Island"])
rownames(annotation_row) = rownames(tmp_data)
annotation_col = as.data.frame(data_frame(Diagnosis = datMeta_CTX$Phenotype))
rownames(annotation_col) = rownames(datMeta_CTX)

pdf(paste0("03_CpG_island_shore_shelf_datMeth_across_samples_",cin,".pdf"), width = 15, height = 10) # width = 10, height = 5 for the first 2; width = 15, height = 10 for the last one on chr6
pheatmap(tmp_data, cluster_rows = F, cluster_cols = F, 
         annotation_row = annotation_row,
         annotation_col = annotation_col) # default scale = "none"
dev.off()
# chr2:85765695-85766983: 1/4 of N_shore gets similar methylation level as N_shelf, while the other 3/4 looks similar to Island. The Island itself can be divided into two regions based on the normalized beta values.
# chr4:100870377-100871994: 1/2 of N_shore looks similar to Island. The island show bands of high and low methylated probes. It is obvious that in some samples the entire column gets hypomethylated.
# chr20:10652573-10655611: The island methylation level is pretty uniform with some bands on the edge.

# Probably better to set distance range and mean methylation level difference to group the probes
tmp$mean = apply(tmp_data, 1, mean)
tmp = tmp[,c(1:5, ncol(tmp), 6:(ncol(tmp)-1))]

pdf(paste0("03_CpG_island_shore_shelf_datMeth_mean_and_distance",cin,".pdf"), width = 8, height = 4)
#plot(x = rownames(tmp), y = tmp$mean)
plot(x = tmp$pos, y = tmp$mean, main = unique(tmp$UCSC_CpG_Islands_Name))
#abline(v = c(85618200, 85619000, 85619160, 85619450, 85619905), col = "red") # chr2:85765695-85766983:
dev.off()
# chr2:85765695-85766983: It would make sense to group CpGs within 300 bp distance and mean beta difference <= 0.05 into CpG clusters.

# Generate heatmap of beta-values scaled by row to check whether these CpG methylation show consistent changes across samples
tmp_data_scaled = t(scale(t(tmp_data)))
pdf(paste0("03_CpG_island_shore_shelf_datMeth_scaled_across_samples",cin,".pdf"), width = 15, height = 10) # width = 10, height = 5 for the first 2; width = 15, height = 10 for the last one on chr6
pheatmap(tmp_data_scaled, cluster_rows = F, cluster_cols = F, 
         annotation_row = annotation_row,
         annotation_col = annotation_col)
dev.off()
# chr2:85765695-85766983: In general it's a mess, but there is some consistency in the changes. eg. the red bar (one sample) at N_short and Island top half; the blue bar (one sample) at Island top half.
# chr4:100870377-100871994: After arrange the CTLs in the middle and ASDs to the left and scale the changes, it's obvious that the two single DMPs show a difference. The red bars are just sparse samples. The entire Island (do not consider the Shore and Shelf) seem to be hypomethylated in ASD as well if we take the average across samples!!!
# chr6:30523775-30525189: It is obvious that the entire Island show hypomethylation in ASD and Dup15q. 
# chr20:10652573-10655611: DMP is cg05799256, but the entire island look hypomethylated in ASD.

## Conclusion, we can potentially average the normalized beta-values from the entire CpG island. For shelves and shores, use them as single probes. This shall give us more power in DMR detection.

## Prepare dataset by deleting chrX, chrY and chrMULTI probes
HM450k_final = HM450k[which(HM450k$Probe %in% rownames(datMeth_CTX) & ! HM450k$chr %in% c("chrX", "chrY", "chrMULTI")),]
HM450k_final = HM450k_final %>% 
  arrange(chr, pos)
datMeth_CTX = datMeth_CTX[match(HM450k_final$Probe, rownames(datMeth_CTX)),]

save(HM450k_final, datMeth_CTX, datMeta_CTX, datMeth_anno, file = "03_Normalized_Beta_anno_datMeta_arrangedByDiagnosisAgeBA.rda")

### (4) Group probes within the same CpG island together, and the rest as themselves
rm(list = ls())
lnames = load("03_Normalized_Beta_anno_datMeta_arrangedByDiagnosisAgeBA.rda")

colnames(HM450k_final)
tmp = HM450k_final[,c("Probe", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]

Probe_ByItself = HM450k_final[which(HM450k_final$Relation_to_UCSC_CpG_Island != "Island" | is.na(HM450k_final$Relation_to_UCSC_CpG_Island)),] # 274231 rows
Probe_Island = HM450k_final[which(HM450k_final$Relation_to_UCSC_CpG_Island == "Island"),]
nrow(Probe_ByItself) + nrow(Probe_Island) == nrow(HM450k_final) # T
Probe_Island_grouped = Probe_Island[, c("Probe", "UCSC_CpG_Islands_Name")] %>%
  group_by(UCSC_CpG_Islands_Name) %>%
  mutate(Probes = paste(Probe, collapse = ",")) %>%
  ungroup()
Probe_Island_grouped = unique(Probe_Island_grouped[,c("Probes", "UCSC_CpG_Islands_Name")]) # 24719 rows

# Decide to limit my differential analysis to CpG islands.
save(Probe_Island_grouped, Probe_Island, Probe_ByItself, file = "04_Probes_inCpGislands_orByItself.rda")

rm(list = ls())
lnames = load("03_Normalized_Beta_anno_datMeta_arrangedByDiagnosisAgeBA.rda") # HM450k_final, datMeth_CTX, datMeta_CTX, datMeth_anno
lnames = load("04_Probes_inCpGislands_orByItself.rda") # Probe_Island_grouped, Probe_Island, Probe_ByItself
rm(Probe_ByItself)

## Limit to CpG islands
all(HM450k_final$Probe == rownames(datMeth_CTX))
idx = which(HM450k_final$Probe %in% Probe_Island$Probe)
HM450k_final = HM450k_final[idx,]
datMeth_CTX = datMeth_CTX[idx,]

datMeth_CTX_island = as.data.frame(matrix(nrow = nrow(Probe_Island_grouped), ncol = ncol(datMeth_CTX)))
colnames(datMeth_CTX_island) = colnames(datMeth_CTX)
rownames(datMeth_CTX_island) = Probe_Island_grouped$UCSC_CpG_Islands_Name

for (i in 1:nrow(datMeth_CTX_island)) {
  if (i %% 1e3 == 0) {print(i)}
  probes = unlist(str_split(Probe_Island_grouped$Probes[i], ","))
  idx_probes = which(rownames(datMeth_CTX) %in% probes)
  datMeth_CTX_island[i,] = apply(datMeth_CTX[idx_probes,], 2, mean)
}

save(datMeth_CTX_island, datMeta_CTX, Probe_Island_grouped, HM450k_final, file = "04_Averaged_normalizedBeta_withinCpGisland_anno_datMeta_arrangedByDiagnosisAgeBA.rda")

### (5) Get the top PCs

## PCA analysis
rm(list = ls())
lnames = load("04_Averaged_normalizedBeta_withinCpGisland_anno_datMeta_arrangedByDiagnosisAgeBA.rda")

norm <- t(scale(t(datMeth_CTX_island),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:61)]) ## Top 61 PCs explain 80% of total variance.
topPC <- PC$rotation[,1:61]

save(PC, varexp, topPC, file = "05_topPCs_of_AveragedNormalizedBetaInIslands.rda")

### (6) Determine the linear mixed model

## Correlation of covariates with the top PCs
# Remove columns with the same values in datMeta_CTX
idx = apply(datMeta_CTX, 2, function(x) which(all(x == x[1])))
idx2 = which(idx == 1)
datMeta_CTX = datMeta_CTX[,-idx2]

# Select meaningful columns
colnames(datMeta_CTX)
# unique(datMeta_CTX$Sentrix_Full) # 153
# unique(datMeta_CTX$DilPlateName) # 4
# unique(datMeta_CTX$DilPlateGrid) # 84
#idx = which(datMeta_CTX$newBrainBankID != datMeta_CTX$BrainBankID)
#tmp = datMeta_CTX[idx,] # "/" becomes "_" in newBrainBankID
# unique(datMeta_CTX$newBrainBankID) # 88 individuals
#idx = which(datMeta_CTX$BrainRegion != datMeta_CTX$BrainRegion_update)
#tmp = datMeta_CTX[idx,] # BA41-42-22 becomes BA41 in BrainRegion_update
# unique(datMeta_CTX$BrainRegion_update) # BA41, BA9
datMeta_CTX$ba9 = ifelse(datMeta_CTX$BrainRegion_update == "BA9", 1, 0)
# sexCode: M=0, F=1
#unique(datMeta_CTX$BrainCentre.M) # "Harvard" "NICHD"   "Oxford" 
# unique(datMeta_CTX$predictedGender) # "male"   "female" "Unsure"
# unique(datMeta_CTX$Sex) # M, F
# idx = which(datMeta_CTX$predictedGender == "male"); all(datMeta_CTX$Sex[idx] == "M") # T
# idx = which(datMeta_CTX$predictedGender == "female"); all(datMeta_CTX$Sex[idx] == "F") # T

Covariates = datMeta_CTX[,c(4,5,9,56,13:16,22,23,26:28,47,49,51:53,55)]
Covariates$Diagnosis = ifelse(Covariates$Phenotype == "CTL", 0, 1)
Covariates$Phenotype = factor(Covariates$Phenotype, levels = c("CTL", "ASD", "dup15q"))
unique(Covariates$Batch) # 1 2
Covariates$batch = ifelse(Covariates$Batch == 1, 0, 1)

mod_mat_expr = paste(c(colnames(Covariates)[-3]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
idx_spearman = ncol(topPC) + c(1:5, 9:11, 13:15)
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
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # DilPlateNamePlate, ba9, No.of.probes.with.pval.0.05(.1), Call.rate, Phenotypedup15q, sexCode, Age, CET, DNAmAge, meanXchromosome, Diagnosis
# PhenotypeASD does not significantly correlate with any topPCs, although Phenotypedup15q and Diagnosis do. Fine.

# Bonferroni correction of the p-val is stringent, use fdr and see how many covaraites are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Now see PhenotypeASD! 16 covariates total

save(topPC, Cor, Covariates, file = "06_Corrplot_topPCofAveragedNormalizedBetaInIslands_Covariates.rda")

## Corrplot (the 2nd plot of each denotes fdr<0.05 by *)
library(corrplot)
SigLev = 0.05

pdf(paste0("06_Corrplot_topPCofAveraged_Covariates.pdf"), height = 35, width = 45) 
# In-between covariates:
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor), (ncol(topPC) + 1):ncol(Cor)], method = "ellipse", tl.pos = "lt", tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor), (ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[(ncol(topPC) + 1):ncol(Cor), (ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = SigLev, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
# TopPCs with covariates:
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], method="ellipse", tl.pos = "lt", tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2) 
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = SigLev, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
dev.off()

# Observations:
# 1. In-between covariates: 
# a. Batch and DilPlateNamePlate ATPDIL16APR14 1 show cor = 1, 
# b. No.of.probes.with.pval.0.05, No.of.probes.with.pval.0.05.1, Call.rate show cor = 1 or -1
# c. Age and DNAmAge show cor = 0.97
# 2. TopPCs with covariates: 
# a. PC1 cor with CET (cor = 0.7)
# b. PC5 cor with ba9 (cor = 0.8)
# c. Some other covariates show cor with topPCs.
# d. dup15q show highest cor with PC10 (cor = 0.4); Diagnosis show highest cor with PC4 (cor = 0.3)

Cov_pool
# [1] "DilPlateNamePlate ATP Autism Brain_Dil 2" "DilPlateNamePlate ATP Autism Brain_Dil 3"
# [3] "ba9"                                      "No.of.probes.with.pval.0.05"      
# [5] "No.of.probes.with.pval.0.05.1"            "Call.rate"                        
# [7] "PhenotypeASD"                             "Phenotypedup15q"                  
# [9] "sexCode"                                  "Age"                              
# [11] "nichd"                                    "oxford"                          
# [13] "CET"                                      "DNAmAge"                         
# [15] "meanXchromosome"                          "Diagnosis"

test_row = which(rownames(Cor_sig) == "PhenotypeASD")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include PhenotypeASD

test_row = which(rownames(Cor_sig) == "Phenotypedup15q")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include Phenotypedup15q

# Consider group ASD and Phenotypedup15q later
test_row = which(rownames(Cor_sig) == "Diagnosis")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PhenotypeASD, Age, oxford, DNAmAge
# Try to include Diagnosis and Age in the same linear model

test_row = which(rownames(Cor_sig) == "CET")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include CET

test_row = which(rownames(Cor_sig) == "Age")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # DNAmAge, harvard, oxford, Diagnosis
# Try to include Diagnosis and Age in the same linear model

test_row = which(rownames(Cor_sig) == "Call.rate")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # No.of.probes.with.pval.0.05, No.of.probes.with.pval.0.05.1
# Include Call.rate

test_row = which(rownames(Cor_sig) == "ba9")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none
# Include ba9

test_row = which(rownames(Cor_sig) == "sexCode")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # meanXchromosome
# Include sexCode

test_row = which(rownames(Cor_sig) == "nichd")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # harvard
# Need to test after putting age in the model.

test_row = which(rownames(Cor_sig) == "DilPlateNamePlate ATP Autism Brain_Dil 2")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # DilPlateNamePlate ATP Autism Brain_Dil 3
# May include Batch and DilPlateNamePlate ATP Autism Brain_Dil 2

## Check VIF
library(olsrr)

expression_lm = "lm(datMeth ~ Diagnosis + CET + Age + ba9 + sexCode + batch + DilPlate2 + Call.rate, data = cur_data)"
# all VIF < 2, Great!

expression_lm = "lm(datMeth ~ Diagnosis + CET + Age + ba9 + sexCode + batch + DilPlate2 + Call.rate + nichd, data = cur_data)"
# all VIF < 2, Great!

i=1
cur_data = datMeth_CTX_island[i,]
cur_data = as.data.frame(cbind(t(cur_data), mod_mat))
colnames(cur_data)[1] = c("datMeth")
colnames(cur_data)[which(colnames(cur_data) == "DilPlateNamePlate ATP Autism Brain_Dil 2")] = "DilPlate2"

fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign

## Final model: 
expression_lm = "lm(datMeth ~ Diagnosis + CET + Age + ba9 + sexCode + batch + DilPlate2 + Call.rate + nichd, data = cur_data)"

save(expression_lm, mod_mat, Covariates, datMeth_CTX_island, HM450k_final, Probe_Island_grouped, file = "06_final_lm.rda")

### (7) Run lme
rm(list = ls())
lnames = load("06_final_lm.rda")

library(nlme)
runlme <- function(thisdat,expression) {
  lm1 <- eval(parse(text=expression));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients$fixed
  lm1.anova = anova(lm1)
  return(list(tabOut, lm1.anova))
}

expression_lm
#expression_model = "lme(logRPKM ~ Diagnosis + Age + Batch + RegionBA38 + RegionBA44_45 + tssenrich.score + FRiP, random=~1|Subject, data = cur_data)"
expression_model = gsub("lm", "lme",expression_lm)
expression_model = gsub(",", ", random=~1|Subject,",expression_model)
Covariates$Subject = Covariates$newBrainBankID
Covariates$DilPlate2 = ifelse(Covariates$DilPlateName == "Plate ATP Autism Brain_Dil 2", 1, 0)

n = length(unlist(str_split(expression_lm, "\\+")))
p = magnitude = matrix(nrow = nrow(datMeth_CTX_island), ncol = n) 

for (i in 1:nrow(datMeth_CTX_island)) {
  if (i %% 5000 == 0) {print(paste0("Done ", i, "th CpG island"))}
  cur_data = datMeth_CTX_island[i,]
  cur_data = as.data.frame(cbind(t(cur_data), Covariates))
  colnames(cur_data)[1] = c("datMeth")
  lm1.out <- try(runlme(cur_data, expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of CpG island', i, rownames(datMeth_CTX_island)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}
# 2023/10/02 5:59-6:04PM
length(which(is.na(p[,1]))) 
# 8 CpG islands show error (Error in lme.formula(datMeth ~ Diagnosis + CET + Age + ba9 + sexCode +  : nlminb problem, convergence error code = 1 message = singular convergence (7)). This is Normal.

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(datMeth_CTX_island)
expr_name = gsub("lm\\(datMeth ~ ", "", expression_lm)
expr_name = gsub(", data = cur_data\\)", "", expr_name)
expr_name = gsub(" \\+ ", "_", expr_name)

save(p, magnitude, expression_model, file = paste0("07_p_magnitude_lme",expr_name,".rda"))

expression_model
# [1] "lme(datMeth ~ Diagnosis + CET + Age + ba9 + sexCode + batch + DilPlate2 + Call.rate + nichd, random=~1|Subject, data = cur_data)"

pdf(paste0("07_HistPval_lme",expr_name,".pdf"), height = 5, width = 8)
hist(p[,1], breaks = 50, xlab = "Diagnosis p-value", main = "Linear model: Diagnosis + CET + Age + ba9 + sexCode\n+ batch + DilPlate2 + Call.rate + nichd, random=~1|Subject")
dev.off()
# Looks amazing, very sharp p<0.05

### (8) Set FDR < 0.05 as threshold for differentially methylated CpG islands (DMRs)
rm(list = ls())
lnames = load("07_p_magnitude_lmeDiagnosis_CET_Age_ba9_sexCode_batch_DilPlate2_Call.rate_nichd.rda") # p, magnitude, expression_model
lnames = load("06_final_lm.rda") # expression_lm, mod_mat, Covariates, datMeth_CTX_island, HM450k_final, Probe_Island_grouped
rm(mod_mat)

p = as.data.frame(p)
fdr = p.adjust(p$Diagnosis, method = "fdr")
range(fdr, na.rm = T) # 0 to 1
length(which(fdr < 0.05)) # 8899 DMRs, Great!

magnitude = as.data.frame(magnitude)

DMR = as.data.frame(cbind(magnitude$Diagnosis, p$Diagnosis, fdr))
rownames(DMR) = rownames(p)
colnames(DMR) = c("magnitude", "pval", "fdr")

Probe_Island_grouped$chr = unname(sapply(Probe_Island_grouped$UCSC_CpG_Islands_Name, function(x) unlist(str_split(x, ":"))[1]))
tmp = sapply(Probe_Island_grouped$UCSC_CpG_Islands_Name, function(x) unlist(str_split(x, ":"))[2])
Probe_Island_grouped$start = as.integer(sapply(tmp, function(x) unlist(str_split(x, "-"))[1]))
Probe_Island_grouped$end = as.integer(sapply(tmp, function(x) unlist(str_split(x, "-"))[2]))
Probe_Island_grouped$Island_length = Probe_Island_grouped$end - Probe_Island_grouped$start

hist(Probe_Island_grouped$Island_length, breaks = 200) # Mostly small a few hundred bp

save(DMR, HM450k_final, Probe_Island_grouped, file = "08_DMR_fdr005.rda")



