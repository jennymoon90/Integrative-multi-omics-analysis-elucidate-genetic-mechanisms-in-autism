# Because I found Hi-C linked ATAC peaks show significantly higher correlateion than random pairs, so to link distal ATAC to gene promoters, in this script I would: 1) Extract Jill_DEG promoter ATAC (TSS +- 2kb) 2) Get all ATAC within promoter Hi-C loops 3) Calculate ATAC-ATAC cor using RPKM.cqn 4) Build null distribution using 10e3 trans ATAC-ATAC cor 5) Calculate p-value, set fdr threshold 0.05

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_corrected20230227")

### 1) Extract Jill_DEG promoter ATAC (TSS -2kb to +100bp)
lnames = load("../../5_DiffATAC/CQN.rda") # "RPKM.cqn"   "cqn.fit"    "Covariates" "peak_info"  "counts_tbl"
rm(list = setdiff(ls(), "RPKM.cqn"))

library(readxl)
Jill_DEG = readxl::read_excel("~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

## Find gene TSS
# BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="grch37.ensembl.org") 
# (att <- listAttributes(mart))
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","ensembl_transcript_id",
             "transcript_start", "transcript_end", "strand","transcript_biotype")
geneAnno <- getBM(attributes = getinfo,filters = c("ensembl_gene_id"), values = Jill_DEG$ensembl_gene_id, mart = mart)
unique(geneAnno$strand) # -1 1
geneAnno$TSS = ifelse(geneAnno$strand == "1", geneAnno$transcript_start, geneAnno$transcript_end)
# geneAnno$TSS = as.integer(geneAnno$TSS) # do it if necessary
geneAnno$TSS_bin = floor(geneAnno$TSS/10e3) * 10e3 + 5e3
geneAnno$chromosome_name = paste0("chr", geneAnno$chromosome_name)
colnames(geneAnno)[3] = "chr"
geneAnno$strand = ifelse(geneAnno$strand == "-1", "-", "+")
geneAnno$start = ifelse(geneAnno$strand == "+", geneAnno$TSS - 2000, geneAnno$TSS - 100)
geneAnno$end = ifelse(geneAnno$strand == "+", geneAnno$TSS + 100, geneAnno$TSS + 2000)

library(GenomicRanges)
promoter_gr = makeGRangesFromDataFrame(geneAnno)
save(promoter_gr, geneAnno, file = "10_3_8_v5_promoter_gr.rda")

library(stringr)
ATAC_location = as.data.frame(str_split_fixed(rownames(RPKM.cqn), "_", 3))
colnames(ATAC_location) = c("chr", "start", "end")
ATAC_location$start = as.integer(ATAC_location$start); ATAC_location$end = as.integer(ATAC_location$end)
ATAC_location$ATACpeak = rownames(RPKM.cqn)
ATAC_gr = makeGRangesFromDataFrame(ATAC_location)

promATACs = as.data.frame(findOverlaps(promoter_gr, ATAC_gr))
promATACs$promoter_ATAC = ATAC_location$ATACpeak[promATACs$subjectHits]
promATACs$ensembl_gene_id = geneAnno$ensembl_gene_id[promATACs$queryHits]
promATACs$chr = geneAnno$chr[promATACs$queryHits]
promATACs$TSS_bin = geneAnno$TSS_bin[promATACs$queryHits]
promATACs = unique(promATACs[,-c(1:2)]) # 29527 promoter ATAC-GE pairs

save(promATACs, file = "promATAC_TSSminus100bptoplus2kb.rda")

### 2) Get all ATAC within consensus promoter Hi-C loops (union of Bulk, NeuNp, NeuNn)
lnames = load("../../../../../../Hi-C/2nd_batch_20samples/data_analysis/results/25_nlme_DiffLoopAnalysis/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_includingOutlier.rda") # "Bulk_logCPM"  "NeuNp_logCPM" "NeuNn_logCPM"

union_loops = unique(c(rownames(Bulk_logCPM), rownames(NeuNp_logCPM), rownames(NeuNn_logCPM)))
union_loops_df1 = as.data.frame(str_split_fixed(union_loops, "_", 3)) # 67092 loops
union_loops_df1$V2 = as.integer(union_loops_df1$V2); union_loops_df1$V3 = as.integer(union_loops_df1$V3)
union_loops_df1$Loop = union_loops

union_loops_df2 = union_loops_df1
colnames(union_loops_df1)[1:3] = c("chr", "TSS_bin", "Distal_bin")
colnames(union_loops_df2)[1:3] = c("chr", "Distal_bin", "TSS_bin")

promATAC_Loop1 = left_join(promATACs, union_loops_df1)
promATAC_Loop2 = left_join(promATACs, union_loops_df2)
promATAC_Loop1 = promATAC_Loop1[complete.cases(promATAC_Loop1),]
promATAC_Loop2 = promATAC_Loop2[complete.cases(promATAC_Loop2),]
promATAC_Loop = unique(rbind(promATAC_Loop1, promATAC_Loop2)) # 111856 rows

Distal_bin = unique(data_frame(chr = promATAC_Loop$chr, start = promATAC_Loop$Distal_bin - 5e3, end = promATAC_Loop$Distal_bin + 5e3)) # 32511
Distal_bin_gr = annoDF2GR(Distal_bin)
Distal_bin$distal_bin = paste0(Distal_bin$chr, "_", Distal_bin$start + 5e3)
promATAC_Loop$distal_bin = paste0(promATAC_Loop$chr, "_", promATAC_Loop$Distal_bin)

hits = as.data.frame(findOverlaps(Distal_bin_gr, ATAC_gr))
hits$distal_bin = Distal_bin$distal_bin[hits$queryHits]
hits$distal_ATAC = ATAC_location$ATACpeak[hits$subjectHits]
hits = unique(hits[,3:4])

promATAC_Loop = left_join(promATAC_Loop, hits)
promATAC_Loop = promATAC_Loop[complete.cases(promATAC_Loop),] # 154403 promoter-distal ATAC pairs in promoter Hi-C loops
all(promATAC_Loop$promoter_ATAC %in% promATACs$promoter_ATAC) # T

Cor_Gene_ATAC = promATAC_Loop
save(Cor_Gene_ATAC, file = "ATACpairsInUnionPL.rda")

### 3) Get all distal ATAC peaks within 30kb of promoter ATAC peaks
Cor_Gene_ATAC_within30kb = matrix(nrow = 0, ncol = 3) # ensembl_gene_id, promoter_ATAC, distal_ATAC, cor (omit tss column to reduce calc burden)

for (g in Jill_DEG$ensembl_gene_id) {
  print(g)
  #gene_symbol = Jill_DEG$external_gene_name[Jill_DEG$ensembl_gene_id == g]
  chr = paste0("chr", Jill_DEG$chromosome_name[Jill_DEG$ensembl_gene_id == g])
  
  gene_info = geneAnno[geneAnno$ensembl_gene_id == g,]
  
  tss = gene_info$TSS
  # Strand = gene_info$strand
  # tss_df = data_frame(chr = chr, start = tss - 2000, end = tss + 2000)
  # tss_gr = makeGRangesFromDataFrame(tss_df)
  # 
  # # promoter ATAC
  # hits = as.data.frame(findOverlaps(tss_gr, ATAC_gr))
  # promATACs = ATAC_location$ATACpeak[unique(hits$subjectHits)]
  promATACs_g = unique(promATACs$promoter_ATAC[promATACs$ensembl_gene_id == g])
  
  # TSS +- 30kb
  tss_gr = makeGRangesFromDataFrame(data_frame(chr = chr, start = tss - 30e3, end = tss + 30e3))
  
  # All possible distal ATAC
  hits = as.data.frame(findOverlaps(tss_gr, ATAC_gr))
  distATACs = ATAC_location$ATACpeak[unique(hits$subjectHits)]
  
  # add to data frame
  cur_df = data_frame(ensembl_gene_id = g, promoter_ATAC = rep(promATACs_g, length(distATACs)), distal_ATAC = rep(distATACs, each = length(promATACs_g))) # tss = tss[hits$queryHits]
  idx_rm = which(cur_df$promoter_ATAC == cur_df$distal_ATAC)
  cur_df = cur_df[-idx_rm, ]
  Cor_Gene_ATAC_within30kb = rbind(Cor_Gene_ATAC_within30kb, cur_df)
}
# 2022/12/06 1:09-1:36AM
# corrected the above script at 2023/02/27 10:05-10:19PM
dim(Cor_Gene_ATAC_within30kb) # 221205 3
save(Cor_Gene_ATAC_within30kb, file = "../v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_corrected20230227/ATACpairs_within30kb.rda")
all(Cor_Gene_ATAC_within30kb$promoter_ATAC %in% promATACs$promoter_ATAC) # T, previously F before the correction

# Combine ATAC pairs in PL and within 30kb
tmp = unique(Cor_Gene_ATAC[,c("ensembl_gene_id", "promoter_ATAC", "distal_ATAC")])
Cor_Gene_ATAC_within30kb = anti_join(Cor_Gene_ATAC_within30kb, tmp) # 190302 3

Cor_Gene_ATAC_within30kb$chr = unname(sapply(Cor_Gene_ATAC_within30kb$promoter_ATAC, function(x) unlist(str_split(x, "_"))[1]))
Cor_Gene_ATAC_within30kb$distal_bin = Cor_Gene_ATAC_within30kb$Distal_bin = Cor_Gene_ATAC_within30kb$TSS_bin = 0
Cor_Gene_ATAC_within30kb$Loop = "Within30kb" # the ones with Hi-C loop is labeled as Loop
Cor_Gene_ATAC_within30kb = Cor_Gene_ATAC_within30kb[,match(colnames(Cor_Gene_ATAC), colnames(Cor_Gene_ATAC_within30kb))]

Cor_Gene_ATAC = rbind(Cor_Gene_ATAC, Cor_Gene_ATAC_within30kb)
save(Cor_Gene_ATAC, file = "ATACpairs_PLandWithin30kb.rda")
dim(Cor_Gene_ATAC) # 344705 8

### 4) Calculate ATAC-ATAC cor using RPKM.cqn 
# Cor_Gene_ATAC$cor = 0
# for (i in 1:nrow(Cor_Gene_ATAC)) {
#   if (i %% 20e3 == 0) {print(i)}
#   p = Cor_Gene_ATAC$promoter_ATAC[i]
#   d = Cor_Gene_ATAC$distal_ATAC[i]
#   p_atac = RPKM.cqn[rownames(RPKM.cqn) == p,]
#   d_atac = RPKM.cqn[rownames(RPKM.cqn) == d,]
#   Cor_Gene_ATAC$cor[i] = cor(p_atac,d_atac)
# }
# save(Cor_Gene_ATAC, file = "Cor_Gene_ATAC_calcedCor.rda")
# 2022/12/06 01:38-01:55PM

# 2023/02/27 after correction:
lnames = load("../v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Cor_Gene_ATAC_calcedCor.rda") # Cor_Gene_ATAC
Cor_Gene_ATAC_old = Cor_Gene_ATAC
lnames = load("ATACpairs_PLandWithin30kb.rda") # Cor_Gene_ATAC
Cor_Gene_ATAC = left_join(Cor_Gene_ATAC, Cor_Gene_ATAC_old)
save(Cor_Gene_ATAC, file = "Cor_Gene_ATAC_calcedCor.rda")

pdf("Cor_Gene_PromDistATAC.pdf", width = 8, height = 5)
hist(Cor_Gene_ATAC$cor)
dev.off()
# Observation: It's a skewed distribution, skewed to the right.

mean(Cor_Gene_ATAC$cor, na.rm = T) # PL and within 30kb: 0.23, higher than all promoter-distal ATAC pairs within 500Mb (0.15). 
range(Cor_Gene_ATAC$cor, na.rm = T) # -0.90 to 0.98

### 5) Build null distribution and calculate p-values
# using 10e3 trans ATAC-ATAC cor for each promoter ATAC
uniq_promATAC = unique(Cor_Gene_ATAC$promoter_ATAC); length(uniq_promATAC) # 23768 promoter ATAC peaks

tmp = as.data.frame(table(ATAC_location$chr))
nrow(ATAC_location) - max(tmp$Freq) # 117001, at max can randomly draw 117k trans ATAC peaks.

Cor_Gene_ATAC$p_perPromATAC = Cor_Gene_ATAC$trans_mean_perPromATAC = NA

set.seed(137)
for (i in 1:length(uniq_promATAC)) {
  if (i %% 5e3 == 0) {print(i)}
  p = uniq_promATAC[i]
  chr = unlist(str_split(p, "_"))[1]
  ATAC_chr_pool = which(ATAC_location$chr != chr)
  idx_sample = sample(ATAC_chr_pool, 10e3) # 20e3 doesn't make much difference
  null_atac = ATAC_location$ATACpeak[idx_sample]
  
  null_RPKM = t(RPKM.cqn[rownames(RPKM.cqn) %in% null_atac,])
  p_RPKM = RPKM.cqn[rownames(RPKM.cqn) == p,]
  null_cor = cor(p_RPKM, null_RPKM)
  # hist(null_cor) # It's also skewed, and sometimes irregular shaped.
  
  m = mean(null_cor) # 0.04 for chrX_99890783_99892116

  idx = which(Cor_Gene_ATAC$promoter_ATAC == p)
  Cor_Gene_ATAC$trans_mean_perPromATAC[idx] = m
  Cor_Gene_ATAC$p_perPromATAC[idx] = sapply(Cor_Gene_ATAC$cor[idx], function(x) length(which(null_cor > x))/10e3) # one-tail test
}
# 2022/12/06 01:56-02:02PM
# 2023/02/27 10:30-10:35PM after correction

hist(Cor_Gene_ATAC$p_perPromATAC, breaks = 50) # wow, very enriched for p< 0.05
Cor_Gene_ATAC$fdr_perPromATAC = p.adjust(Cor_Gene_ATAC$p_perPromATAC, method = "fdr")
save(Cor_Gene_ATAC, file = "Cor_Gene_ATAC_calcedCorP.rda")

length(which(Cor_Gene_ATAC$fdr_perPromATAC < 0.05)) # PL + 30kb: 9130 pairs
length(which(Cor_Gene_ATAC$fdr_perPromATAC < 0.01)) # PL + 30kb: 1780
union_loops_df1$loop_distance = union_loops_df1$Distal_bin - union_loops_df1$TSS_bin
range(union_loops_df1$loop_distance) # 20kb to 500kb, lacking ATAC peaks within 30kb of promoters!

sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.05,] # 9130 ATAC-ATAC pairs

## 6) Overlap with known eQTL, Hi-C loop, or other public databases.
rm(list = ls())
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/GTEx_eQTL/processed/Union_eQTL_CortexAndBA9_v8_hg19.rda") # Variant_Location_hg19_gr Union_CortexAndBA9_hg19, generated by 10_3_8_LinkATACtoGene_v3RPKMcqnDistalPromoterATACcor.R

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_corrected20230227/")
lnames = load("Cor_Gene_ATAC_calcedCorP.rda") # Cor_Gene_ATAC

# Try different ATAC cor FDR threshold
sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.05,] # 9130
sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.1,]; range(sig_CorATAC$cor) # 25184; 0.336 - 0.98. 
sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.15,]; range(sig_CorATAC$cor) # 44365; 0.29 - 0.98.
sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.2,]; range(sig_CorATAC$cor) # 63990; 0.23 - 0.98.

sig_DistalATAC_location = as.data.frame(str_split_fixed(sig_CorATAC$distal_ATAC, "_",3))
colnames(sig_DistalATAC_location) = c("chr", "start", "end")
sig_DistalATAC_location$start = as.integer(sig_DistalATAC_location$start)
sig_DistalATAC_location$end = as.integer(sig_DistalATAC_location$end)
sig_DistalATAC_gr = annoDF2GR(sig_DistalATAC_location)

hits = as.data.frame(findOverlaps(sig_DistalATAC_gr, Variant_Location_hg19_gr))
hits$distal_ATAC = sig_CorATAC$distal_ATAC[hits$queryHits]
hits$variant_id = Union_CortexAndBA9_hg19$variant_id[hits$subjectHits]
hits$ensembl_gene_id = Union_CortexAndBA9_hg19$ensembl_gene_id[hits$subjectHits]

library(tidyverse)
common = inner_join(unique(hits[, c("distal_ATAC", "ensembl_gene_id")]), unique(sig_CorATAC[, c("distal_ATAC", "ensembl_gene_id")])) 
nrow(common) # 523 -> 1367 -> 2251 -> 3084 rows 
(a = nrow(unique(hits[, c("distal_ATAC", "ensembl_gene_id")]))) # 3825 -> 8173 -> 11900 -> 14671
(b = nrow(unique(sig_CorATAC[, c("distal_ATAC", "ensembl_gene_id")]))) # 8043 -> 21067 -> 35700 -> 49727
a/b # 47.6% -> 38.8% -> 33.3% -> 29.5%

## Is this a significant overlap? 

# a. all eQTLs that were paired with the nearest gene were removed
# the following part only needs to be done once:
# ----
# library(readxl)
# Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")
# library(biomaRt)
# mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
#                 dataset="hsapiens_gene_ensembl",
#                 host="https://grch37.ensembl.org") 
# # (att <- listAttributes(mart))
# getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","ensembl_transcript_id",
#              "transcript_start", "transcript_end", "strand","transcript_biotype")
# geneAnno <- getBM(attributes = getinfo,filters = c("ensembl_gene_id"), values = Jill_DEG$ensembl_gene_id, mart = mart)
# unique(geneAnno$strand) # -1 1
# geneAnno$TSS = ifelse(geneAnno$strand == "1", geneAnno$transcript_start, geneAnno$transcript_end)
# save(geneAnno, file = "10_3_8_v5_geneAnno_JillDEG.rda")
# ----
lnames = load("10_3_8_v5_geneAnno_JillDEG.rda") # geneAnno

TSS = data_frame(chr = paste0("chr", geneAnno$chromosome_name), start = geneAnno$TSS, end = geneAnno$TSS, ensembl_gene_id = geneAnno$ensembl_gene_id)
TSS_gr = annoDF2GR(TSS)

hits2 = as.data.frame(distanceToNearest(Variant_Location_hg19_gr, TSS_gr))
hits2$variant_id = Union_CortexAndBA9_hg19$variant_id[hits2$queryHits]
hits2$ensembl_gene_id = TSS$ensembl_gene_id[hits2$subjectHits]
hits2$pair = paste0(hits2$variant_id, "_", hits2$ensembl_gene_id)

Union_CortexAndBA9_hg19$pair = paste0(Union_CortexAndBA9_hg19$variant_id, "_", Union_CortexAndBA9_hg19$ensembl_gene_id)
idx_rm2 = which(Union_CortexAndBA9_hg19$pair %in% hits2$pair) # 395599 rows

# b. remove eQTL more than 500kb away from TSS
TSS_500kb = data_frame(chr = paste0("chr", geneAnno$chromosome_name), start = geneAnno$TSS - 500e3, end = geneAnno$TSS + 500e3, ensembl_gene_id = geneAnno$ensembl_gene_id)
TSS_500kb_gr = annoDF2GR(TSS_500kb)

hits1 = as.data.frame(findOverlaps(Variant_Location_hg19_gr, TSS_500kb_gr))
idx_rm1 = which(! rownames(Union_CortexAndBA9_hg19) %in% hits1$queryHits) # 2946 eQTL-gene pairs over 500kb away from TSS, 1Mb away reaches memory limit. xxx should be removed

idx_rm = unique(c(idx_rm1, idx_rm2)) # 397750 rows
Union_CortexAndBA9_hg19_filtered = Union_CortexAndBA9_hg19[-idx_rm2,] # 1425318 eQTL-gene pairs remain
Union_CortexAndBA9_hg19_filtered_gr = annoDF2GR(Union_CortexAndBA9_hg19_filtered[,-4])

# c. created 1000 random peak-to-gene link sets by randomly assigning these peaks to any gene within 500kb of the peak
sig_DistalATAC_location_uniq = unique(sig_DistalATAC_location) # 5793
sig_DistalATAC_location_uniq$ATACpeak = paste0(sig_DistalATAC_location_uniq$chr, "_", sig_DistalATAC_location_uniq$start, "_", sig_DistalATAC_location_uniq$end)
sig_DistalATAC_uniq_gr = annoDF2GR(sig_DistalATAC_location_uniq)

hits3 = as.data.frame(findOverlaps(sig_DistalATAC_uniq_gr, TSS_500kb_gr))
hits3$distal_ATAC = sig_DistalATAC_location_uniq$ATACpeak[hits3$queryHits]
hits3$ensembl_gene_id = TSS_500kb$ensembl_gene_id[hits3$subjectHits]
hits3 = unique(hits3[,3:4]) # 122703 random distal ATAC-gene pairs within 500kb
hits3$pairs = paste0(hits3$distal_ATAC, ":", hits3$ensembl_gene_id)

sig_CorATAC_uniq = unique(sig_CorATAC[,c("distal_ATAC", "ensembl_gene_id")]) # 8043 unique distal ATAC-GE pairs with significant correlation between distal and promoter ATAC intensity
sig_CorATAC_uniq$pairs = paste0(sig_CorATAC_uniq$distal_ATAC, ":", sig_CorATAC_uniq$ensembl_gene_id)
sig_CorATAC_uniq_within500kb = sig_CorATAC_uniq[which(sig_CorATAC_uniq$pairs %in% hits3$pairs), ] # 7872 unique distal ATAC-gene pairs within 500kb

Overlap_n = c()
set.seed(294)
maxB = 1000
for (B in 1:maxB) {
  random_idx = sample(1:nrow(hits3), nrow(sig_CorATAC_uniq_within500kb))
  random_CorATAC = hits3[random_idx,] # 7872 rows
  
  # how many of these overlap eQTL-GE pairs
  distalATAC_location = as.data.frame(str_split_fixed(random_CorATAC$distal_ATAC, "_", 3))
  colnames(distalATAC_location) = c("chr", "start", "end")
  distalATAC_location$start = as.integer(distalATAC_location$start); distalATAC_location$end = as.integer(distalATAC_location$end)
  distalATAC_location$distal_ATAC = random_CorATAC$distal_ATAC
  distalATAC_location = unique(distalATAC_location) # 3803 rows
  distalATAC_location_gr = annoDF2GR(distalATAC_location)
  
  hits_random = as.data.frame(findOverlaps(distalATAC_location_gr, Union_CortexAndBA9_hg19_filtered_gr))
  hits_random$distal_ATAC = distalATAC_location$distal_ATAC[hits_random$queryHits]
  hits_random$ensembl_gene_id = Union_CortexAndBA9_hg19_filtered$ensembl_gene_id[hits_random$subjectHits]
  hits_random = unique(hits_random[,3:4]) # 2618 rows
  
  Overlap_B = inner_join(random_CorATAC, hits_random) # 169
  Overlap_n = c(Overlap_n, nrow(Overlap_B))
}

#hist(Overlap_n, breaks = 20)

# d. calculated the z-score and enrichment of eQTL-gene pairs in our determined peak-to-gene links compared to the randomized peak sets. 
# how many of sig_CorATAC_uniq overlap eQTL-GE pairs
distalATAC_location = as.data.frame(str_split_fixed(sig_CorATAC_uniq$distal_ATAC, "_", 3))
colnames(distalATAC_location) = c("chr", "start", "end")
distalATAC_location$start = as.integer(distalATAC_location$start); distalATAC_location$end = as.integer(distalATAC_location$end)
distalATAC_location$distal_ATAC = sig_CorATAC_uniq$distal_ATAC
distalATAC_location = unique(distalATAC_location) # 5793 rows
distalATAC_location_gr = annoDF2GR(distalATAC_location)

hits_sigATAC = as.data.frame(findOverlaps(distalATAC_location_gr, Union_CortexAndBA9_hg19_filtered_gr))
hits_sigATAC$distal_ATAC = distalATAC_location$distal_ATAC[hits_sigATAC$queryHits]
hits_sigATAC$ensembl_gene_id = Union_CortexAndBA9_hg19_filtered$ensembl_gene_id[hits_sigATAC$subjectHits]
hits_sigATAC = unique(hits_sigATAC[,3:4]) # 3213 rows

Overlap_sigATAC = inner_join(sig_CorATAC_uniq, hits_sigATAC)
nrow(Overlap_sigATAC) # 297 -> 800 -> 1313 ->
# Significant overlap between significant distal-promoter ATAC pairs and eQTL-GE pairs :-)

(eQTL_distalATAC_p = length(which(Overlap_n > nrow(Overlap_sigATAC)))/maxB)
# ATAC cor < 0.05: p < 1/1000

pdf("10_3_8_v5_corrected20230308_SignificantOverlapBw_eQTL_DistalATAC_ATACcorFDR02.pdf", height = 4, width = 6)

hist(Overlap_n, breaks = 20, main = paste0("p < 1e-3"), xlim = c(min(Overlap_n), max(nrow(Overlap_sigATAC), max(Overlap_n))), xlab = "Number of distal ATAC-gene pairs overlap eQTL-gene pairs") # FDR < 0.05 to 0.2
abline(v = nrow(Overlap_sigATAC), col = "red")

dev.off()
## Great, in all FDR conditions, the distal ATAC-GE pairs have significant overlap with eQTL-GE pairs. 

## How about ATAC-GE logFC?
# a.Distal ATAC - linked with promoter ATAC
lnames = load("../../5_DiffATAC/DiffATAC.rda") # "DiffATAC" 

library(readxl)
Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

sig_CorATAC$promoterATAC_logFC = DiffATAC$ATAC_logFC[match(sig_CorATAC$promoter_ATAC, rownames(DiffATAC))]
sig_CorATAC$promoterATAC_FDR = DiffATAC$ATAC_FDR[match(sig_CorATAC$promoter_ATAC, rownames(DiffATAC))]
sig_CorATAC$distalATAC_logFC = DiffATAC$ATAC_logFC[match(sig_CorATAC$distal_ATAC, rownames(DiffATAC))]
sig_CorATAC$distalATAC_FDR = DiffATAC$ATAC_FDR[match(sig_CorATAC$distal_ATAC, rownames(DiffATAC))]
sig_CorATAC$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
sig_CorATAC$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

library(ggplot2)

sig_CorATAC_uniq = unique(sig_CorATAC[,c("distal_ATAC", "ensembl_gene_id", "distalATAC_logFC", "distalATAC_FDR", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")]) # 8043
idx_plot = which(sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05)
lm_res = summary(lm(WholeCortex_ASD_logFC ~ 0 + distalATAC_logFC, data = sig_CorATAC_uniq[idx_plot, ]))
lm_coef = lm_res$coefficients[1]
lm_p = lm_res$coefficients[4]
lm_r2 = lm_res$r.squared

#pdf("ATAC_GE_logFC_SEPdistalORpromoterATAC_FDRcorATAC005.pdf", height = 4, width = 6)
#pdf("ATAC_GE_logFC_SEPdistalORpromoterATAC_FDRcorATAC01.pdf", height = 4, width = 6)
#pdf("ATAC_GE_logFC_SEPdistalORpromoterATAC_FDRcorATAC015.pdf", height = 4, width = 6)
#pdf("ATAC_GE_logFC_SEPdistalORpromoterATAC_FDRcorATAC02.pdf", height = 4, width = 6)
sig_CorATAC_uniq[idx_plot, ] %>%
  ggplot(aes(x = distalATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point() + # col = "grey"
  #geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = lm_coef, col = "red") +
  theme_bw() + 
  ggtitle("Distal ATAC significantly correlated with promoter ATAC") +
  annotate("text", x = 0, y = range(sig_CorATAC_uniq$WholeCortex_ASD_logFC[idx_plot])[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2),2), "\n", ifelse(lm_p < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p, format = "e", digits = 1)))), col = "red")
# much better than eQTL ATAC

sig_CorATAC[sig_CorATAC$distalATAC_FDR < 0.05 & sig_CorATAC$WholeCortex_ASD_FDR < 0.05, ] %>%
  ggplot(aes(x = promoterATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point() + # col = "grey"
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() + 
  ggtitle("Promoter ATAC with significantly correlated differential distal ATAC")
# much better than eQTL ATAC

sig_CorATAC[sig_CorATAC$promoterATAC_FDR < 0.05 & sig_CorATAC$WholeCortex_ASD_FDR < 0.05, ] %>%
  ggplot(aes(x = promoterATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point() + # col = "grey"
  geom_smooth(method='lm', formula= y~x) +
  theme_bw() + 
  ggtitle("Promoter ATAC with significantly correlated distal ATAC")
# nope, negative cor

dev.off()
# corATAC fdr < 0.2: cor = 0.12 and p-val=6.9e-4
# corATAC fdr < 0.15: cor = 0.11 and p-val=8.6e-3
# corATAC fdr < 0.1: cor = 0.24 and p-val=1.3e-5 * looks the best
# corATAC fdr < 0.05: cor = 0.28 and p-val=2.3e-03

save(sig_CorATAC, sig_CorATAC_uniq, file = "Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC005.rda")
save(sig_CorATAC, sig_CorATAC_uniq, file = "Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda")
save(sig_CorATAC, sig_CorATAC_uniq, file = "Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC015.rda")
save(sig_CorATAC, sig_CorATAC_uniq, file = "Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC02.rda")
load("Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda")
nrow(sig_CorATAC_uniq) # 21067

## How do distal ATAC-GE logFC look like, if I do not filter for the significantly correlated ATAC-ATAC pairs?
lnames = load("Cor_Gene_ATAC_calcedCor.rda") # Cor_Gene_ATAC

Cor_Gene_ATAC$promoterATAC_logFC = DiffATAC$ATAC_logFC[match(Cor_Gene_ATAC$promoter_ATAC, rownames(DiffATAC))]
Cor_Gene_ATAC$promoterATAC_FDR = DiffATAC$ATAC_FDR[match(Cor_Gene_ATAC$promoter_ATAC, rownames(DiffATAC))]
Cor_Gene_ATAC$distalATAC_logFC = DiffATAC$ATAC_logFC[match(Cor_Gene_ATAC$distal_ATAC, rownames(DiffATAC))]
Cor_Gene_ATAC$distalATAC_FDR = DiffATAC$ATAC_FDR[match(Cor_Gene_ATAC$distal_ATAC, rownames(DiffATAC))]
Cor_Gene_ATAC$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(Cor_Gene_ATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
Cor_Gene_ATAC$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(Cor_Gene_ATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

Cor_Gene_ATAC_uniq = unique(Cor_Gene_ATAC[,c("distal_ATAC", "ensembl_gene_id", "distalATAC_logFC", "distalATAC_FDR", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")]) # 195359
idx_plot = which(Cor_Gene_ATAC_uniq$distalATAC_FDR < 0.05 & Cor_Gene_ATAC_uniq$WholeCortex_ASD_FDR < 0.05)

lm_res = summary(lm(WholeCortex_ASD_logFC ~ 0 + distalATAC_logFC, data = Cor_Gene_ATAC_uniq[idx_plot, ]))
lm_coef = lm_res$coefficients[1]
lm_p = lm_res$coefficients[4]
lm_r2 = lm_res$r.squared

pdf("ATAC_GE_logFC_SEPdistalORpromoterATAC_FDRcorATAC1.pdf", height = 4, width = 6)
Cor_Gene_ATAC_uniq[idx_plot,] %>%
  ggplot(aes(x = distalATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point() + # col = "grey"
  #geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = lm_coef, col = "red") +
  theme_bw() + 
  ggtitle("Distal ATAC linked with promoter ATAC\nby Hi-C loop or within 30kb of gene TSS") +
  annotate("text", x = 0, y = range(Cor_Gene_ATAC_uniq$WholeCortex_ASD_logFC[idx_plot])[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2),2), "\n", ifelse(lm_p < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p, format = "e", digits = 1)))), col = "red")
dev.off()
# corATAC fdr <= 1: cor = 0.05 and p-val=1.9e-2. The correlation indeed improved by filtering for correlated ATAC-ATAC pairs (FDR < 0.05 to 0.2). Yay!

## b. All promoter ATAC-GE pairs
# Load all promoter ATAC peaks
lnames = load("promATAC_TSSminus100bptoplus2kb.rda") # promATACs
colnames(promATACs)[1] = "ATACpeak"
promATAC = unique(promATACs[,1:2]) # 28250 rows
length(unique(promATAC$ensembl_gene_id)) # 18610

# promoter ATAC-GE logFC
promATAC$ATAC_logFC = DiffATAC$ATAC_logFC[match(promATAC$ATACpeak, rownames(DiffATAC))]
promATAC$ATAC_FDR = DiffATAC$ATAC_FDR[match(promATAC$ATACpeak, rownames(DiffATAC))]
promATAC$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(promATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
promATAC$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(promATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

lm_res_prom = summary(lm(WholeCortex_ASD_logFC ~ 0 + ATAC_logFC, data = promATAC[promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05, ]))
lm_coef_prom = lm_res_prom$coefficients[1]
lm_p_prom = lm_res_prom$coefficients[4]
lm_r2_prom = lm_res_prom$r.squared

pdf("ATAC_GE_logFC_AllPromoterATAC.pdf", height = 4, width = 6)
promATAC[promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05, ] %>%
  ggplot(aes(x = ATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point() + # col = "grey"
  #geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = lm_coef_prom, col = "red") +
  theme_bw() + 
  ggtitle("Promoter ATAC (gene TSS -2kb to +100bp)") +
  annotate("text", x = 0, y = range(promATAC$WholeCortex_ASD_logFC[promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2_prom),2), "\n", ifelse(lm_p_prom < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p_prom, format = "e", digits = 1)))), col = "red")
dev.off()
# Many in Q4, although positive and significant cor in general
# Promoter ATAC as TSS+-2kb get cor=0.13, p-val = 8.3e-03
# Promoter ATAC as TSS-2kb to TSS get cor = 0.17 and p-val=1.4e-3
# Promoter ATAC as TSS-2kb to TSS+100bp get cor = 0.17 and p-val=1.2e-3

# Label promoter ATAC that also serve as distal ATAC
promATAC$AsDistalATAC = ifelse(promATAC$ATACpeak %in% sig_CorATAC_uniq$distal_ATAC, "yes", "no")
idx_plot = which(promATAC$ATAC_FDR < 0.05 & promATAC$WholeCortex_ASD_FDR < 0.05)

pdf("ATAC_GE_logFC_AllPromoterATAC_colPromATACasDistalATAC.pdf", height = 4, width = 6)
promATAC[idx_plot, ] %>%
  ggplot(aes(x = ATAC_logFC, y = WholeCortex_ASD_logFC)) +
  geom_point(aes(col = AsDistalATAC)) + # col = "grey"
  #geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = lm_coef_prom, col = "black") +
  theme_bw() + 
  ggtitle("Promoter ATAC (gene TSS -2kb to +100bp)") +
  annotate("text", x = 0, y = range(promATAC$WholeCortex_ASD_logFC[idx_plot], na.rm = T)[2] * 0.8, label = paste0("cor = ", round(sqrt(lm_r2_prom),2), "\n", ifelse(lm_p_prom < 2e-16, "p < 2e-16", paste0("p = ", formatC(lm_p_prom, format = "e", digits = 1)))), col = "black")
dev.off()

## Save all results
load("Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda")
save(sig_CorATAC, sig_CorATAC_uniq, promATAC, file = "Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda")
# before -> after correction
nrow(sig_CorATAC_uniq) # 23594 -> 21067
nrow(promATAC) # 28250 -> no change
nrow(DiffATAC) # 128036
length(unique(c(sig_CorATAC_uniq$distal_ATAC, promATAC$ATACpeak))) # 32906 -> 31797 unique ATAC peaks, 25% of all ATAC peaks

# ---- for spread sheet ----
lnames = load("Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC005.rda") # sig_CorATAC, sig_CorATAC_uniq
lnames = load("Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda") # sig_CorATAC, sig_CorATAC_uniq, promATAC
nrow(sig_CorATAC_uniq)
nrow(sig_CorATAC_uniq[which(sig_CorATAC_uniq$distalATAC_FDR < 0.05 & sig_CorATAC_uniq$WholeCortex_ASD_FDR < 0.05),])
# --------------------------

## Decision: 
# 1) Use the results from this script: Promoter ATAC as TSS -2kb to +100bp, Distal ATAC significantly correlated with promoter ATAC (FDR < 0.1), either Hi-C linked or within 30kb of TSS. 
# 2) For Number of distalATAC-GE pairs, cor and p of diffATAC-DEG logFC, see excel sheet 10_3_8_ATACcor_summary.xlsx

## Conclusion: 
# 1) Both differential promoter and distal ATAC peaks show significant positive correlation with differential gene expression.
# 2) Distal ATAC changes may correlate with GE changes even better than promoter ATAC changes, as there are few dots in Q4 for distal ATAC.

# Next step: 
# 1) Subset for CTCF-bound or BATF/FOS TF cluster-bound ATAC, watch the logFC cor between diffATAC and DEG.


