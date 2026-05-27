# Modified from 20_1_HighConfidenceLoops_EachCellType_VaryingPercentSamples.R

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)

### How does UpSetR look like if I remove loops that anchor at segmental duplicate regions and remove loops <= 10kb distance?
# setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/20_HighConfidenceLoops_varyingPercentSamples_diffPromoterLoops/")
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/22_autosome_AllFithic2SettingsSame/q001_AllFithic2SettingsSame")
library(UpSetR)

# Consensus loops defined as those present in 8 Bulk, 7 NeuNn, 9 NeuNp samples. 
load("Bulk_consensus_Autosome_loops_inAtLeast8outof36samples.rda") # consensus_Bulk
load("NeuNp_consensus_Autosome_loops_inAtLeast9outof15samples.rda") # consensus_NeuNp
load("NeuNn_consensus_Autosome_loops_inAtLeast7outof15samples.rda") # consensus_NeuNn

range(consensus_Bulk$distance) # 10 kb - 242.89 Mb, the max distance is unreasonable
range(consensus_NeuNp$distance) # 10 kb - 242.87 Mb, the max distance is unreasonable
range(consensus_NeuNn$distance) # 10 kb - 223.48 Mb, the max distance is unreasonable

Bulk_loops = consensus_Bulk[,1:5]; Bulk_loops$Bulk = 1
NeuNp_loops = consensus_NeuNp[,1:5]; NeuNp_loops$NeuNp = 1
NeuNn_loops = consensus_NeuNn[,1:5]; NeuNn_loops$NeuNn = 1

Upset_tbl = full_join(Bulk_loops, NeuNp_loops)
Upset_tbl = full_join(Upset_tbl, NeuNn_loops)
Upset_tbl[is.na(Upset_tbl)] = 0
Upset_tbl = as.data.frame(Upset_tbl)
rownames(Upset_tbl) = Upset_tbl$loop

## Filter out loops that anchor at segmental duplicates
seg_dup = read.table("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_Dup15q/Hi-C/data_analysis/0_SegDups/SeqDup_build37_fracMatch_merged.bed")

seg_dup_bins = c()
for (i in 1:nrow(seg_dup)) {
  cur_segdup = seq(floor(seg_dup$V2[i]/10000) * 10000 + 5000, 
                   floor(seg_dup$V3[i]/10000) * 10000 + 5000,
                   10000)
  cur_segdupbin = paste0(seg_dup$V1[i], "_", cur_segdup)
  seg_dup_bins = c(seg_dup_bins, cur_segdupbin)
}
length(seg_dup_bins) # 22881 bins

Upset_tbl$start_anchor = paste0(Upset_tbl$V1, "_", Upset_tbl$V2)
Upset_tbl$end_anchor = paste0(Upset_tbl$V1, "_", Upset_tbl$V3)
idx_segdup = which(Upset_tbl$start_anchor %in% seg_dup_bins | Upset_tbl$end_anchor %in% seg_dup_bins) # 18207 loops anchor at segmental duplicates

hist(Upset_tbl$distance, breaks = 100, xlim = c(0,1e8))
hist(Upset_tbl$distance[idx_segdup], breaks = 100, xlim = c(0,1e8)) # relatively enriched at 25Mb and 60Mb distance

Upset_tbl = Upset_tbl[-idx_segdup,]
# After filtering, we have 185739 loops remain.

## Keep loops that are [20kb, 5Mb]
Upset_tbl = Upset_tbl[Upset_tbl$distance >= 20000,]
# After filtering, we have 16943 loops remain.
hist(Upset_tbl$distance, breaks = 1000, xlim = c(0,8e6))
Upset_tbl = Upset_tbl[Upset_tbl$distance <= 5e6,]

save(Upset_tbl, file = "UpsetTbl_8Bulk7NeuNn9NeuNp_LoopFilter_20kbto5MbRmSegdup.rda")

## Plot UpSetR
Upset_tbl2 = Upset_tbl[,6:8]

pdf("Upset_consensus8Bulk7NeuNn9NeuNp_AutosomeLoopsAtleast20kbAtmost5MbRmSegdup_bwCelltypes_AllSamplesRmTechRepAndOutliers.pdf", height = 6, width = 14) #
upset(Upset_tbl2, sets = colnames(Upset_tbl2),
      sets.bar.color = c("darkseagreen2", "pink", "#56B4E9"), 
      order.by = "freq", 
      text.scale = 2, show.numbers = "yes",
      # set_size.scale_max = max(c(nrow(Bulk_loops), nrow(NeuNp_loops), nrow(NeuNn_loops))) + 1000) #, 
      set_size.scale_max = 85e3) #,
#mainbar.y.label = "Number of shared TAD boundaries",
#sets.x.label = "Number of high-confidence\n TAD boundaries")
#sets.bar.color = c(rep(c("", ""), each = 2), rep("", 2))
dev.off()

# Interested in know whether differential loops in Bulk and NeuNp are coordinated.

# Dan suggests to draw a pie chart for each cell-type to show what proportion is unique, what proportion is shared with other cell-types.
# Run Figure3.R

### Loop cell-type specificity regarding gene expression and function
## Get loops in NeuNp(n) consensus set but never present in samples of the other celltype, ie. NeuNn(p).
rm(list = ls())
lnames = load("UpsetTbl_8Bulk7NeuNn9NeuNp_LoopFilter_20kbto5MbRmSegdup.rda") # Upset_tbl
lnames = load("Bulk_AutosomeLoop_by_sample_presence_matrix_AllSamplesRmTechRepAndOutliers.rda") # df_Bulk
lnames = load("NeuNp_AutosomeLoop_by_sample_presence_matrix_AllSamplesRmTechRepAndOutliers.rda") # df_NeuNp
lnames = load("NeuNn_AutosomeLoop_by_sample_presence_matrix_AllSamplesRmTechRepAndOutliers.rda") # df_NeuNn

df_NeuNp$n_samples = apply(df_NeuNp[,6:20], 1, sum)
df_NeuNn$n_samples = apply(df_NeuNn[,6:20], 1, sum)

NeuNp_consens = Upset_tbl$loop[Upset_tbl$NeuNp == 1] # 79504
NeuNn_consens = Upset_tbl$loop[Upset_tbl$NeuNn == 1] # 80876
NeuNp_specific = NeuNp_consens[which(NeuNp_consens %in% df_NeuNn$loop[df_NeuNn$n_samples == 0])] # 13671
NeuNn_specific = NeuNn_consens[which(NeuNn_consens %in% df_NeuNp$loop[df_NeuNp$n_samples == 0])] # 6265

## Intersect promoter
# Copied from 20_6_CelltypeSpecificPromoterLoops_enrichment.R
library(EnsDb.Hsapiens.v75)
library(Repitools)
#library(diffloop)

biotype = c("lincRNA", "miRNA", "protein_coding", "antisense", "snRNA", "sense_intronic", "snoRNA", "rRNA", "sense_overlapping", "3prime_overlapping_ncrna")
txs <- transcripts(EnsDb.Hsapiens.v75, columns = c("uniprot_id", "tx_biotype", "gene_id", "gene_name"), filter = c(TxBiotypeFilter(biotype), SeqNameFilter(1:22))) # do it without filter = AnnotationFilter(~ tx_biotype %in% biotype) first to select biotypes.
#txs = addchr(txs) 
length(txs) # 190463 non-pseudo genes in the genome
names(txs) = 1:length(txs)
txs_df = annoGR2DF(txs)
txs_df$chr = paste0("chr", txs_df$chr)
unique(txs_df$tx_biotype)

txs_df$tss = ifelse(txs_df$strand == "+", txs_df$start, txs_df$end)
res = 10000
txs_df$promoter_bin = floor(txs_df$tss/res) * res + res/2
txs_df$chr_promoterbin = paste0(txs_df$chr, "_", txs_df$promoter_bin)

save(txs_df, file = "../txs_df.rda")
load("../txs_df.rda")
promoter_bins = unique(txs_df$chr_promoterbin) # 43300 promoter bins

NeuNp_specific_location = as.data.frame(str_split_fixed(NeuNp_specific, "_", 3))
NeuNn_specific_location = as.data.frame(str_split_fixed(NeuNn_specific, "_", 3))
NeuNp_specific_location$anchor1 = paste0(NeuNp_specific_location$V1, "_", NeuNp_specific_location$V2)
NeuNp_specific_location$anchor2 = paste0(NeuNp_specific_location$V1, "_", NeuNp_specific_location$V3)
NeuNn_specific_location$anchor1 = paste0(NeuNn_specific_location$V1, "_", NeuNn_specific_location$V2)
NeuNn_specific_location$anchor2 = paste0(NeuNn_specific_location$V1, "_", NeuNn_specific_location$V3)
rownames(NeuNp_specific_location) = NeuNp_specific
rownames(NeuNn_specific_location) = NeuNn_specific

NeuNp_specific_promoter_loops = NeuNp_specific_location[which(NeuNp_specific_location$anchor1 %in% promoter_bins | NeuNp_specific_location$anchor2 %in% promoter_bins),] # 5980 NeuNp-specific promoter loops
NeuNn_specific_promoter_loops = NeuNn_specific_location[which(NeuNn_specific_location$anchor1 %in% promoter_bins | NeuNn_specific_location$anchor2 %in% promoter_bins),] # 3880 NeuNn-specific promoter loops

save(NeuNp_specific_promoter_loops, NeuNn_specific_promoter_loops, file = "../loop_specificity/Celltype_specific_consensus_promoter_loops.rda")
lnames = load("../loop_specificity/Celltype_specific_consensus_promoter_loops.rda")
  
## Get the genes at the cell-type specific promoter loop anchors
NeuNp_specificPL_anchors = c(NeuNp_specific_promoter_loops$anchor1, NeuNp_specific_promoter_loops$anchor2)
NeuNn_specificPL_anchors = c(NeuNn_specific_promoter_loops$anchor1, NeuNn_specific_promoter_loops$anchor2)
NeuNp_specificPL_genes = txs_df[which(txs_df$chr_promoterbin %in% NeuNp_specificPL_anchors), c("gene_name", "gene_id")] # 18559 genes at NeuNp-specific promoter loops
NeuNn_specificPL_genes = txs_df[which(txs_df$chr_promoterbin %in% NeuNn_specificPL_anchors), c("gene_name", "gene_id")] # 22961 genes at NeuNn-specific promoter loops

# How many genes overlap
ov_genes = intersect(NeuNp_specificPL_genes$gene_id, NeuNn_specificPL_genes$gene_id) # Only 433 genes (2%) overlap, no need to tease them out.

## Cell Specificity Enrichment Analysis

## 1) GO enrichment analysis for cell-type specific promoter loops
library(gprofiler2)
library(readxl)
library(gridExtra)

Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")

gost_NeuNp_genes = gost(query = unique(NeuNp_specificPL_genes$gene_id), organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_NeuNp_genes = gost_NeuNp_genes$result # nervous system development, synaptic signaling
gost_NeuNn_genes = gost(query = unique(NeuNn_specificPL_genes$gene_id), organism = "hsapiens", source = "GO:BP", custom_bg = Jill_DEG$ensembl_gene_id, correction_method = "fdr")
GO_NeuNn_genes = gost_NeuNn_genes$result # metabolic process

# GOterms_NeuNp = GO_NeuNp_genes; GOterms_NeuNn = GO_NeuNn_genes
GOterms_NeuNp = GO_NeuNp_genes[GO_NeuNp_genes$term_size < 10000,]
GOterms_NeuNn = GO_NeuNn_genes[GO_NeuNn_genes$term_size < 10000,]
#GOterms_NeuNp = GO_NeuNp_genes[(! GO_NeuNp_genes$term_name %in% GO_NeuNn_genes$term_name) & GO_NeuNp_genes$term_size < 10000,] # no need
#GOterms_NeuNn = GO_NeuNn_genes[(! GO_NeuNn_genes$term_name %in% GO_NeuNp_genes$term_name) & GO_NeuNn_genes$term_size < 10000,] # no need

single_plot = function(GO_dataframe, GO_title) {
  plot_nrow = min(nrow(GO_dataframe), 10)
  plot_up <- GO_dataframe[1:plot_nrow,] %>%
    ggplot(aes(x = term_name, y = logFDR)) + # x = term.name for gProfileR, x = term_name for gprofiler2
    geom_col(fill = "#00BFC4") +
    geom_hline(yintercept = -log10(0.05), color = "red", size = 1) +
    coord_flip() +
    ggtitle(GO_title) +
    xlab(NULL) + ylab("-log10FDR") +
    theme_bw() +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          title = element_text(size = 30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 1))
}

GOterms_NeuNp$logFDR = -log10(GOterms_NeuNp$p_value)
GOterms_NeuNp = GOterms_NeuNp[order(GOterms_NeuNp$p_value),]
GOterms_NeuNp$term_name = factor(GOterms_NeuNp$term_name, levels = rev(GOterms_NeuNp$term_name))
plot_NeuNp_GOterms = 
  single_plot(GOterms_NeuNp,"Genes at NeuNp-specific promoter loops")

GOterms_NeuNn$logFDR = -log10(GOterms_NeuNn$p_value)
GOterms_NeuNn = GOterms_NeuNn[order(GOterms_NeuNn$p_value),]
GOterms_NeuNn$term_name = factor(GOterms_NeuNn$term_name, levels = rev(GOterms_NeuNn$term_name))
plot_NeuNn_GOterms =
  single_plot(GOterms_NeuNn,"Genes at NeuNn-specific promoter loops")

# png("../loop_specificity/GOenrichment_genes_at_NeuNp_and_NeuNn_specific_promoterLoops_20kbto5MbRmSegdup.png", width = 2200, height = 800) # 
png("../loop_specificity/GOenrichment_genes_at_NeuNp_and_NeuNn_specific_promoterLoops_20kbto5MbRmSegdup_GOtermsizeMax10000.png", width = 2200, height = 800) # 
grid.arrange(plot_NeuNp_GOterms, plot_NeuNn_GOterms, nrow = 1)
dev.off()

# Observation:
# Do not limit term size: NeuNp-specific promoter loops are enriched for genes involved in nervous system development, signaling, etc. NeuNn-specific promoter loops are enriched for genes involved in metabolic process, etc. See some overlapping terms including: biological_process, cellular rocess.
# Limit term size to 10,000: now looks more specific.

## b) pSI enrichment analysis for cell-type specific promoter loops
lnames = load("~/Documents/Documents/Geschwind_lab/LAB/Database_download/Zhang_2015_pSI/Zhang2015pSI_usingMatureOfAllAges.rda")
colnames(pSI_res)

# In this case I should use Fisher's exact test (for association test), although Gokul uses logistic regression (for prediction).
Genes = data_frame(gene_id = c(NeuNp_specificPL_genes$gene_id, NeuNn_specificPL_genes$gene_id),
                   gene_name = c(NeuNp_specificPL_genes$gene_name, NeuNn_specificPL_genes$gene_name),
                   Loop = c(rep("NeuNp", nrow(NeuNp_specificPL_genes)),
                            rep("NeuNn", nrow(NeuNn_specificPL_genes))))

idx_Astro = which(Genes$gene_name %in% rownames(pSI_res)[which(pSI_res$Astrocyte<0.05)]) #
idx_Neuron = which(Genes$gene_name %in% rownames(pSI_res)[which(pSI_res$Neuron<0.05)]) #
idx_Oligo = which(Genes$gene_name %in% rownames(pSI_res)[which(pSI_res$Oligodendrocyte<0.05)]) #
idx_Microglia = which(Genes$gene_name %in% rownames(pSI_res)[which(pSI_res$Microglia<0.05)]) #
idx_Endoth = which(Genes$gene_name %in% rownames(pSI_res)[which(pSI_res$Endothelial<0.05)]) # 

Genes$SI_Astro = 
  Genes$SI_Neuron = 
  Genes$SI_Oligo = 
  Genes$SI_Microglia = 
  Genes$SI_Endoth = 
  Genes$SI_Glia = "non-specific"

# Genes$SI_Astro[idx_Astro] = "specific"
# Genes$SI_Neuron[idx_Neuron] = "specific"
# Genes$SI_Oligo[idx_Oligo] = "specific"
# Genes$SI_Microglia[idx_Microglia] = "specific"
# Genes$SI_Endoth[idx_Endoth] = "specific"
# Genes$SI_Glia[unique(c(idx_Astro, idx_Oligo, idx_Microglia))] = "specific"

Genes$SI_Endoth[idx_Endoth] = "specific"
Genes$SI_Microglia[idx_Microglia] = "specific"
Genes$SI_Oligo[idx_Oligo] = "specific"
Genes$SI_Astro[idx_Astro] = "specific"
Genes$SI_Neuron[idx_Neuron] = "specific"
Genes$SI_Glia[unique(c(idx_Astro, idx_Oligo, idx_Microglia))] = "specific"

p_Fisher = OR_Fisher = as.data.frame(matrix(nrow = 2, ncol = 6))
df = Genes
df$Quadrants = df$Loop
Quadrants = c("NeuNp", "NeuNn")
colnames(df)
k = which(colnames(df) == "SI_Glia") - 1

for (i in 1:2) { # Loop celltype
  for (j in 1:ncol(p_Fisher)) { # pSI Celltypes
    Q_Cell = length(which(df$Quadrants == Quadrants[i] & df[,k + j] == "specific"))
    Q_nonCell = length(which(df$Quadrants == Quadrants[i] & df[,k + j] != "specific"))
    nonQ_Cell = length(which(df$Quadrants != Quadrants[i] & df[,k + j] == "specific"))
    nonQ_nonCell = length(which(df$Quadrants != Quadrants[i] & df[,k + j] != "specific"))
    (Res = fisher.test(matrix(c(Q_Cell,Q_nonCell,nonQ_Cell,nonQ_nonCell), nrow = 2)))
    p_Fisher[i,j] =  format(Res$p.value, scientific = T, digits = 2)
    OR_Fisher[i,j] = format(Res$estimate, scientific = F, digits = 2)
  } 
}

OR_Fisher = as.data.frame(apply(OR_Fisher, 2, as.numeric))
p_Fisher = as.data.frame(apply(p_Fisher, 2, as.numeric))
colnames(p_Fisher) = colnames(OR_Fisher) = sapply(colnames(df)[(k+1):(k+6)], function(x) substr(x, 4, nchar(x)))
rownames(p_Fisher) = rownames(OR_Fisher) = Quadrants

OR_Fisher$Quadrant = rownames(OR_Fisher)
library(reshape2)
OR_Fisher_melt = reshape2::melt(OR_Fisher, value.name = "Odds_ratio", variable.name = "Celltype")
p_Fisher$Quadrant = rownames(p_Fisher)
p_Fisher_melt = reshape2::melt(p_Fisher, value.name = "P_value", variable.name = "Celltype")
Enrichment_df = left_join(OR_Fisher_melt, p_Fisher_melt)
Enrichment_df$FDR = p.adjust(Enrichment_df$P_value, method = "fdr") # All 1
Enrichment_df$FDR_chr = format(Enrichment_df$FDR, scientific = T, digits = 2)
Enrichment_df$label_text = paste0(Enrichment_df$Odds_ratio, "\n(", Enrichment_df$FDR_chr, ")")
Enrichment_df$label_text[which(Enrichment_df$Odds_ratio <= 1 | Enrichment_df$FDR > 0.1)] = NA
save(Enrichment_df, file = "../loop_specificity/NeuNp_and_NeuNn_specific_promoterLoops_20kbto5MbRmSegdup_Genes_pSIenrichment.rda")
# Only loops with Odds ratio > 1 & FDR <= 0.1 are shown

pdf("../loop_specificity/pSIenrichment_Genes_at_NeuNp_and_NeuNn_specific_promoterLoops_20kbto5MbRmSegdup.pdf", width = 6, height = 4)
Enrichment_df[Enrichment_df$Celltype != "Glia", ] %>% 
  ggplot(aes(x = Celltype, y = Quadrant)) +
  geom_tile(aes(fill = Odds_ratio)) + # , color = "black", size = 2
  scale_fill_gradient2(low = "white", high = "red", name = "Odds ratio", midpoint = 1, na.value = "red") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"), #panel.background = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank(),
        axis.text = element_text(size = 12)) +
  geom_text(aes(label=label_text)) +
  ggtitle("Enrichment of cell-type specifically expressed genes\nat NeuNp- or NeuNn-specific promoter loops")
dev.off()

## Observations:
# NeuNp-specific promoter loops are significantly enriched for neuron and endothelial -specifically expression genes 
# NeuNn-specific promoter loops are significantly enriched for oligodendrocyte, astrocyte,  and microglia cell-specifically expression genes 



