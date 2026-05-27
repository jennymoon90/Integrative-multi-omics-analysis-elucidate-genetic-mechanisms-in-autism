rm(list = ls())
options(stringsAsFactors = F)
setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_3_11_OverlapScRNAregulon/")

### Load Lucy's scRNA regulon
library(readxl)
Lucy = read.csv("~/Documents/Documents/Geschwind_lab/LAB/Database_download/36_LucyBicks_ASDscRNA_regulons/Whole_regulons.csv")
Lucy = Lucy[,-1]

# Lucy$AR # lots of blanks. 
Lucy$AR[Lucy$AR != ""] # Looks good
# To extract a regulon, do Lucy$AR[Lucy$AR != ""]

all_Lucy = unique(c(as.matrix(Lucy)))
all_Lucy = all_Lucy[all_Lucy != ""] # 12348 genes in Lucy's table

Jill_DEG = read_excel("/Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx")
length(intersect(all_Lucy, Jill_DEG$external_gene_name)) # 11791 genes

### Load TF bound ATAC target genes
lnames = load("../10_3_10_OverlapATACwTFbound/save_all.rda")

### Reasoning on how to assess overlap ###
# Lucy's regulon is calculated using SCENIC (Aibar 2017 Nat Methods), which first identifies gene co-expression modules using GENIE3, next identifies TFBS over-represented in each gene module and subset their predicted candidate target genes using RcisTarget, finally scores the activity of each regulons in each cell with AUCell.
# So, DEG doesn't matter. I shall overlap Lucy's regulon (containing co-expression info) with TF-bound ATAC target genes. And if there's a significant overlap, then I can see in which cell clusters the regulon is up/down-reg in ASD.

### 1) Overlap target genes with promoter ATAC bound by BATF/FOSL1/BACH1/NFE2 with Lucy's regulon
# Lucy's regulon does not have BATF, FOSL1, NFE2 as head of regulon. There is BACH1 regulon.
Genes_promboundby_BACH1 = unique(promATAC$external_gene_name[which(promATAC$BACH1 == "yes")]) # 4493 genes
Lucy_BACH1 = Lucy$BACH1[Lucy$BACH1 != ""] # 431 genes.
both = intersect(Genes_promboundby_BACH1, Lucy_BACH1) # 144 genes, 30% of Lucy's BACH1 regulon. But this is just saying 30%  of Lucy's BACH1 regulon have promoter binding of BACH1.

## Fisher's exact test for overlap significance
# [incorrect]
#Fisher_res = fisher.test(matrix(c(length(both), length(Lucy_BACH1) - length(both), length(Genes_promboundby_BACH1) - length(both), 0), nrow = 2))
#Fisher_p = formatC(Fisher_res$p.value, format = "e", digits = 1) # 0, there is significant overlap between 

# [correct]
All = intersect(all_Lucy, Jill_DEG$external_gene_name) # 11791 genes 
LucyOnly = intersect(setdiff(Lucy_BACH1, both), All) # 272 genes
ATAConly = intersect(setdiff(Genes_promboundby_BACH1, both), All) # 2871 genes
Neither = length(All) - length(both) - length(LucyOnly) - length(ATAConly) # 8504 genes
Fisher_res = fisher.test(matrix(c(length(both), length(LucyOnly), length(ATAConly), Neither), nrow = 2))
Fisher_p = formatC(Fisher_res$p.value, format = "e", digits = 1) # 2.8e-5. As expected, there is significant overlap between Lucy's BACH1 regulon and my BACH1-bound promoter ATAC targets

## UpsetR plot
library(UpSetR)
list_BACH1bound = list(Genes_promboundby_BACH1 = Genes_promboundby_BACH1, Lucy_BACH1_regulon = Lucy_BACH1)

pdf("UpSetR_Overlap_BACH1boundATACpromTargetGenes_LucyRegulon.pdf", width = 8, height = 5, onefile = F)
upset(fromList(list_BACH1bound), order.by = "freq")
dev.off()

## Venn diagram
library(VennDiagram)
library(gridExtra)

save(Genes_promboundby_BACH1, All, Lucy_BACH1, both, Fisher_res, Fisher_p, file = "VennDiagram_Overlap_BACH1boundATACpromTargetGenes_LucyRegulon.rda")

png("VennDiagram_Overlap_BACH1boundATACpromTargetGenes_LucyRegulon.png", height = 500, width = 500)
grid.newpage()
my_venndiag <- draw.pairwise.venn(length(intersect(Genes_promboundby_BACH1, All)), length(intersect(Lucy_BACH1, All)), length(both), fill = c("orange","chartreuse"), alpha = 0.6, category = c("BACH1-bound promoter\nATAC target genes", "BACH1\nscRNA-seq\nregulon"), cat.pos = c(0,5), cat.dist = c(0.05, 0.06), cex = 2, cat.cex = 1.5, fontfamily = "Arial", cat.fontfamily = "Arial")
grid.arrange(gTree(children = my_venndiag), # Add title & subtitle
             #top = "My Main Title",
             bottom = textGrob(paste0("Significant overlap, p = ", Fisher_p), gp = gpar(fontsize=20)))
dev.off()
dev.off()

# Conclusion for BACH1 targets:See significant overlap between promoter BACH1 target genes and Lucy's BACH1 regulon, indicating BACH1 binding at the promoter regulates target gene co-expression.

### 2) Overlap target genes with promoter and distal ATAC bound by CTCF with Lucy's CTCF regulon
Genes_promboundby_CTCF = unique(promATAC$external_gene_name[which(promATAC$CTCF_MA0139.1 == "yes" | promATAC$CTCF_MA1929.1 == "yes" | promATAC$CTCF_MA1930.1 == "yes")]) # 17123 genes
Genes_distalboundby_CTCF = unique(sig_CorATAC_uniq$external_gene_name[which(sig_CorATAC_uniq$CTCF_MA0139.1 == "yes" | sig_CorATAC_uniq$CTCF_MA1929.1 == "yes" | sig_CorATAC_uniq$CTCF_MA1930.1 == "yes")]) # 8035 genes
Lucy_CTCF = Lucy$CTCF[Lucy$CTCF != ""] # 234 genes.
both_prom = intersect(Genes_promboundby_CTCF, Lucy_CTCF) # 208 genes, 89% of Lucy's CTCF regulon. 
both_distal = intersect(Genes_distalboundby_CTCF, Lucy_CTCF) # 130 genes, 56% of Lucy's CTCF regulon. 
both = intersect(union(Genes_promboundby_CTCF, Genes_distalboundby_CTCF), Lucy_CTCF) # 210 genes, 90% of Lucy's CTCF regulon. Since Lucy's regulon uses cis-TFBS, distal sites won't be in the regulon. 

## Fisher's exact test for overlap significance
# Test the union of CTCF promoter and distal target genes
All = intersect(all_Lucy, Jill_DEG$external_gene_name) # 11791 genes 
LucyOnly = intersect(setdiff(Lucy_CTCF, both), All) # 13 genes
ATAConly = intersect(setdiff(union(Genes_promboundby_CTCF, Genes_distalboundby_CTCF), both), All) # 10717 genes
Neither = length(All) - length(both) - length(LucyOnly) - length(ATAConly) # 851 genes
Fisher_res = fisher.test(matrix(c(length(both), length(LucyOnly), length(ATAConly), Neither), nrow = 2))
Fisher_p = formatC(Fisher_res$p.value, format = "e", digits = 1) # 0.44. There is no significant overlap between Lucy's CTCF regulon and my CTCF-bound promoter and distal ATAC targets. The reason is that there are too many CTCF-bound ATAC target genes.

# Test promoter target genes
All = intersect(all_Lucy, Jill_DEG$external_gene_name) # 11791 genes 
LucyOnly = intersect(setdiff(Lucy_CTCF, both_prom), All) # 15 genes
ATAConly = intersect(setdiff(Genes_promboundby_CTCF, both_prom), All) # 10569 genes
Neither = length(All) - length(both_prom) - length(LucyOnly) - length(ATAConly) # 999 genes
Fisher_res_prom = fisher.test(matrix(c(length(both_prom), length(LucyOnly), length(ATAConly), Neither), nrow = 2))
Fisher_p_prom = formatC(Fisher_res_prom$p.value, format = "e", digits = 1) # 0.4. There is no significant overlap between Lucy's CTCF regulon and my CTCF-bound promoter ATAC targets. The reason is that there are too many CTCF-bound promoter ATAC target genes that do not belong to Lucy's list.

# Test distal target genes
All = intersect(all_Lucy, Jill_DEG$external_gene_name) # 11791 genes 
LucyOnly = intersect(setdiff(Lucy_CTCF, both_distal), All) # 93 genes
ATAConly = intersect(setdiff(Genes_distalboundby_CTCF, both_distal), All) # 5075 genes
Neither = length(All) - length(both_distal) - length(LucyOnly) - length(ATAConly) # 6493 genes
Fisher_res = fisher.test(matrix(c(length(both_distal), length(LucyOnly), length(ATAConly), Neither), nrow = 2))
Fisher_p = formatC(Fisher_res$p.value, format = "e", digits = 1) # 2.2e-5. Wow! There is significant overlap between Lucy's CTCF regulon and my CTCF-bound DISTAL ATAC targets.

# --- Calculate intersected genes stats (Upon Dan's suggestion) ----
# Test the intersection of CTCF promoter and distal target genes
All = intersect(all_Lucy, Jill_DEG$external_gene_name) # 11791 genes 
both_intersect = intersect(intersect(Genes_promboundby_CTCF, Genes_distalboundby_CTCF), Lucy_CTCF) # 128 genes
LucyOnly = intersect(setdiff(Lucy_CTCF, both_intersect), All) # 95 genes
ATAConly = intersect(setdiff(intersect(Genes_promboundby_CTCF, Genes_distalboundby_CTCF), both_intersect), All) # 4927 genes
Neither = length(All) - length(both_intersect) - length(LucyOnly) - length(ATAConly) # 6641 genes
Fisher_res_intersect = fisher.test(matrix(c(length(both_intersect), length(LucyOnly), length(ATAConly), Neither), nrow = 2))
Fisher_res_intersect$estimate # 1.8
Fisher_p_intersect = formatC(Fisher_res_intersect$p.value, format = "e", digits = 1) # 1e-5. There is significant overlap between Lucy's CTCF regulon and the intersection of my CTCF-bound promoter and distal ATAC targets. 
# ------------------------------------------------------------------

## UpsetR diagram
list_CTCFbound = list(Genes_distalboundby_CTCF = Genes_distalboundby_CTCF, Genes_promboundby_CTCF = Genes_promboundby_CTCF, Lucy_CTCF_regulon = Lucy_CTCF)

pdf("UpSetR_Overlap_CTCFboundATACpromANDdistalTargetGenes_LucyRegulon.pdf", width = 8, height = 5, onefile = F)
upset(fromList(list_CTCFbound), order.by = "freq")
dev.off()
# Most of distal CTCF-bound ATAC target genes are included in the promoter CTCF-bound ATAC target genes.

## Venn diagram
save(Genes_distalboundby_CTCF, All, Genes_promboundby_CTCF, Lucy_CTCF, Fisher_res, Fisher_p, Fisher_res_prom, Fisher_p_prom, Fisher_res_intersect, Fisher_p_intersect, file = "VennDiagram_Overlap_CTCFboundATACpromANDdistalTargetGenes_LucyRegulon.rda")

png("VennDiagram_Overlap_CTCFboundATACpromANDdistalTargetGenes_LucyRegulon.png", height = 500, width = 500)
#grid.newpage()
my_venndiag <- draw.triple.venn(length(intersect(Genes_distalboundby_CTCF, All)), length(intersect(Genes_promboundby_CTCF, All)), length(intersect(Lucy_CTCF, All)), n12 = length(intersect(intersect(Genes_distalboundby_CTCF, Genes_promboundby_CTCF), All)), n13 = length(intersect(intersect(Genes_distalboundby_CTCF, Lucy_CTCF), All)), n23 = length(intersect(intersect(Genes_promboundby_CTCF, Lucy_CTCF), All)), n123 = length(intersect(intersect(Genes_distalboundby_CTCF, Genes_promboundby_CTCF), Lucy_CTCF)), fill = c("orange", "dodgerblue", "chartreuse"), alpha = 0.6, category = c("CTCF-bound distal\nATAC target genes", "CTCF-bound promoter\nATAC target genes", "CTCF scRNA-seq regulon"), cat.pos = c(-10,10,0), cat.dist = c(0.1, 0.1, -0.45), cex = 2, cat.cex = 1.5, fontfamily = "Arial", cat.fontfamily = "Arial")
grid.arrange(gTree(children = my_venndiag), # Add title & subtitle
             #top = "My Main Title",
             bottom = textGrob(paste0("Significant overlap only between CTCF distal\ntarget genes and scRNA-seq regulon, p = ", Fisher_p), gp = gpar(fontsize=20)), heights = c(30,1)) # bottom = textGrob(..., (vjust = 0.7)) has title cut-off. Use heights = c(n,1) to adjust the venn diagram size and allow more space between the bottom title with the diagram
dev.off()

# Conclusion for CTCF targets: Only see significant overlap between distal CTCF targets and Lucy's CTCF regulon, indicating CTCF binding at distal enhancers regulates target gene co-expression.


