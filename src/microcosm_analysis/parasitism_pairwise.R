#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(DESeq2)
library(viridis)
library(tidyverse)
library(pheatmap)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]
annotation_file <- snakemake@input[["annotation_file"]]
parasitism_tpms_file <- snakemake@input[["parasitism_tpms_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

dds$pc1_sign <- factor(paste(dds$PC1_sign_together))
dds$Location <- factor(paste(dds$Location))
dds$Parasitism <- factor(paste(dds$Parasitism))
dds$Parasitism <- relevel(dds$Parasitism, ref="Undetected")

dds_parasitism <- copy(dds)
design(dds_parasitism) <- ~pc1_sign+Location+Parasitism
dds_parasitism <- DESeq(dds_parasitism)
saveRDS(dds_parasitism, snakemake@output[["dds_file"]])

resultsNames(dds_parasitism)
res_group <- results(dds_parasitism, lfcThreshold = 1, alpha = 0.05)
summary(res_group)

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, snakemake@output[["res_table"]])

##sig with annots
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
# read in annots and merge
annotations <- fread(annotation_file, na.strings="")
sig_annots <- merge(ordered_sig_res_group_table, annotations, by.x="rn", by.y="#gene_id", all.x=TRUE)
# read TPMs and merge
parasitism_tpms <- fread(parasitism_tpms_file)
sig_annots_tpm <- merge(sig_annots, parasitism_tpms, by.x="rn", by.y="#gene_id")
fwrite(sig_annots_tpm, snakemake@output[["sig_annots"]])

#############
## HEATMAP ##
#############

##vst transform
vst <- varianceStabilizingTransformation(dds_parasitism, blind=TRUE)
vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)
##subset for DEGs
vst_degs <- subset(vst_assay_dt, rn %in% sig_annots$rn)
##turn first row back to row name
vst_degs <- vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
vst_degs_plot <- vst_degs[,c(1,2,6,12,14,16,17,19,21,22,24,27,28,32,37,38,41, # Dun para
                             42,45,50,60, #Ru para
                             3:5,7:10,13,15,18,20,23,25,26,29:31,33:36,39,40, # Dun unpara
                             43,44,46:49,51:59, 61:79)] #ru unpara
##get location label info
sample_to_label <- data.table(data.frame(colData(dds_parasitism)[,c("Location", "Parasitism", "sample_name")]))
sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Location = c(Dunedin="#fcfdbf", Ruakura="#fc8961"),
                         Parasitism=c(Parasitised = "#b73779", Undetected="#51127c"))

##plot ordering
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

##plot
##not clustered by sample
pdf(snakemake@output[["row_clustered_heatmap"]])
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, clustering_callback=callback, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

pdf(snakemake@output[["col_clustered_heatmap"]])
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, clustering_callback=callback, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()



# write log
sessionInfo()