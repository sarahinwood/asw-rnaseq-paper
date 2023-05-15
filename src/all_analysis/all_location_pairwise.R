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
location_tpms <- snakemake@input[["location_tpms"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

dds$pc1_sign <- factor(paste(dds$PC1_sign_together))
dds$Location <- factor(paste(dds$Location))
dds$Experiment <- factor(paste(dds$Experiment))
dds$Parasitism <- factor(paste(dds$Parasitism))

dds_location <- copy(dds)
design(dds_location) <- ~Experiment+Parasitism+pc1_sign+Location
dds_location <- DESeq(dds_location)
saveRDS(dds_location, snakemake@output[["dds_file"]])

res_group <- results(dds_location, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
resLFC <- lfcShrink(dds_location, coef=5, res=res_group)

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
## merge with shrunk LFCs
shrunk_LFCs <- data.table(data.frame(resLFC), keep.rownames = T)
gene_to_shrunk_LFC <- shrunk_LFCs[,c(1,3,4)]
setnames(gene_to_shrunk_LFC, old=c("log2FoldChange", "lfcSE"), new=c("shrunk_log2FoldChange", "shrunk_lfcSE"))
full_res_group_table <- merge(ordered_res_group_table, gene_to_shrunk_LFC)
ordered_full_res_group_table <- full_res_group_table[order(full_res_group_table$padj),]
fwrite(ordered_full_res_group_table, snakemake@output[["res_table"]])

##sig with annots
ordered_sig_res_group_table <- subset(ordered_full_res_group_table, padj < 0.05)
ordered_sig_res_group_table$abs_shrunk_LFC_above_1 <- ifelse(abs(ordered_sig_res_group_table$shrunk_log2FoldChange>1), TRUE, FALSE)
# read annotations & merge
annotations <- fread(annotation_file, na.strings="")
sig_annots <- merge(ordered_sig_res_group_table, annotations, by.x="rn", by.y="#gene_id", all.x=TRUE)
# read TPMs and merge
location_tpms <- fread(location_tpms)
sig_annots_tpm <- merge(sig_annots, location_tpms, by.x="rn", by.y="#gene_id")
fwrite(sig_annots_tpm, snakemake@output[["sig_annots"]])

#############
## HEATMAP ##
#############

##vst transform
vst <- varianceStabilizingTransformation(dds_location, blind=TRUE)
vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)
##subset for DEGs
vst_degs <- subset(vst_assay_dt, rn %in% sig_annots$rn)
##turn first row back to row name
vst_degs <- vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
vst_degs_plot <- vst_degs[,c(1:53, 66:103, 54:65)]

##get location label info
sample_to_label <- data.table(data.frame(colData(dds_location)[,c("Location", "Experiment", "sample_name")]))
sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Location = c(Dunedin="#fcfdbf", Ruakura="#fc8961"),
                         Experiment=c(Evasion = "#b73779", Exposure="#51127c"))

##plot
##not clustered by sample
pdf(snakemake@output[["row_clustered_heatmap"]])
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

pdf(snakemake@output[["col_clustered_heatmap"]])
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()
# write log
sessionInfo()