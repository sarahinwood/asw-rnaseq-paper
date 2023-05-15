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
attack_location_tpms <- snakemake@input[["attack_location_tpms"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

dds$pc1_sign <- factor(paste(dds$PC1_sign_together))
dds$Location <- factor(paste(dds$Location))
dds$Parasitism <- factor(paste(dds$Parasitism))
dds$Parasitism <- relevel(dds$Parasitism, ref="Undetected")
# add attack status
dds$Attacked <- factor(ifelse(dds$Attack=="ANO", "N", "Y")) # ANO = attack not observed
dds$Attacked <- relevel(dds$Attacked, ref="Y")

dds_attacked_int <- copy(dds)
design(dds_attacked_int) <- ~pc1_sign+Parasitism+Location+Attacked+Location:Attacked
dds_attacked_int <- DESeq(dds_attacked_int)
saveRDS(dds_attacked_int, snakemake@output[["dds_file"]])

resultsNames(dds_attacked_int)
res_group <- results(dds_attacked_int, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
resLFC <- lfcShrink(dds_attacked_int, coef=6, res=res_group)

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
# read in annots and merge
annotations <- fread(annotation_file, na.strings="")
sig_annots <- merge(ordered_sig_res_group_table, annotations, by.x="rn", by.y="#gene_id", all.x=TRUE)
# read TPMs and merge
attack_location_tpms <- fread(attack_location_tpms)
sig_annots_tpm <- merge(sig_annots, attack_location_tpms, by.x="rn", by.y="#gene_id")
fwrite(sig_annots_tpm, snakemake@output[["sig_annots"]])

#############
## HEATMAP ##
#############

##vst transform
#vst <- varianceStabilizingTransformation(dds_attacked_int, blind=TRUE)
#vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)
##subset for DEGs
#vst_degs <- subset(vst_assay_dt, rn %in% sig_annots$rn)
##turn first row back to row name
#vst_degs <- vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
#vst_degs_plot <- vst_degs[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,27,28,32,37,38,41, # Dun attacked
#                             23,25,26,29,30,31,33,34,35,36,39,40, # Dun ANO
#                             42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61, #Ru attacked
#                             62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79)] #Ru ANO
##get location label info
#sample_to_label <- data.table(data.frame(colData(dds_attacked_int)[,c("Location", "Attacked", "sample_name")]))
#sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")

#location_colours <- list(Location = c(Dunedin="#fcfdbf", Ruakura="#fc8961"),
#                         Attacked=c(N = "#b73779", Y="#51127c"))

##plot
##not clustered by sample
#pdf(snakemake@output[["row_clustered_heatmap"]])
#pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
#         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
#         show_colnames = FALSE, border_color=NA, color=viridis(50))
#dev.off()

#pdf(snakemake@output[["col_clustered_heatmap"]])
#pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
#         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
#         show_colnames = FALSE, border_color=NA, color=viridis(50))
#dev.off()

# write log
sessionInfo()