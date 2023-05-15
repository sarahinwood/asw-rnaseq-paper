library(data.table)
library(DESeq2)
library(viridis)
library(EnhancedVolcano)
library(VennDiagram)

##should be controlling for sequencing run

asw_dds <- readRDS("output/deseq2/asw_all_dds.rds")
##filter out genes with less than 10 reads mapping to them - helpful when transcriptome is large

asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign_together))
asw_dds$experiment <- factor(paste(asw_dds$experiment))
asw_dds$parasitism <- factor(paste(asw_dds$parasitism))
asw_dds$location <- factor(paste(asw_dds$location))
asw_dds_pc1 <- copy(asw_dds)
design(asw_dds_pc1) <- ~parasitism+experiment+location+pc1_sign
asw_dds_pc1 <- DESeq(asw_dds_pc1)

res_group <- results(asw_dds_pc1, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/pc1_all/pc1_res_group.csv")

##volcano plot w/viridis colours
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##sig with annots
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
trinotate <- fread("data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings="")
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
sig_annots$edited_rn <- tstrsplit(sig_annots$rn, "ASW_", keep=c(2))
fwrite(sig_annots, "output/deseq2/pc1_all/pc1_sig_annots.csv")

exposure_pc1 <- fread("../asw-exposed-rnaseq/output/deseq2/asw/pc1/sig_annots.csv")

vd1 <- venn.diagram(x = list("Exposure only"=exposure_pc1$rn,
                             "Both experiments"=sig_annots$edited_rn),
                    fill=c("#440154FF", "#FDE725FF"), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)
