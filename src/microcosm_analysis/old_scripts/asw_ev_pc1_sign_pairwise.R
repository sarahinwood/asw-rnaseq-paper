library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/evasion_deseq2/asw_evasion_dds.rds")
##filter out genes with less than 10 reads mapping to them - helpful when transcriptome is large
keep <- rowSums(counts(asw_dds)) >= 10
asw_dds_fil <- asw_dds[keep,]

asw_dds_fil$pc1_sign <- factor(paste(asw_dds_fil$PC1_sign))
asw_dds_pc1 <- copy(asw_dds_fil)

asw_dds_pc1 <- asw_dds_pc1[,asw_dds_pc1$parasitism == "undetected"]

design(asw_dds_pc1) <- ~pc1_sign

asw_dds_pc1 <- DESeq(asw_dds_pc1)

res_group <- results(asw_dds_pc1, contrast = c("pc1_sign", "negative", "positive"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/evasion_deseq2/res_group.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/evasion_deseq2/sig_degs.csv")

trinotate_report <- fread("data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/evasion_deseq2/sig_annots.csv")
