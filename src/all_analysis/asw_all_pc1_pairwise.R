library(data.table)
library(DESeq2)
library(viridis)
library(EnhancedVolcano)
library(VennDiagram)

asw_dds <- readRDS("output/deseq2/asw_all_dds.rds")
##filter out genes with less than 10 reads mapping to them - helpful when transcriptome is large
keep <- rowSums(counts(asw_dds)) >= 10
asw_dds_fil <- asw_dds[keep,]

asw_dds_fil$pc1_sign <- factor(paste(asw_dds_fil$PC1_sign_together))
asw_dds_pc1 <- copy(asw_dds_fil)
design(asw_dds_pc1) <- ~pc1_sign
asw_dds_pc1 <- DESeq(asw_dds_pc1)

res_group <- results(asw_dds_pc1, contrast = c("pc1_sign", "positive", "negative"), lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/pc1_res_group.csv")

##volcano plot w/viridis colours
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##sig with annots
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
trinotate <- fread("data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/pc1_sig_annots.csv")

