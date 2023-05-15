library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)
library(tidyverse)
library(readxl)

asw_dds <- readRDS("output/triple_deseq2/asw_dds.rds")

##different factors for analysis
asw_dds$PC1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign_together))
asw_dds$RQN <- factor(paste(asw_dds$RQN))
asw_dds$concentration <- factor(paste(asw_dds$conc))
asw_dds$PC1_experiment <- factor(paste(asw_dds$PC1_sign, asw_dds$experiment, sep=", "))
asw_dds$parasitism <- factor(paste(asw_dds$parasitism))
asw_dds$location <- factor(paste(asw_dds$location))
asw_dds$attack <- factor(paste(asw_dds$attack))

##transformation
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(asw_vst, intgroup=c("location"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=location))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Location")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

## which genes contribute to PCA
ntop=500
rv <- rowVars(assay(asw_vst))
#returns row numbers for top 500 genes
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# top 500 gene names
gene_names <- data.frame(rownames(asw_vst))
gene_names$row_num <- seq.int(nrow(gene_names))
top500_genes <- subset(gene_names, row_num %in% select)
# top 500 rv values
rv_dt <- data.frame(rv)
rv_dt$row_num <- seq.int(nrow(rv_dt))
top500_rv <- subset(rv_dt, row_num %in% select)
#merge
top500 <- merge(top500_genes, top500_rv, by="row_num")
#drop row_num
top500 <- top500[,-1]
# DEG list
sex_pc1 <- read_excel("output/deseq2/pc1_all/pc1_sig_annots.xlsx")
sex_pc1$edited_rn <- paste("ASW", sex_pc1$rn, sep="_")
top_500_annots <- merge(sex_pc1, top500, by.x="edited_rn", by.y="rownames.asw_vst.")
fwrite(top_500_annots, "output/deseq2/pc1_all/pca_top_500_annots.csv")
##470 in sex DEGs - enrichment within this list? - clusterprofiler

##PCA with concentration and RQN
asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("RQN"), returnData=TRUE))

asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.character)
asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.numeric)

asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.character)
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.numeric)

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=RQN))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis()+
  labs(colour="RQN")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_classic()

## plot multiple PCs from analysis ##

asw_dds_location <- readRDS("output/deseq2/asw_location_all/location_dds.rds")
##pull out results from dds_group
res_group <- results(asw_dds_location, lfcThreshold = 1, alpha = 0.05)
##only keep genes where padj is not NA as prcomp can't deal with them
kept_genes <- rownames(subset(res_group, !is.na(padj)))
vst_asssay<- assay(asw_vst)[kept_genes,]
##perform PCA on vst data matrix
pc <- prcomp(t(vst_asssay), center = TRUE, scale = TRUE)
##generate data table of results
pc_wide <- data.table(pc$x, keep.rownames = TRUE)
pc_pd <- melt(pc_wide)

##merge with sample info
sample_data <- fread("data/sample_table.csv")
pc_pd_sample_data <- merge(pc_pd, sample_data, by.x="rn", by.y="sample_name")

##reduce table to first 5 PCs
pc1to6 <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
pc6 <- subset(pc_pd_sample_data, variable %in% pc1to6)

##plot first 5 pcs
ggplot(pc6, aes(x=rn, y=value, colour=attack))+
  scale_colour_viridis(discrete=TRUE)+
  geom_point()+
  theme_bw()+
  theme(axis.text.x=element_blank())+
  facet_wrap(~variable)

