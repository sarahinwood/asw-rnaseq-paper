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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

##different factors for analysis
dds$pc1_sign <- factor(paste(dds$PC1_sign_together))
dds$RQN <- factor(paste(dds$RQN))
dds$concentration <- factor(paste(dds$conc))
dds$Parasitism <- factor(paste(dds$Parasitism))
dds$Location <- factor(paste(dds$Location))
dds$Attack <- factor(paste(dds$Attack))

##transformation
vst <- varianceStabilizingTransformation(dds, blind=TRUE)

## para & Location with 3 diff aspect ratios
pca_plot1 <- plotPCA(vst, intgroup=c("Location", "Parasitism"), returnData=TRUE)
percentVar1 <- round(100 * attr(pca_plot1, "percentVar")) 
pdf(snakemake@output[["PCA_paraloc_1"]])
ggplot(pca_plot1, aes(x=PC1, y=PC2, shape=Location, colour=Parasitism))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Location", colour="Parasitism status")+
  xlab(paste("PC1:", percentVar1[1], "% variance")) + 
  ylab(paste("PC2:", percentVar1[2], "% variance")) + 
  coord_fixed(ratio=1)+
  theme_bw()
dev.off()

pdf(snakemake@output[["PCA_paraloc_2"]])
ggplot(pca_plot1, aes(x=PC1, y=PC2, shape=Location, colour=Parasitism))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Location", colour="Parasitism status")+
  xlab(paste("PC1:", percentVar1[1], "% variance")) + 
  ylab(paste("PC2:", percentVar1[2], "% variance")) + 
  coord_fixed(ratio=2)+
  theme_bw()
dev.off()

pdf(snakemake@output[["PCA_paraloc_3"]])
ggplot(pca_plot1, aes(x=PC1, y=PC2, shape=Location, colour=Parasitism))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Location", colour="Parasitism status")+
  xlab(paste("PC1:", percentVar1[1], "% variance")) + 
  ylab(paste("PC2:", percentVar1[2], "% variance")) + 
  coord_fixed(ratio=3)+
  theme_bw()
dev.off()

## Attack & Location
pca_plot2 <- plotPCA(vst, intgroup=c("Location", "Attack"), returnData=TRUE)
percentVar2 <- round(100 * attr(pca_plot2, "percentVar")) 
pdf(snakemake@output[["PCA_attloc"]])
ggplot(pca_plot2, aes(x=PC1, y=PC2, color=Location, shape=Attack))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Location", shape = "Attack status")+
  xlab(paste("PC1:", percentVar2[1], "% variance")) + 
  ylab(paste("PC2:", percentVar2[2], "% variance")) + 
  coord_fixed(ratio=2)+
  theme_bw()
dev.off()

## concentration
pca_plot3 <- plotPCA(vst, intgroup=c("concentration"), returnData=TRUE)
percentVar3 <- round(100 * attr(pca_plot3, "percentVar")) 
pdf(snakemake@output[["PCA_concentration"]])
ggplot(pca_plot3, aes(x=PC1, y=PC2, color=concentration))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Concentration")+
  xlab(paste("PC1:", percentVar3[1], "% variance")) + 
  ylab(paste("PC2:", percentVar3[2], "% variance")) + 
  coord_fixed(ratio=2)+
  theme_bw()
dev.off()

## RQN
pca_plot4 <- plotPCA(vst, intgroup=c("RQN"), returnData=TRUE)
percentVar4 <- round(100 * attr(pca_plot4, "percentVar")) 
pdf(snakemake@output[["PCA_RQN"]])
ggplot(pca_plot4, aes(x=PC1, y=PC2, color=RQN))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="RQN")+
  xlab(paste("PC1:", percentVar4[1], "% variance")) + 
  ylab(paste("PC2:", percentVar4[2], "% variance")) + 
  coord_fixed(ratio=2)+
  theme_bw()
dev.off()

## which genes contribute to PCA
rv <- rowVars(assay(vst))
rv_dt <- data.frame(rv)
setorder
rv_dt$row_num <- seq.int(nrow(rv_dt))
# link to gene names
gene_names <- data.frame(rownames(vst))
gene_names$row_num <- seq.int(nrow(gene_names))
rv_genes <- merge(rv_dt, gene_names, by="row_num")
rv_genes <- rv_genes[,-1]
setorder(rv_genes, -rv)
fwrite(rv_genes, snakemake@output[["PCA_weightings"]])
## take top 500
top500 <- rv_genes[seq_len(min(500, length(rv_genes)))]
fwrite(top500, snakemake@output[["PCA_top500"]])

# write log
sessionInfo()