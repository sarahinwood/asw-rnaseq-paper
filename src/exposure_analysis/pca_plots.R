library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

dds <- readRDS("output/exposed_deseq2/exposed_dds.rds")
dds$RQN <- factor(paste(dds$RQN))
dds$concentration <- factor(paste(dds$conc))
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(vst, intgroup=c("parasitism"))

pca_table <- data.table(plotPCA(vst, intgroup=c("RQN"), returnData=TRUE))

pca_table$concentration <- sapply(pca_table$concentration, as.character)
pca_table$concentration <- sapply(pca_table$concentration, as.numeric)

pca_table$RIN <- sapply(pca_table$RQN, as.character)
pca_table$RIN <- sapply(pca_table$RQN, as.numeric)


ggplot(pca_table, aes(x=PC1, y=PC2, color=RQN))+
  geom_point(alpha=0.5, size=3)+
  scale_color_viridis()




