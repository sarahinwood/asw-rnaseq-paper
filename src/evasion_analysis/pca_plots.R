library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

#########
## ASW ##
#########

asw_dds <- readRDS("output/evasion_deseq2/asw_evasion_dds.rds")
asw_dds$para_pcr <- factor(paste(asw_dds$parasitism))
asw_dds$RIN <- factor(paste(asw_dds$RQN))
asw_dds$concentration <- factor(paste(asw_dds$conc))
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("parasitism"))

asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("concentration"), returnData=TRUE))

asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.numeric)

asw_pca_table$RIN <- sapply(asw_pca_table$RIN, as.numeric)
##remove sample where RQN=0
asw_pca_table <- asw_pca_table[-58,]

ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=concentration))+
  geom_point(alpha=0.5, size=3)+
  scale_color_viridis()

########
## Mh ##
########

mh_dds <- readRDS("output/evasion_deseq2/mh_evasion_dds.rds")
mh_dds$para_pcr <- factor(paste(mh_dds$parasitism))
mh_dds$RIN <- factor(paste(mh_dds$RQN))
mh_dds$concentration <- factor(paste(mh_dds$conc))
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(mh_vst, intgroup=c("parasitism"))

mh_pcr_counts <- fread("output/evasion_deseq2/mh_pcr_counts.csv")

mh_pca_table <- data.table(plotPCA(mh_vst, intgroup=c("parasitism"), returnData=TRUE))
mh_pca_plot <- merge(mh_pca_table, mh_pcr_counts, by.x="name", by.y="rn")

#and para because of high outliers for para
ggplot(mh_pca_plot, aes(x=PC1, y=PC2, color=count))+
  geom_point(alpha=0.5, size=3)+
  scale_color_viridis()
