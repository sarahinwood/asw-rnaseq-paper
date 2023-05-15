library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

#########
## ASW ##
#########

asw_dds <- readRDS("output/evasion_deseq2/asw_evasion_dds.rds")
asw_dds$para_pcr <- factor(paste(asw_dds$parasitism))
asw_dds$RQN <- factor(paste(asw_dds$RQN))
asw_dds$concentration <- factor(paste(asw_dds$conc))
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(asw_vst, intgroup=c("location"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar"))

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=location))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="ASW source\nlocation")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_classic()


##PCA with concentration and RQN
asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("concentration"), returnData=TRUE))

asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.character)
asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.numeric)

asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.character)
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.numeric)

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=concentration))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis()+
  labs(colour="RNA\nconcentration\n(ng/uL)")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_classic()
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
