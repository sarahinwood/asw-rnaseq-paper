library(data.table)
library(pheatmap)
library(tidyverse)
library(DESeq2)
library(viridis)

dds_parasitism <- readRDS("output/triple_deseq2/Mh/microcosm_samples/all_microcosm_parasitism/para_dds.rds")
mh_venom_degs <- fread("output/triple_deseq2/Mh/microcosm_samples/all_microcosm_parasitism/Mh_degs_venom_paper_annots.csv")
  
##vst transform
vst <- varianceStabilizingTransformation(dds_parasitism, blind=TRUE)
vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)
##subset for DEGs
vst_degs <- subset(vst_assay_dt, rn %in% mh_venom_degs$`Gene ID`)
##turn first row back to row name
vst_degs <- vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
vst_degs_plot <- vst_degs[,c(1,2,6,12,14,16,17,19,21,22,24,27,28,32,37,38,41, # Dun para
                             42,45,50,60, #Ru para
                             3:5,7:10,13,15,18,20,23,25,26,29:31,33:36,39,40, # Dun unpara
                             43,44,46:49,51:59, 61:79)] #ru unpara
##get location label info
sample_to_label <- data.table(data.frame(colData(dds_parasitism)[,c("Location", "Parasitism", "sample_name")]))
sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Location = c(Dunedin="#fcfdbf", Ruakura="#fc8961"),
                         Parasitism=c(Parasitised = "#b73779", Undetected="#51127c"))

##plot
##not clustered by sample
pdf("output/triple_deseq2/Mh/microcosm_samples/all_microcosm_parasitism/Mh_venom_unclustered.pdf")
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

pdf("output/triple_deseq2/Mh/microcosm_samples/all_microcosm_parasitism/Mh_venom_clustered.pdf")
pheatmap(vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_label, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

# expression in unpara samples
venom_unpara <- vst_degs[,c(3:5,7:10,13,15,18,20,23,25,26,29:31,33:36,39,40, # Dun unpara
                             43,44,46:49,51:59, 61:79)]  #ru unpara

#venom gene TPMs
all_tpms <- fread("output/tpms/sample_TPMs.csv")
Mh_tpms <- subset(all_tpms, grepl("MH_TRINITY", all_tpms$rn))
Mh_venom_tpms <- subset(Mh_tpms, rn %in% mh_venom_degs$`Gene ID`)
Mh_venom_tpms_t <- data.table(t(Mh_venom_tpms), keep.rownames=T)
colnames(Mh_venom_tpms_t) <- as.character(Mh_venom_tpms_t[1,])
Mh_venom_tpms_t <- Mh_venom_tpms_t[-1,]

# subset for only unpara
sample_table <- fread("data/sample_table.csv")
unpara <- subset(sample_table, Parasitism=="Undetected")
unpara_names <- unpara$sample_name
Mh_venom_tpms_unpara <- subset(Mh_venom_tpms_t, rn %in% unpara_names)
fwrite(Mh_venom_tpms_unpara, "output/triple_deseq2/Mh/microcosm_samples/all_microcosm_parasitism/Mh_venom_tpms_unpara.csv")




