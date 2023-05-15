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
library(tidyverse)
library(ggplot2)
library(viridis)
library(cowplot)

###########
# GLOBALS #
###########

##DEG_list             
pfam_file <- snakemake@input[["pfam_file"]]
go_file <- snakemake@input[["go_file"]]

########
# MAIN #
########

pfam_enrich <- fread(pfam_file)
pfam_enrich$mod <- ifelse(pfam_enrich$ID=="PF13499.6", "7",
                              ifelse(pfam_enrich$ID=="PF13833.6", "8",
                                     ifelse(pfam_enrich$ID=="PF00036.32", "1",
                                            ifelse(pfam_enrich$ID=="PF13202.6", "5", ""))))
pfam_enrich$Description <- paste(pfam_enrich$Description, pfam_enrich$mod, sep=" ")
pfam_enrich$Description <- factor(pfam_enrich$Description, levels=pfam_enrich$Description[order(pfam_enrich$GeneRatio, pfam_enrich$Description, decreasing=F)])
pfam_enrich$GeneRatio_percent <- pfam_enrich$GeneRatio*100

pfam_enrich$BG_DEG_no <- tstrsplit(pfam_enrich$BgRatio, "/", keep=1)
pfam_enrich$BG_total_no <- tstrsplit(pfam_enrich$BgRatio, "/", keep=2)
pfam_enrich$BgRatio <- as.numeric(pfam_enrich$BG_DEG_no)/as.numeric(pfam_enrich$BG_total_no)
pfam_enrich$DEG_to_BG_ratio <- pfam_enrich$GeneRatio/pfam_enrich$BgRatio
pfam_enrich$log_padj <- -log(pfam_enrich$p.adjust)
pfam_enrich$colour <- paste("#440154FF")
pfam_enrich$Kind <- paste("Pfam")
pfam_enrich$plot_label <- paste("Pfam domain")

go_enrich <- fread(go_file)
go_enrich$mod <- ifelse(go_enrich$ID=="GO:0004553", ", hydrolyzing O-glycosyl compounds", "")
go_enrich$Description <- paste(go_enrich$Description, go_enrich$mod, sep="")
go_enrich$Description <- factor(go_enrich$Description, levels=go_enrich$Description[order(go_enrich$pathway_kind, go_enrich$GeneRatio, go_enrich$Description, decreasing=T)])
go_enrich$GeneRatio_percent <- go_enrich$GeneRatio*100

go_enrich$BG_DEG_no <- tstrsplit(go_enrich$BgRatio, "/", keep=1)
go_enrich$BG_total_no <- tstrsplit(go_enrich$BgRatio, "/", keep=2)
go_enrich$BgRatio <- as.numeric(go_enrich$BG_DEG_no)/as.numeric(go_enrich$BG_total_no)
go_enrich$DEG_to_BG_ratio <- go_enrich$GeneRatio/go_enrich$BgRatio
go_enrich$log_padj <- -log(go_enrich$p.adjust)
go_enrich$Kind <- paste("GO")
go_enrich$colour <- ifelse(go_enrich$pathway_kind=="biological process", paste("#21918c"), ifelse(go_enrich$pathway_kind=="cellular component", paste("#5ec962"), paste("#fde725")))
go_enrich$plot_label <- paste("GO", go_enrich$pathway_kind)

full_enrich <- full_join(pfam_enrich, go_enrich)
full_enrich$pathway_kind <- factor(full_enrich$pathway_kind, levels=c("GO biological process", "GO cellular component", "GO molecular function", "Pfam domain"))

pfam_old <- ggplot(pfam_enrich, aes(x=Description, y=GeneRatio_percent)) +
  geom_col(aes(fill="#440154FF"))+
  labs(x="Pfam domain", y="Percentage of DTU genes with associated term", size="Leading\nedge size") +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  theme(legend.position = "none")


go_old <- ggplot(go_enrich, aes(Description, GeneRatio_percent)) +
  geom_col(aes(fill=pathway_kind))+
  labs(x="Gene ontology terms", y="Percentage of DTU genes with associated term",
       fill="GO domain", size="Leading\nedge size") +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE, begin=0.5)+
  theme_bw()#+
  #theme(legend.position = "none",
  #      axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())

pfam <- ggplot(pfam_enrich, aes(Description, log_padj)) +
  geom_point(aes(size=Count, colour=colour))+
  labs(x="Pfam domains", y="-log(adj. p-value)",
       size="Leading\nedge size") +
  coord_flip() +
  scale_size(limits=c(1, 30), range=c(5,10), breaks = c(5,10,15,20,25,30))+
  scale_y_continuous(limits=c(0,10))+
  scale_colour_manual(values=c("#440154FF"="#440154FF"), guide='none')+
  theme_bw()

go <- ggplot(go_enrich, aes(Description, log_padj)) +
  geom_point(aes(colour=pathway_kind, size=Count))+
  labs(x="Gene ontology terms", y="-log(adj. p-value)",
       colour="GO domain", size="Leading\nedge size") +
  coord_flip() +
  scale_colour_viridis(discrete=TRUE, begin=0.5)+
  scale_size(limits=c(1, 30), range=c(5,10), breaks = c(5,10,15,20,25,30))+
  scale_y_continuous(limits=c(0,10))+
  theme_bw()

go_length <- length(go_enrich$Description)
pfam_length <- length(pfam_enrich$Description)
plot_length <- (go_length+pfam_length)/4.272727
go_plot_height <-go_length/pfam_length

pdf(snakemake@output[["enrich_plot"]], width=11, height=plot_length) # default 7x7
plot_grid(go, pfam, ncol=1, align="v",
          rel_heights=c(go_plot_height, 1))
dev.off()

full_plot <- ggplot(full_enrich, aes(Description, log_padj)) +
  geom_point(aes(size=Count, colour=plot_label))+
  labs(x="Pfam domains", y="-log(adj. p-value)",
       colour="Enriched term category",
       size="No. genes \ndriving enrichment") +
  coord_flip() +
  scale_size(limits=c(1, 30), range=c(5,10), breaks = c(5,10,15,20,25,30))+
  scale_y_continuous(limits=c(0,10))+
  scale_colour_viridis(discrete=T, direction=-1)+
  theme_bw()


pdf(snakemake@output[["full_plot"]], width=11, height=3) # default 7x7
full_plot
dev.off()

# write log
sessionInfo()