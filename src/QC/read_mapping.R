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
library(DESeq2)
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

asw_dds_file <- snakemake@input[["asw_dds_file"]]
mh_dds_file <- snakemake@input[["mh_dds_file"]]
MhV1_dds_file <- snakemake@input[["MhV1_dds_file"]]

########
# MAIN #
########

## ASW reads mapped
asw_dds <- readRDS(asw_dds_file)
counts_table_asw <- (data.table(counts(asw_dds)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("sample_name", "readpairs_mapped_ASW"))

## Mh reads mapped
mh_dds <- readRDS(mh_dds_file)
counts_table_mh <- (data.table(counts(mh_dds)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("sample_name", "readpairs_mapped_MH"))

## MhV1 reads mapped
MhV1_dds <- readRDS(MhV1_dds_file)
counts_table_mhv1 <- (data.table(counts(MhV1_dds)))
counts_colSums_mhv1 <- setDT(data.frame(colSums(counts_table_mhv1, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mhv1, old=c("rn", "colSums.counts_table_mhv1..na.rm...TRUE."), new=c("sample_name", "readpairs_mapped_MhV1"))

##all reads mapped - hoping multiqc will collate a good table for me
read_mapping <- merge(counts_colSums_asw, merge(counts_colSums_mh, counts_colSums_mhv1, by="sample_name"), by="sample_name")
read_mapping$total_mapped_readpairs <- (read_mapping$readpairs_mapped_ASW+read_mapping$readpairs_mapped_MH+read_mapping$readpairs_mapped_MhV1)

##mapping %s
read_mapping$`%_ofmapped_ASW` <- (read_mapping$readpairs_mapped_ASW/read_mapping$total_mapped_readpairs)*100
read_mapping$`%_ofmapped_MH` <- (read_mapping$readpairs_mapped_MH/read_mapping$total_mapped_readpairs)*100
read_mapping$`%_ofmapped_MhV1` <- (read_mapping$readpairs_mapped_MhV1/read_mapping$total_mapped_readpairs)*100
fwrite(read_mapping, snakemake@output[["read_mapping_table"]])

## merge with para status
coldata <- data.table(data.frame(colData(asw_dds)), keep.rownames=T)
sample_para <- coldata[,c(2,9,10,11)]
read_mapping_para_status <- merge(read_mapping, sample_para)
percent_map_para <- read_mapping_para_status[,c(6,7,8,10)]
long <- melt(percent_map_para)
summary <- aggregate(value ~ Parasitism+variable, long, mean)
fwrite(summary, snakemake@output[["read_mapping_summary"]])

## all stacked bar - Parasitism
pdf(snakemake@output[["all_stacked"]])
ggplot(summary, aes(fill=variable, y=value, x=Parasitism))+
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_y_continuous(name="Percentage of mapped reads", breaks = scales::pretty_breaks(n = 10))+
  scale_fill_viridis(discrete=T, name="Species", labels=c("ASW", expression(italic("M. hyperodae")), "MhV1"))+
  theme(legend.text.align = 0)+
  xlab("Parasitism status")
dev.off()

## ASW mapping - Parasitism
pdf(snakemake@output[["asw_box"]])
ggplot(read_mapping_para_status, aes(x=Parasitism, y=`%_ofmapped_ASW`))+
  geom_boxplot(aes(fill=Parasitism, colour=Parasitism), alpha=0.7, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("Parasitism status")+
  ylab(expression(paste("% of normalised counts to ASW transcriptome")))
dev.off()

## Mh mapping - Parasitism
pdf(snakemake@output[["mh_box"]])
ggplot(read_mapping_para_status, aes(x=Parasitism, y=`%_ofmapped_MH`))+
  geom_boxplot(aes(fill=Parasitism, colour=Parasitism), alpha=0.7, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("Parasitism status")+
  ylab(expression(paste("% of normalised counts to ", italic("M. hyperodae"), " transcriptome")))
dev.off()

## MhV1 mapping - Parasitism
pdf(snakemake@output[["MhV1_box"]])
ggplot(read_mapping_para_status, aes(x=Parasitism, y=`%_ofmapped_MhV1`))+
  geom_boxplot(aes(fill=Parasitism, colour=Parasitism), alpha=0.7, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("Parasitism status")+
  ylab(expression(paste("% of normalised counts to MhV1")))
dev.off()

## Mh parasitism PCR target gene plot counts
Mh_PCR_counts <- plotCounts(mh_dds, "MH_TRINITY_DN7604_c0_g1", intgroup = c("Parasitism_PCR_result"), returnData = TRUE)
pdf(snakemake@output[["Mh_PCR_target"]])
ggplot(Mh_PCR_counts, aes(x=Parasitism_PCR_result, y=count, colour=Parasitism_PCR_result))+
  geom_point(size=3, alpha=0.7)+
  labs(y="Normalised counts", x="Parasitism PCR result")+
  scale_colour_viridis(discrete=TRUE)+
  scale_y_continuous(trans="log10")+
  theme_bw()+
  theme(legend.position="none")
dev.off()

##t test for ASW mapping between para and undetected
asw_ttest_para <- t.test(`%_ofmapped_ASW` ~ Parasitism, data=read_mapping_para_status)
paste("ASW mapping by parasitism status", asw_ttest_para)
##t test for Mh mapping between para and undetected
mh_ttest_para <- t.test(`%_ofmapped_MH` ~ Parasitism, data=read_mapping_para_status)
paste("Mh mapping between parasitism", mh_ttest_para)
##t test for MhV1 mapping between para and undetected
MhV1_ttest_para <- t.test(`%_ofmapped_MhV1` ~ Parasitism, data=read_mapping_para_status)
paste("MhV1 mapping by parasitism status", MhV1_ttest_para)

##t test for Mh mapping between Inv para and undetected
inv_para_all_unpara <- subset(read_mapping_para_status, !(Parasitism=="Parasitised"&Location=="Ruakura"))
inv_para_mh <- t.test(`%_ofmapped_MH` ~ Parasitism, data=inv_para_all_unpara)
paste("Inv ASW mapping by parasitism status", inv_para_mh)

#t test for mapping rates of parasitised ASW on location
para_mapping <- subset(read_mapping_para_status, Parasitism=="Parasitised")
para_loc_mh <- t.test(`%_ofmapped_MH` ~Location, data=para_mapping)
paste("parasitised ASW mapping rates by location", para_loc_mh)

# write log
sessionInfo()

