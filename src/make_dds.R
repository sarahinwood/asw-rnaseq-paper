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

library(tximport)
library(data.table)
library(DESeq2)

###########
# GLOBALS #
###########

gene2tx_file <- snakemake@input[["gene2tx_file"]]
sample_data_file <- snakemake@input[["sample_data_file"]]

########
# MAIN #
########

gene2tx <- fread(gene2tx_file, header = T)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_MhV1_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)

##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##make table of TPMs
txi_abundance <- data.table(txi$abundance, keep.rownames=TRUE)
fwrite(txi_abundance, snakemake@output[["salmon_tpm"]])

##Import table describing samples
sample_data <- fread(sample_data_file, header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##estimate size factors on whole dataset
dds <- estimateSizeFactors(dds)
##save whole dds object
saveRDS(dds, snakemake@output[["mh_asw_MhV1_dds"]])

################################################
##subset gene table into ASW, Mh & MhV1 genes ##
################################################

##asw gene list
asw_tx <- subset(tx2gene, grepl("ASW_", V1))
asw_gene <- unique(asw_tx$V1)
##subset asw dds
asw_dds <- dds[asw_gene,]
saveRDS(asw_dds, snakemake@output[["asw_dds"]])

##mh gene list
mh_tx <- subset(tx2gene, grepl("MH_", V1))
mh_gene <- unique(mh_tx$V1)
##subset mh dds
mh_dds <- dds[mh_gene,]
saveRDS(mh_dds, snakemake@output[["mh_dds"]])

##MhV1 gene list
mhv1_tx <- subset(tx2gene, grepl("ORF", V1))
mhv1_gene <- unique(mhv1_tx$V1)
##subset MhV dds
mhv1_dds <- dds[mhv1_gene,]
saveRDS(mhv1_dds, snakemake@output[["MhV1_dds"]])

# write log
sessionInfo()