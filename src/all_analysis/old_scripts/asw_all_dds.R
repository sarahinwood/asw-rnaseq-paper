#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(tximport)
library(data.table)
library(DESeq2)

gene_trans_map <- snakemake@input[['asw_mh_gene_trans_map']]
gene2tx <- fread(gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
all_quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)

##assign names to quant files from folder name
names(all_quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", all_quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(all_quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_table <- snakemake@input[['sample_table']]
sample_data <- fread(sample_table, header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
dds <- estimateSizeFactors(dds)

all_dds <- snakemake@output[['all_dds']]
saveRDS(dds, all_dds)

total_reads_per_sample <- colSums(counts(dds))
table_counts <- data.table(data.frame(total_reads_per_sample), keep.rownames=TRUE)
fwrite(table_counts, "output/deseq2/asw_total_counts.csv")

## subset into ASW and MH ##
##asw gene list
asw_tx <- subset(tx2gene, grepl("ASW_", V1))
asw_gene <- unique(asw_tx$V1)
##mh gene list
mh_tx <- subset(tx2gene, grepl("MH_", V1))
mh_gene <- unique(mh_tx$V1)

##subset asw dds
asw_dds <- dds[asw_gene,]
saveRDS(asw_dds, snakemake@output[["asw_dds"]])
##subset mh dds
mh_dds <- dds[mh_gene,]
saveRDS(mh_dds, snakemake@output[["mh_dds"]])

# log
sessionInfo()


