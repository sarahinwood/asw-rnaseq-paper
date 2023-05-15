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

genes_file <- snakemake@input[["genes_file"]]
trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

genes <- fread(genes_file)
trinotate <- fread(trinotate_file)

genes_annots <- merge(genes, trinotate, by.x="rownames.vst.", by.y="#gene_id", all.x=T, all.y=F)
fwrite(genes_annots, snakemake@output[["genes_annots"]])

# write log
sessionInfo()