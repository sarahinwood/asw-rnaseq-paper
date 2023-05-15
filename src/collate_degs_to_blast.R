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
library(seqinr)

###########
# GLOBALS #
###########

asw_deg_dir <- snakemake@params[["asw_deg_dir"]]
deg_file_suffix <- snakemake@params[["deg_file_suffix"]]
fasta <- snakemake@input[["fasta"]]

########
# MAIN #
########

# list of DEGs
deg_files <- list.files(asw_deg_dir, deg_file_suffix, full.names = T, recursive=T)
deg_tables <-  Map(read_csv, deg_files)
deg_tables_dt <- data.table(do.call(plyr::rbind.fill, deg_tables))
degs_list <- unique(deg_tables_dt$transcript_id)

# read and subset fasta
fastafile <- read.fasta(file = fasta, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
degs_to_blast <- fastafile[c(which(names(fastafile) %in% degs_list))]

write.fasta(degs_to_blast, names=names(degs_to_blast), file.out=snakemake@output[["degs_to_blast"]])

# write log
sessionInfo()