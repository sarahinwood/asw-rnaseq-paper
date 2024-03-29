#!/usr/bin/env Rscript

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

library(DESeq2)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

exposure_dds <- dds[,dds$Experiment=="Exposure"]

saveRDS(exposure_dds, snakemake@output[["exposure_dds"]])

# log
sessionInfo()