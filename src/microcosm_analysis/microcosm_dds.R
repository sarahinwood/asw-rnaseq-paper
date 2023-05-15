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

microcosm_dds <- dds[,dds$Experiment=="Evasion"]
saveRDS(microcosm_dds, snakemake@output[["microcosm_dds"]])

dunedin_microcosm_dds <- dds[,dds$Location=="Dunedin"]
saveRDS(dunedin_microcosm_dds, snakemake@output[["dunedin_microcosm_dds"]])

ruakura_microcosm_dds <- dds[,dds$Location=="Ruakura"]
saveRDS(ruakura_microcosm_dds, snakemake@output[["ruakura_microcosm_dds"]])


# log
sessionInfo()