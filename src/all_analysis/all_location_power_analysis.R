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

library(RnaSeqSampleSize)
library(DESeq2)
library(tidyverse)
library(data.table)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

asw_dds <- readRDS(dds_file)
counts <- data.frame(counts(asw_dds, normalized=F))

# make table to fill with res
power_res <- data.frame(matrix(ncol = 4, nrow = 1))
col_names <- c("LFC1", "LFC2", "LFC5", "LFC10")
colnames(power_res) <- col_names

# location
location_res <- power_res
loc_disp <- est_count_dispersion(counts, group=c(rep(0,53), rep(1,50)))
location_res$LFC1 <- est_power_distribution(n=103, f=0.05, rho=1, distributionObject=loc_disp, repNumber=10000)
location_res$LFC2 <- est_power_distribution(n=103, f=0.05, rho=2, distributionObject=loc_disp, repNumber=10000)
location_res$LFC5 <- est_power_distribution(n=103, f=0.05, rho=5, distributionObject=loc_disp, repNumber=10000)
location_res$LFC10 <- est_power_distribution(n=103, f=0.05, rho=10, distributionObject=loc_disp, repNumber=10000)
rownames(location_res) <- "location"

location_res$LFC1 <- 100*location_res$LFC1
location_res$LFC2 <- 100*location_res$LFC2
location_res$LFC5 <- 100*location_res$LFC5
location_res$LFC10 <- 100*location_res$LFC10
location_res <- data.table(location_res, keep.rownames = T)
fwrite(location_res, snakemake@output[["power_res"]])

# write log
sessionInfo()