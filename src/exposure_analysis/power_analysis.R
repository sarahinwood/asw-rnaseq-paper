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

## location
location_res <- power_res
loc_disp <- est_count_dispersion(counts, group=c(rep(0,12), rep(1,12)))
location_res$LFC1 <- est_power_distribution(n=24, f=0.05, rho=1, distributionObject=loc_disp, repNumber=10000)
location_res$LFC2 <- est_power_distribution(n=24, f=0.05, rho=2, distributionObject=loc_disp, repNumber=10000)
location_res$LFC5 <- est_power_distribution(n=24, f=0.05, rho=5, distributionObject=loc_disp, repNumber=10000)
location_res$LFC10 <- est_power_distribution(n=24, f=0.05, rho=10, distributionObject=loc_disp, repNumber=10000)
rownames(location_res) <- "location"
location_res <- data.table(location_res, keep.rownames=T)

## exposure
ex_counts <- counts[,c(1:6, 13:18, 7:12, 19:24)]
exposure_res <- power_res
ex_disp <- est_count_dispersion(ex_counts, group=c(rep(0, 12), rep(1, 12)))
exposure_res$LFC1 <- est_power_distribution(n=24, f=0.05, rho=1, distributionObject=ex_disp, repNumber=10000)
exposure_res$LFC2 <- est_power_distribution(n=24, f=0.05, rho=2, distributionObject=ex_disp, repNumber=10000)
exposure_res$LFC5 <- est_power_distribution(n=24, f=0.05, rho=5, distributionObject=ex_disp, repNumber=10000)
exposure_res$LFC10 <- est_power_distribution(n=24, f=0.05, rho=10, distributionObject=ex_disp, repNumber=10000)
rownames(exposure_res) <- "exposure"
exposure_res <- data.table(exposure_res, keep.rownames=T)

## interaction - ru ex vs rest
ruex_counts <- counts[,c(13:18, 1:12, 19:24)]
ruex_res <- power_res
ruex_disp <- est_count_dispersion(ruex_counts, group=c(rep(0, 6), rep(1, 18)))
ruex_res$LFC1 <- est_power_distribution(n=24, f=0.05, rho=1, distributionObject=ruex_disp, repNumber=10000)
ruex_res$LFC2 <- est_power_distribution(n=24, f=0.05, rho=2, distributionObject=ruex_disp, repNumber=10000)
ruex_res$LFC5 <- est_power_distribution(n=24, f=0.05, rho=5, distributionObject=ruex_disp, repNumber=10000)
ruex_res$LFC10 <- est_power_distribution(n=24, f=0.05, rho=10, distributionObject=ruex_disp, repNumber=10000)
rownames(ruex_res) <- "ruex"
ruex_res <- data.table(ruex_res, keep.rownames=T)

## interaction - inv ex vs rest
invex_res <- power_res
invex_disp <- est_count_dispersion(counts, group=c(rep(0, 6), rep(1, 18)))
invex_res$LFC1 <- est_power_distribution(n=24, f=0.05, rho=1, distributionObject=invex_disp, repNumber=10000)
invex_res$LFC2 <- est_power_distribution(n=24, f=0.05, rho=2, distributionObject=invex_disp, repNumber=10000)
invex_res$LFC5 <- est_power_distribution(n=24, f=0.05, rho=5, distributionObject=invex_disp, repNumber=10000)
invex_res$LFC10 <- est_power_distribution(n=24, f=0.05, rho=10, distributionObject=invex_disp, repNumber=10000)
rownames(invex_res) <- "invex"
invex_res <- data.table(invex_res, keep.rownames=T)

full_power_table <- full_join(location_res, full_join(exposure_res, full_join(ruex_res, invex_res)))
full_power_table$LFC1 <- 100*full_power_table$LFC1
full_power_table$LFC2 <- 100*full_power_table$LFC2
full_power_table$LFC5 <- 100*full_power_table$LFC5
full_power_table$LFC10 <- 100*full_power_table$LFC10

fwrite(full_power_table, snakemake@output[["power_res"]])

# write log
sessionInfo()