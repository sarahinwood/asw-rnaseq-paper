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
loc_disp <- est_count_dispersion(counts, group=c(rep(0,49), rep(1,30)))
location_res$LFC1 <- est_power_distribution(n=79, f=0.05, rho=1, distributionObject=loc_disp, repNumber=10000)
location_res$LFC2 <- est_power_distribution(n=79, f=0.05, rho=2, distributionObject=loc_disp, repNumber=10000)
location_res$LFC5 <- est_power_distribution(n=79, f=0.05, rho=5, distributionObject=loc_disp, repNumber=10000)
location_res$LFC10 <- est_power_distribution(n=79, f=0.05, rho=10, distributionObject=loc_disp, repNumber=10000)
rownames(location_res) <- "location"
location_res <- data.table(location_res, keep.rownames=T)

# attack
attack_res <- power_res
att_counts <- counts[,c(1:22,24,27,28,32,37,38,41, # attacked
                        42:61, 23,25,26,29:31,33:36,39,40, 62:79)] # not
att_disp <- est_count_dispersion(att_counts, group=c(rep(0,41), rep(1,38)))
attack_res$LFC1 <- est_power_distribution(n=79, f=0.05, rho=1, distributionObject=att_disp, repNumber=10000)
attack_res$LFC2 <- est_power_distribution(n=79, f=0.05, rho=2, distributionObject=att_disp, repNumber=10000)
attack_res$LFC5 <- est_power_distribution(n=79, f=0.05, rho=5, distributionObject=att_disp, repNumber=10000)
attack_res$LFC10 <- est_power_distribution(n=79, f=0.05, rho=10, distributionObject=att_disp, repNumber=10000)
rownames(attack_res) <- "attack"
attack_res <- data.table(attack_res, keep.rownames=T)

# Ru attack vs others
ruatt_res <- power_res
ruatt_counts <- counts[,c(42:61, # ru attack
                          1:41, # all inv
                          62:79)] # ru no attack
ruatt_disp <- est_count_dispersion(ruatt_counts, group=c(rep(0,20), rep(1,59)))
ruatt_res$LFC1 <- est_power_distribution(n=79, f=0.05, rho=1, distributionObject=ruatt_disp, repNumber=10000)
ruatt_res$LFC2 <- est_power_distribution(n=79, f=0.05, rho=2, distributionObject=ruatt_disp, repNumber=10000)
ruatt_res$LFC5 <- est_power_distribution(n=79, f=0.05, rho=5, distributionObject=ruatt_disp, repNumber=10000)
ruatt_res$LFC10 <- est_power_distribution(n=79, f=0.05, rho=10, distributionObject=ruatt_disp, repNumber=10000)
rownames(ruatt_res) <- "ru attack"
ruatt_res <- data.table(ruatt_res, keep.rownames=T)

# Inv attack vs others
invatt_res <- power_res
invatt_counts <- counts[,c(1:22,24,27,28,32,37,38,41, # inv attack
                           23,25,26,29:31,33:36,39,40,42:79)] # others
invatt_disp <- est_count_dispersion(invatt_counts, group=c(rep(0,29), rep(1,50)))
invatt_res$LFC1 <- est_power_distribution(n=79, f=0.05, rho=1, distributionObject=invatt_disp, repNumber=10000)
invatt_res$LFC2 <- est_power_distribution(n=79, f=0.05, rho=2, distributionObject=invatt_disp, repNumber=10000)
invatt_res$LFC5 <- est_power_distribution(n=79, f=0.05, rho=5, distributionObject=invatt_disp, repNumber=10000)
invatt_res$LFC10 <- est_power_distribution(n=79, f=0.05, rho=10, distributionObject=invatt_disp, repNumber=10000)
rownames(invatt_res) <- "inv attack"
invatt_res <- data.table(invatt_res, keep.rownames=T)

# parasitism
para_res <- power_res
para_counts <- counts[,c(1,2,6,12,14,16,17,19,21,22,24,27,28,32,37,38,41,42,45,50,60, # para
                         3:5,7:11,13,15,18,20,23,25,26,29:31,33:36,39,40,43,44,46:49,51:59,61:79)] # unpara
para_disp <- est_count_dispersion(para_counts, group=c(rep(0,21), rep(1,58)))
para_res$LFC1 <- est_power_distribution(n=79, f=0.05, rho=1, distributionObject=para_disp, repNumber=10000)
para_res$LFC2 <- est_power_distribution(n=79, f=0.05, rho=2, distributionObject=para_disp, repNumber=10000)
para_res$LFC5 <- est_power_distribution(n=79, f=0.05, rho=5, distributionObject=para_disp, repNumber=10000)
para_res$LFC10 <- est_power_distribution(n=79, f=0.05, rho=10, distributionObject=para_disp, repNumber=10000)
rownames(para_res) <- "parasitism"
para_res <- data.table(para_res, keep.rownames=T)

full_power_table <- full_join(location_res, full_join(attack_res, full_join(ruatt_res, full_join(invatt_res, para_res))))
full_power_table$LFC1 <- 100*full_power_table$LFC1
full_power_table$LFC2 <- 100*full_power_table$LFC2
full_power_table$LFC5 <- 100*full_power_table$LFC5
full_power_table$LFC10 <- 100*full_power_table$LFC10

fwrite(full_power_table, snakemake@output[["power_res"]])

# write log
sessionInfo()