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
library(viridis)

###########
# GLOBALS #
###########

tpm_file <- snakemake@input[["tpm_file"]]

########
# MAIN #
########

tpm_sample_table <- fread(tpm_file)

exposure_tpms <- data.table()
exposure_tpms$`#gene_id` <- tpm_sample_table$rn
exposure_tpms$exposure_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9", "RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9")))
exposure_tpms$control_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8", "RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7")))
exposure_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8", "RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7", "IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9", "RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9")))
fwrite(exposure_tpms, snakemake@output[["exposure_tpms"]])

exposure_location_tpms <- data.table()
exposure_location_tpms$`#gene_id` <- tpm_sample_table$rn
exposure_location_tpms$dunedin_exposed_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9")))
exposure_location_tpms$ruakura_exposed_tpm <- rowMeans(subset(tpm_sample_table, select = c("RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9")))
exposure_location_tpms$dunedin_control_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8")))
exposure_location_tpms$ruakura_control_tpm <- rowMeans(subset(tpm_sample_table, select = c("RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7")))
exposure_location_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7", "IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8", "RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9", "IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9")))
fwrite(exposure_location_tpms, snakemake@output[["exposure_location_tpms"]])

location_tpms <- data.table()
location_tpms$`#gene_id` <- tpm_sample_table$rn
location_tpms$ruakura_tpm <- rowMeans(subset(tpm_sample_table, select = c("RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9", "RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7", "RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C", "RN1", "RN8", "RN12", "RN20")))
location_tpms$dunedin_tpm <- rowMeans(subset(tpm_sample_table, select = c("IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9", "IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8", "DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27", "DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
location_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("RALM_Ex_2", "RALM_Ex_5", "RALM_Ex_6", "RALM_Ex_7", "RALM_Ex_8", "RALM_Ex_9", "RALM_NC_2", "RALM_NC_3", "RALM_NC_4", "RALM_NC_5", "RALM_NC_6", "RALM_NC_7", "RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C", "RN1", "RN8", "RN12", "RN20", "IALM_Ex_2", "IALM_Ex_4", "IALM_Ex_5", "IALM_Ex_6", "IALM_Ex_8", "IALM_Ex_9", "IALM_NC_2", "IALM_NC_4", "IALM_NC_5", "IALM_NC_6", "IALM_NC_7", "IALM_NC_8", "DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27", "DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
fwrite(location_tpms, snakemake@output[["location_tpms"]])

parasitism_tpms <- data.table()
parasitism_tpms$`#gene_id` <- tpm_sample_table$rn
parasitism_tpms$parasitised_tpm <- rowMeans(subset(tpm_sample_table, select = c("DY1", "DN1", "RN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "RN8", "DN9", "DN10C", "DY12C", "RN12", "DN14", "DY15C", "DY16", "RN20", "DY25", "DN25C")))
parasitism_tpms$undetected_tpm <- rowMeans(subset(tpm_sample_table, select = c("RY1", "DY2C", "RY2", "RN2C", "DY3", "RY3", "DN3", "RN3", "RN4C", "RY5", "RN5C", "DY6", "DN6", "RN6C", "DY7C", "RY7", "RN7", "RY8C", "DN8C", "RY9", "RN9", "DY10", "RY10C", "RN10", "RY11", "DN11", "RN11C", "RY12C", "DN12", "DY13C", "RY13C", "DN13C", "RN13", "DY14", "RY14C", "RN14", "RY15", "DN15", "RN15", "RY16", "DN16", "RN16", "DY20", "RY20", "DN20C", "DY21C", "RY21", "DN21", "RN21", "DY22", "RY22C", "DN22", "RN22", "RY25", "RN25C", "DY26", "DN26", "DY27")))
parasitism_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("RY1", "DY2C", "RY2", "RN2C", "DY3", "RY3", "DN3", "RN3", "RN4C", "RY5", "RN5C", "DY6", "DN6", "RN6C", "DY7C", "RY7", "RN7", "RY8C", "DN8C", "RY9", "RN9", "DY10", "RY10C", "RN10", "RY11", "DN11", "RN11C", "RY12C", "DN12", "DY13C", "RY13C", "DN13C", "RN13", "DY14", "RY14C", "RN14", "RY15", "DN15", "RN15", "RY16", "DN16", "RN16", "DY20", "RY20", "DN20C", "DY21C", "RY21", "DN21", "RN21", "DY22", "RY22C", "DN22", "RN22", "RY25", "RN25C", "DY26", "DN26", "DY27", "DY1", "DN1", "RN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "RN8", "DN9", "DN10C", "DY12C", "RN12", "DN14", "DY15C", "DY16", "RN20", "DY25", "DN25C")))
fwrite(parasitism_tpms, snakemake@output[["parasitism_tpms"]])

parasitism_location_tpms <- data.table()
parasitism_location_tpms$`#gene_id` <- tpm_sample_table$rn
parasitism_location_tpms$dunedin_parasitised_tpm <- rowMeans(subset(tpm_sample_table, select = c("DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
parasitism_location_tpms$ruakura_parasitised_tpm <- rowMeans(subset(tpm_sample_table, select = c("RN1", "RN8", "RN12", "RN20")))
parasitism_location_tpms$dunedin_undetected_tpm <- rowMeans(subset(tpm_sample_table, select = c("DY2C", "DY3", "DN3", "DY6", "DN6", "DY7C", "DN8C", "DY10", "DN11", "DN12", "DY13C", "DN13C", "DY14", "DN15", "DN16", "DY20", "DN20C", "DY21C", "DN21", "DY22", "DN22", "DY26", "DN26", "DY27")))
parasitism_location_tpms$ruakura_undetected_tpm <- rowMeans(subset(tpm_sample_table, select = c("RY1", "RY2", "RN2C", "RY3", "RN3", "RN4C", "RY5", "RN5C", "RN6C", "RY7", "RN7", "RY8C", "RY9", "RN9", "RY10C", "RN10", "RY11", "RN11C", "RY12C", "RY13C", "RN13", "RY14C", "RN14", "RY15", "RN15", "RY16", "RN16", "RY20", "RY21", "RN21", "RY22C", "RN22", "RY25", "RN25C")))
parasitism_location_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("RY1", "RY2", "RN2C", "RY3", "RN3", "RN4C", "RY5", "RN5C", "RN6C", "RY7", "RN7", "RY8C", "RY9", "RN9", "RY10C", "RN10", "RY11", "RN11C", "RY12C", "RY13C", "RN13", "RY14C", "RN14", "RY15", "RN15", "RY16", "RN16", "RY20", "RY21", "RN21", "RY22C", "RN22", "RY25", "RN25C", "DY2C", "DY3", "DN3", "DY6", "DN6", "DY7C", "DN8C", "DY10", "DN11", "DN12", "DY13C", "DN13C", "DY14", "DN15", "DN16", "DY20", "DN20C", "DY21C", "DN21", "DY22", "DN22", "DY26", "DN26", "DY27", "RN1", "RN8", "RN12", "RN20", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
fwrite(parasitism_location_tpms, snakemake@output[["parasitism_location_tpms"]])

attack_tpms  <- data.table()
attack_tpms$`#gene_id` <- tpm_sample_table$rn
attack_tpms$attacked_tpm <- rowMeans(subset(tpm_sample_table, select = c("DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C", "RN1", "RN8", "RN12", "RN20")))
attack_tpms$not_detected_tpm <- rowMeans(subset(tpm_sample_table, select = c("DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27", "RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25")))
attack_tpms$all_tpm <- rowMeans(subset(tpm_sample_table, select = c("DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27", "RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25", "DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C", "RN1", "RN8", "RN12", "RN20")))
fwrite(attack_tpms, snakemake@output[["attack_tpms"]])

attack_location_tpms <- data.table()
attack_location_tpms$`#gene_id` <- tpm_sample_table$rn
attack_location_tpms$dunedin_attacked <- rowMeans(subset(tpm_sample_table, select = c("DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
attack_location_tpms$ruakura_attacked <- rowMeans(subset(tpm_sample_table, select = c("RN1", "RN8", "RN12", "RN20", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C")))
attack_location_tpms$dunedin_undetected <- rowMeans(subset(tpm_sample_table, select = c("DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27")))
attack_location_tpms$ruakura_undetected <- rowMeans(subset(tpm_sample_table, select = c("RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25")))
attack_location_tpms$all_undetected <- rowMeans(subset(tpm_sample_table, select = c("RY1", "RY2", "RY3", "RY5", "RY7", "RY8C", "RY9", "RY10C", "RY11", "RY12C", "RY13C", "RY14C", "RY15", "RY16", "RY20", "RY21", "RY22C", "RY25", "DY2C", "DY3", "DY6", "DY7C", "DY10", "DY13C", "DY14", "DY20", "DY21C", "DY22", "DY26", "DY27", "RN1", "RN8", "RN12", "RN20", "RN2C", "RN3", "RN4C", "RN5C", "RN6C", "RN7", "RN9", "RN10", "RN11C", "RN13", "RN14", "RN15", "RN16", "RN21", "RN22", "RN25C", "DN3", "DN6", "DN8C", "DN11", "DN12", "DN13C", "DN15", "DN16", "DN20C", "DN21", "DN22", "DN26", "DY1", "DN1", "DN2C", "DY4", "DN4", "DY5", "DN5C", "DN7", "DY8", "DN9", "DN10C", "DY12C", "DN14", "DY15C", "DY16", "DY25", "DN25C")))
fwrite(attack_location_tpms, snakemake@output[["attack_location_tpms"]])

# write log
sessionInfo()
