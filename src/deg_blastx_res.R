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
library(dplyr)

###########
# GLOBALS #
###########

asw_deg_dir <- snakemake@params[["asw_deg_dir"]]
deg_file_suffix <- snakemake@params[["deg_file_suffix"]]
degs_blastx <- snakemake@input[["degs_blastx"]]

########
# MAIN #
########

degs_blastx <- fread(degs_blastx)
setnames(degs_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "blast_hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "nr_blastx_annotation"))
# remove uncharacterized/hypothetical hits
fil_blast <- data.table(dplyr::filter(degs_blastx, !grepl('uncharacterized|GSCOCG0000|hypothetical|unnamed protein', nr_blastx_annotation, ignore.case=TRUE)))
# keep best blast hit
setorder(fil_blast, transcript_id, evalue, -bit_score)
blastx_res <- fil_blast[,.SD[which.min(evalue)], by=transcript_id]
# hits for those filtered out by hypo etc
fil_out <- subset(degs_blastx, !(transcript_id %in% blastx_res$transcript_id))

# for loop to generate DEG tables to save
deg_files <- list.files(asw_deg_dir, deg_file_suffix, full.names = T, recursive=T)
for (x in deg_files) {
  deg_table <- fread(x)
  deg_table_blastx <- merge(deg_table, blastx_res, by="transcript_id", all.x=T)
  myfile <- file.path(gsub(deg_file_suffix, "sig_annots_blastx.csv", x))
  fwrite(deg_table_blastx, file=myfile)
}

# write log
sessionInfo()