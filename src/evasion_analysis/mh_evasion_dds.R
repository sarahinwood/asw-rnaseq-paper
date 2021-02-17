library("tximport")
library("data.table")
library("DESeq2")

mh_gene_trans_map <- "data/asw_mh_transcriptome/output/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map"
gene2tx <- fread(mh_gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
all_quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##filter out file paths for exposed samples - should NOT match pattern therefore invert=TRUE
ev_samples <- grep(pattern='output/asw_mh_concat_salmon/.*_NC|Ex_.*_quant/quant.sf', all_quant_files, value=TRUE, invert=TRUE)

##assign names to quant files from folder name
names(ev_samples) <- gsub(".*/(.+)_quant/.*", "\\1", ev_samples)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
txi <- tximport(ev_samples, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, "output/evasion_deseq2/mh_evasion_dds.rds")

dds$group <- factor(paste(dds$parasitism,sep="_"))
##plot counts for PCR target gene
plotCounts(dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))
mh_pcr_counts <- plotCounts(dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"), returnData = TRUE)
fwrite(data.table(mh_pcr_counts, keep.rownames=TRUE), "output/evasion_deseq2/mh_pcr_counts.csv")
