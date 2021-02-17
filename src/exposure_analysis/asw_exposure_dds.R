library("tximport")
library("data.table")
library("DESeq2")

asw_gene_trans_map <- "data/asw_mh_transcriptome/output/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map"
gene2tx <- fread(asw_gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
all_quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##filter out file paths for exposed samples
ex_samples <- grep(pattern='output/asw_mh_concat_salmon/.*_NC|Ex_.*_quant/quant.sf', all_quant_files, value=TRUE)

##assign names to quant files from folder name
names(ex_samples) <- gsub(".*/(.+)_quant/.*", "\\1", ex_samples)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
txi <- tximport(ex_samples, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, "output/exposed_deseq2/exposed_dds.rds")

