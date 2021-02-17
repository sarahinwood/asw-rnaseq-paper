library("tximport")
library("data.table")
library("DESeq2")

mh_gene_trans_map <- "data/asw_mh_transcriptome/output/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map"
gene2tx <- fread(mh_gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
all_quant_files <- list.files(path="output/asw_mh_concat_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)

##assign names to quant files from folder name
names(all_quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", all_quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(all_quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, "output/deseq2/mh_all_dds.rds")

total_reads_per_sample <- colSums(counts(dds))
table_counts <- data.table(data.frame(total_reads_per_sample), keep.rownames=TRUE)
fwrite(table_counts, "output/deseq2/mh_total_counts.csv")


