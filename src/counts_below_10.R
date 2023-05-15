library(data.table)

counts <- fread("output/tpms/sample_counts.csv")
counts$rowsums <- rowSums(counts[,c(2:104)])

counts_rowsums <- counts[,c(1, 105)]
counts_rowsums_below_10 <- subset(counts_rowsums, rowsums<10)
fwrite(counts_rowsums_below_10, "output/tpms/genes_counts_below_10.csv")

counts_rowsums_above_10 <- subset(counts_rowsums, rowsums>10)
fwrite(counts_rowsums_below_10, "output/tpms/genes_counts_above_10.csv")

# keeps only 54% of genes
(length(counts_rowsums_above_10$rn)/length(counts_rowsums$rn))*100
# but 99.99% of reads
(sum(counts_rowsums_above_10$rowsums)/sum(counts_rowsums$rowsums))*100
