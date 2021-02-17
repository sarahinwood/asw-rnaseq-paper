library(data.table)
library(fgsea)
library(ggplot2)

trinotate_report <- fread("/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
res_group <- fread("output/exposed_deseq2/res_group.csv")

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
###negative stat at top of table??
setorder(res_group, stat)
ranks <- res_group[!is.na(stat), stat]
names(ranks) <- res_group[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sum(sorted_fgsea_res$padj<0.05)
##subset into only sig terms and merge w/annotations
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/exposed_deseq2/sig_annot_fgsea_pfam.csv")

