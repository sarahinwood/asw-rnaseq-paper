library(data.table)

trinotate_file = fread('data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv')

## exposure:location
FGSEA_res_exloc = fread('output/triple_deseq2/ASW/FGSEA/exposure-loc-int_FGSEA_enrichment.csv')

nuc_binding <- subset(FGSEA_res_exloc, pathway=="GO:0003676")
nuc_binding_LE <- tstrsplit(nuc_binding$leadingEdge, "|", fixed=T)
annots_nuc_binding_LE <- subset(trinotate_file, `#gene_id` %in% nuc_binding_LE)
fwrite(annots_nuc_binding_LE, "output/triple_deseq2/ASW/FGSEA/leadingedge/exloc_nuc_binding_annots.csv")

sig_trans <- subset(FGSEA_res_exloc, pathway=="GO:0035556")
sig_trans_LE <- tstrsplit(sig_trans$leadingEdge, "|", fixed=T)
annots_sig_trans_LE <- subset(trinotate_file, `#gene_id` %in% sig_trans_LE)
fwrite(annots_sig_trans_LE, "output/triple_deseq2/ASW/FGSEA/leadingedge/exloc_sig_trans_annots.csv")

## attack:location
FGSEA_res_attloc = fread('output/triple_deseq2/ASW/FGSEA/att_loc_int_FGSEA_enrichment.csv')

reg_trans <- subset(FGSEA_res_attloc, pathway=="GO:0006355")
reg_trans_LE <- tstrsplit(reg_trans$leadingEdge, "|", fixed=T)
annots_reg_trans_LE <- subset(trinotate_file, `#gene_id` %in% reg_trans_LE)
fwrite(annots_reg_trans_LE, "output/triple_deseq2/ASW/FGSEA/leadingedge/attloc_reg_trans_annots.csv")

ss_DNA_bind <- subset(FGSEA_res_attloc, pathway=="GO:0043565")
ss_DNA_bind_LE <- tstrsplit(ss_DNA_bind$leadingEdge, "|", fixed=T)
annots_ss_DNA_bind_LE <- subset(trinotate_file, `#gene_id` %in% ss_DNA_bind_LE)
fwrite(annots_ss_DNA_bind_LE, "output/triple_deseq2/ASW/FGSEA/leadingedge/attloc_ss_DNA_bind_annots.csv")

