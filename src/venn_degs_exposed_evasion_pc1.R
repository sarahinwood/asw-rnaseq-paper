library(data.table)
library(VennDiagram)
library(dplyr)

evasion_pc1_degs <- fread("output/evasion_deseq2/sig_annots.csv")
exposed_pc1_degs <- fread("output/exposed_deseq2/sig_annots.csv")

vd <- venn.diagram(x = list("Exposed"=exposed_pc1_degs$rn, "Evasion"=evasion_pc1_degs$rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

exposed_evasion <- merge(evasion_pc1_degs, exposed_pc1_degs, by="#gene_id")
fwrite(exposed_evasion, "output/deseq2/exposed_evasion_shared_degs.csv")

sex_degs <- dplyr::filter(exposed_evasion, grepl('sex|sperm|testis|ovary|vitellogenin', sprot_Top_BLASTX_hit.x, ignore.case = TRUE))
fwrite(sex_degs, "output/deseq2/shared_degs_sex.csv")
