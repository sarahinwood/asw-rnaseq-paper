#!/usr/bin/env python3
import peppy

#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['r1', 'r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

all_species = ["ASW", "Mh", "MhFV"]
all_locations = ["Dunedin", "Ruakura"]

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
salmon_container = 'docker://combinelab/salmon:1.5.1'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.14:0.0.1' # DESeq2 v1.34
bioconductor_cowplot_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'
multiqc_container = 'docker://ewels/multiqc:v1.11'
blast_container= 'docker://ncbi/blast:2.13.0'

#########
# RULES #
#########

rule target:
    input:
            ## QC
        'output/multiqc/multiqc_report.html',
        'output/triple_deseq2/read_mapping/table.csv',
        expand('output/triple_deseq2/{species}/PCA/parasitism_location_1.pdf', species=["ASW", "Mh", "MhFV"]),
        'output/triple_deseq2/ASW/PCA/top500_genes_annots.csv',
        'output/triple_deseq2/ASW/location/power_res.csv',
        'output/tpms/location_tpms.csv',
        'output/tpms/exposure_tpms.csv',
        'output/tpms/exposure_location_tpms.csv',
        'output/tpms/attack_tpms.csv',
        'output/tpms/attack_location_tpms.csv',
        'output/tpms/parasitism_tpms.csv',
        'output/tpms/parasitism_location_tpms.csv',
            ## DGE - location
        'output/triple_deseq2/ASW/location/location_sig_annots.csv',
           ## DGE - exposure
        'output/triple_deseq2/ASW/exposure_samples/power_res.csv',
        'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots.csv',
        'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots.csv',
            ## DGE - microcosm
        expand('output/triple_deseq2/{species}/microcosm_samples/{species}_all_dds.rds', species=["ASW"]),
        'output/triple_deseq2/ASW/microcosm_samples/power_res.csv',
                    # attack
        expand('output/triple_deseq2/{species}/microcosm_samples/all_microcosm_attacked/attack_sig_annots.csv', species=["ASW"]),
        expand('output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/sig_annots.csv', species=["ASW"]),
                    # parasitism
        expand('output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/para_sig_annots.csv', species=["ASW", "Mh", "MhFV"]),
        expand('output/triple_deseq2/{species}/microcosm_samples/microcosm_parasitism_{location}/para_sig_annots.csv', species=["ASW", "Mh", "MhFV"], location=["Dunedin"]),
                    # blastx degs
        'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_sig_annots_blastx.csv',
                    # overrep analysis
        expand('output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_BlastP_GO_overrep.csv', analysis=["para"]), # "exposure", "location" # "attack", "att_loc_int", "exposure-loc-int" - no sig GO
        expand('output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_Pfam_overrep.csv', analysis=["para"]), # "exposure", "attack" #  "exposure-loc-int", "att_loc_int", "location"  - no sig Pfam #
        expand('output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_full_plot.pdf',
            analysis=["para"]), #"exposure", # "location", "exposure-loc-int", "att_loc_int",  -  no DEGs with Pfam/GO,
                                            # "Dunedin_para", "attack" - no GO enrichment
        expand('output/triple_deseq2/ASW/FGSEA/{analysis}_FGSEA_enrichment.pdf', analysis=["para"]) #"location", "exposure", "exposure-loc-int", "att_loc_int", "attack"

#############################
## term overrepresentation ##
#############################

rule plot_overrepresentation:
    input:
        pfam_file = 'output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_Pfam_overrep.csv',
        go_file  = 'output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_BlastP_GO_overrep.csv'
    output:
        enrich_plot = 'output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_overrep_plot.pdf',
        full_plot = 'output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_full_plot.pdf'
    singularity:
        bioconductor_cowplot_container
    log:
        'output/logs/term_overrepresentation/plot_{analysis}_overrepresentation.log'
    script:
        'src/term_overrepresentation/plot_overrepresentation.R'

rule term_overrepresentation:
    input:
        DEG_file = 'output/triple_deseq2/ASW/linked_deg_files/{analysis}_sig_annots.csv',
        term_annot_table_file = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_annots.csv',
        term_to_gene_file = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_to_gene.csv',
        term_to_name_file = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_to_name.csv'
    output:
        enrichment_table = 'output/triple_deseq2/ASW/term_overrepresentation/{analysis}/{analysis}_{term}_overrep.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/term_overrepresentation/{analysis}_analysis_{term}_enrichment.log'
    script:
        'src/term_overrepresentation/{wildcards.term}_overrepresentation.R'

rule term_to_gene:
    input:
        trinotate_file = 'data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt'
    output:
        term_annot_table = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_annots.csv',
        term_to_gene = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_to_gene.csv',
        term_to_name = 'output/triple_deseq2/ASW/term_overrepresentation/{term}_annots/{term}_to_name.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/term_overrepresentation/{term}_term_to_gene.log'
    script:
        'src/term_overrepresentation/{wildcards.term}_to_geneID.R'

rule link_deg_files:
    input:
        degs1 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots.csv',
        degs2 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots.csv',
        degs3 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_attacked/attack_sig_annots.csv',
        degs4 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_att_loc_int/sig_annots.csv',
        degs5 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_sig_annots.csv',
        degs6 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/microcosm_parasitism_Dunedin/para_sig_annots.csv',
        degs7 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/location/location_sig_annots.csv'
    output:
        degs1 = 'output/triple_deseq2/ASW/linked_deg_files/exposure_sig_annots.csv',
        degs2 = 'output/triple_deseq2/ASW/linked_deg_files/exposure-loc-int_sig_annots.csv',
        degs3 = 'output/triple_deseq2/ASW/linked_deg_files/attack_sig_annots.csv',
        degs4 = 'output/triple_deseq2/ASW/linked_deg_files/att_loc_int_sig_annots.csv',
        degs5 = 'output/triple_deseq2/ASW/linked_deg_files/para_sig_annots.csv',
        degs6 = 'output/triple_deseq2/ASW/linked_deg_files/Dunedin_para_sig_annots.csv',
        degs7 = 'output/triple_deseq2/ASW/linked_deg_files/location_sig_annots.csv'
    shell:
        'ln -s {input.degs1} {output.degs1} & '
        'ln -s {input.degs2} {output.degs2} & '
        'ln -s {input.degs3} {output.degs3} & '
        'ln -s {input.degs4} {output.degs4} & '
        'ln -s {input.degs5} {output.degs5} & '
        'ln -s {input.degs6} {output.degs6} & '
        'ln -s {input.degs7} {output.degs7} & '

######################
## FGSEA enrichment ##
######################

rule FGSEA_analysis:
    input:
        trinotate_file = 'data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_annotation_report.txt',
        res_group_file = 'output/triple_deseq2/ASW/linked_deg_files/{analysis}_res_group.csv'
    output:
        FGSEA_res = 'output/triple_deseq2/ASW/FGSEA/{analysis}_FGSEA_enrichment.csv',
        FGSEA_plot = 'output/triple_deseq2/ASW/FGSEA/{analysis}_FGSEA_enrichment.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/FGSEA_{analysis}_analysis_enrichment.log'
    script:
        'src/FGSEA_enrichment.R'

rule link_res_group_files:
    input:
        degs1 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/exposure_samples/exposure/exposure_res_group.csv',
        degs2 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_res_group.csv',
        degs3 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_attacked/attack_res_group.csv',
        degs4 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_att_loc_int/res_group.csv',
        degs5 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_res_group.csv',
        degs6 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/microcosm_samples/microcosm_parasitism_Dunedin/para_res_group.csv',
        degs7 = '/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-rnaseq-pc1/output/triple_deseq2/ASW/location/location_res_group.csv'
    output:
        degs1 = 'output/triple_deseq2/ASW/linked_deg_files/exposure_res_group.csv',
        degs2 = 'output/triple_deseq2/ASW/linked_deg_files/exposure-loc-int_res_group.csv',
        degs3 = 'output/triple_deseq2/ASW/linked_deg_files/attack_res_group.csv',
        degs4 = 'output/triple_deseq2/ASW/linked_deg_files/att_loc_int_res_group.csv',
        degs5 = 'output/triple_deseq2/ASW/linked_deg_files/para_res_group.csv',
        degs6 = 'output/triple_deseq2/ASW/linked_deg_files/Dunedin_para_res_group.csv',
        degs7 = 'output/triple_deseq2/ASW/linked_deg_files/location_res_group.csv'
    shell:
        'ln -s {input.degs1} {output.degs1} & '
        'ln -s {input.degs2} {output.degs2} & '
        'ln -s {input.degs3} {output.degs3} & '
        'ln -s {input.degs4} {output.degs4} & '
        'ln -s {input.degs5} {output.degs5} & '
        'ln -s {input.degs6} {output.degs6} & '
        'ln -s {input.degs7} {output.degs7} & '

################
## blast DEGs ##
################

rule deg_blastx_res:
    input:
        degs_blastx = 'output/triple_deseq2/ASW/degs_BlastX/blastx.outfmt6',
        degs1 = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots.csv',
        degs2 = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots.csv',
        degs3 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_attacked/attack_sig_annots.csv',
        degs4 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_att_loc_int/sig_annots.csv',
        degs5 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_sig_annots.csv',
        degs6 = 'output/triple_deseq2/ASW/microcosm_samples/microcosm_parasitism_Dunedin/para_sig_annots.csv'
    output:
        degs1 = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots_blastx.csv',
        degs2 = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots_blastx.csv',
        degs3 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_attacked/attack_sig_annots_blastx.csv',
        degs4 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_att_loc_int/sig_annots_blastx.csv',
        degs5 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_sig_annots_blastx.csv',
        degs6 = 'output/triple_deseq2/ASW/microcosm_samples/microcosm_parasitism_Dunedin/para_sig_annots_blastx.csv'
    params:
        asw_deg_dir = 'output/triple_deseq2/ASW',
        deg_file_suffix = 'sig_annots.csv'
    log:
        'output/logs/deg_blastx_res.log'
    script:
        'src/deg_blastx_res.R'

rule blastx_degs:
    input:
        degs_to_blast = 'output/triple_deseq2/ASW/degs_BlastX/degs_to_blast.fasta'
    output:
        blastx_res = 'output/triple_deseq2/ASW/degs_BlastX/blastx.outfmt6'
    params:
        blastdb = 'bin/blast_db/nr_2023_03_20/nr'
    threads:
        40
    log:
        'output/logs/blastx_degs.log'
    singularity:
        blast_container
    shell:
        'blastx '
        '-query {input.degs_to_blast} '
        '-db {params.blastdb} '
        '-num_threads {threads} '
        '-evalue 1e-5 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '

rule collate_degs_to_blast:
    input:
        degs1 = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots.csv',
        degs2 = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots.csv',
        degs3 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_attacked/attack_sig_annots.csv',
        degs4 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_att_loc_int/sig_annots.csv',
        degs5 = 'output/triple_deseq2/ASW/microcosm_samples/all_microcosm_parasitism/para_sig_annots.csv',
        degs6 = 'output/triple_deseq2/ASW/microcosm_samples/microcosm_parasitism_Dunedin/para_sig_annots.csv',
        fasta = 'data/asw_mh_MhFV.fasta'
    output:
        degs_to_blast = 'output/triple_deseq2/ASW/degs_BlastX/degs_to_blast.fasta'
    params:
        asw_deg_dir = 'output/triple_deseq2/ASW',
        deg_file_suffix = 'sig_annots.csv'
    log:
        'output/logs/collate_degs_to_blast.log'
    script:
        'src/collate_degs_to_blast.R'

###############
## microcosm ##
###############

## parasitism analyses - Ruakura para ASW likely para for much longer and so using this instead of interaction?
rule parasitism_locsp_pairwise: # error on Ruakura MhFV as not enough to cluster - why aren't more genes DE?
    input:
        dds_file ='output/triple_deseq2/{species}/microcosm_samples/{location}_dds.rds',
        annotation_file = 'data/annotations/{species}_annots.csv',
        parasitism_location_tpms = 'output/tpms/parasitism_location_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/{species}/microcosm_samples/microcosm_parasitism_{location}/para_dds.rds',
        res_table = 'output/triple_deseq2/{species}/microcosm_samples/microcosm_parasitism_{location}/para_res_group.csv',
        sig_annots = 'output/triple_deseq2/{species}/microcosm_samples/microcosm_parasitism_{location}/para_sig_annots.csv',
        col_clustered_heatmap = 'output/triple_deseq2/{species}/microcosm_samples/microcosm_parasitism_{location}/col_clust_heatmap.pdf'
    log:
        'output/logs/deseq2_parasitism_{species}_{location}.log'
    singularity:
        bioconductor_container
    script:
        'src/microcosm_analysis/parasitism_loc-sp_pairwise.R' 

rule parasitism_pairwise:
    input:
        dds_file ='output/triple_deseq2/{species}/microcosm_samples/{species}_all_dds.rds',
        annotation_file = 'data/annotations/{species}_annots.csv',
        parasitism_tpms_file = 'output/tpms/parasitism_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/para_dds.rds',
        res_table = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/para_res_group.csv',
        sig_annots = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/para_sig_annots.csv',
        row_clustered_heatmap = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/row_clust_heatmap.pdf',
        col_clustered_heatmap = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_parasitism/col_clust_heatmap.pdf'
    log:
        'output/logs/deseq2_parasitism_{species}_all.log'
    singularity:
        bioconductor_container
    script:
        'src/microcosm_analysis/parasitism_pairwise.R'

## attacked analyses - could do with only unpara too!!! <----------
rule attacked_location_interaction:
    input:
        dds_file ='output/triple_deseq2/{species}/microcosm_samples/{species}_all_dds.rds',
        annotation_file = 'data/annotations/{species}_annots.csv',
        attack_location_tpms = 'output/tpms/attack_location_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/dds.rds',
        res_table = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/res_group.csv',
        sig_annots = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/sig_annots.csv',
        #row_clustered_heatmap = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/row_clust_heatmap.pdf',
        #col_clustered_heatmap = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_att_loc_int/col_clust_heatmap.pdf'
    log:
        'output/logs/deseq2_attacked_location_interaction_{species}_all.log'
    singularity:
        bioconductor_container
    script:
        'src/microcosm_analysis/attacked_location_interaction.R'

rule attacked_pairwise: # no DEGs so no heatmap
    input:
        dds_file ='output/triple_deseq2/{species}/microcosm_samples/{species}_all_dds.rds',
        annotation_file = 'data/annotations/{species}_annots.csv',
        attack_tpms = 'output/tpms/attack_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_attacked/attack_dds.rds',
        res_table = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_attacked/attack_res_group.csv',
        sig_annots = 'output/triple_deseq2/{species}/microcosm_samples/all_microcosm_attacked/attack_sig_annots.csv'
    log:
        'output/logs/deseq2_attacked_{species}_all.log'
    singularity:
        bioconductor_container
    script:
        'src/microcosm_analysis/attacked_pairwise.R'

rule microcosm_power_analysis:
    input:
        dds_file = 'output/triple_deseq2/ASW/microcosm_samples/ASW_all_dds.rds'
    output:
        power_res = 'output/triple_deseq2/ASW/microcosm_samples/power_res.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/microcosm_power_analysis.log'
    script:
        'src/microcosm_analysis/power_analysis.R'

rule microcosm_dds:
    input:
        dds_file = 'output/triple_deseq2/{species}/{species}_dds.rds'
    output:
        microcosm_dds = 'output/triple_deseq2/{species}/microcosm_samples/{species}_all_dds.rds',
        dunedin_microcosm_dds = 'output/triple_deseq2/{species}/microcosm_samples/Dunedin_dds.rds',
        ruakura_microcosm_dds = 'output/triple_deseq2/{species}/microcosm_samples/Ruakura_dds.rds'
    singularity:
        bioconductor_container
    log:
        'output/logs/microcosm_dds_{species}.log'
    script:
        'src/microcosm_analysis/microcosm_dds.R'

##############
## exposure ##
##############

rule exposure_location_interaction:
    input:
        dds_file = 'output/triple_deseq2/ASW/exposure_samples/ASW_dds.rds',
        annotation_file = 'data/annotations/ASW_annots.csv',
        exposure_location_tpms = 'output/tpms/exposure_location_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_dds.rds',
        res_table = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_res_group.csv',
        sig_annots = 'output/triple_deseq2/ASW/exposure_samples/exposure-loc-int/int_sig_annots.csv'
    log:
        'output/logs/deseq2_exposure_location_interaction.log'
    singularity:
        bioconductor_container
    script:
        'src/exposure_analysis/exposure_location_interaction.R'

rule exposure_pairwise:
    input:
        dds_file = 'output/triple_deseq2/ASW/exposure_samples/ASW_dds.rds',
        annotation_file = 'data/annotations/ASW_annots.csv',
        exposure_tpms = 'output/tpms/exposure_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_dds.rds',
        res_table = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_res_group.csv',
        sig_annots = 'output/triple_deseq2/ASW/exposure_samples/exposure/exposure_sig_annots.csv'
    log:
        'output/logs/deseq2_exposure.log'
    singularity:
        bioconductor_container
    script:
        'src/exposure_analysis/exposure_pairwise.R'

rule exposure_power_analysis:
    input:
        dds_file = 'output/triple_deseq2/ASW/exposure_samples/ASW_dds.rds'
    output:
        power_res = 'output/triple_deseq2/ASW/exposure_samples/power_res.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/exposure_power_analysis.log'
    script:
        'src/exposure_analysis/power_analysis.R'

rule exposure_dds:
    input:
        dds_file = 'output/triple_deseq2/ASW/ASW_dds.rds'
    output:
        exposure_dds = 'output/triple_deseq2/ASW/exposure_samples/ASW_dds.rds'
    singularity:
        bioconductor_container
    log:
        'output/logs/exposure_dds.log'
    script:
        'src/exposure_analysis/asw_exposure_dds.R'

##################################
## comparisons with all samples ##
##################################

rule location_pairwise:
    input:
        dds_file = 'output/triple_deseq2/ASW/ASW_dds.rds',
        annotation_file = 'data/annotations/ASW_annots.csv',
        location_tpms = 'output/tpms/location_tpms.csv'
    output:
        dds_file = 'output/triple_deseq2/ASW/location/location_dds.rds',
        res_table = 'output/triple_deseq2/ASW/location/location_res_group.csv',
        sig_annots = 'output/triple_deseq2/ASW/location/location_sig_annots.csv',
        row_clustered_heatmap = 'output/triple_deseq2/ASW/location/row_clust_heatmap.pdf',
        col_clustered_heatmap = 'output/triple_deseq2/ASW/location/col_clust_heatmap.pdf'
    log:
        'output/logs/deseq2_location_ASW.log'
    singularity:
        bioconductor_container
    script:
        'src/all_analysis/all_location_pairwise.R'

rule location_power_analysis:
    input:
        dds_file = 'output/triple_deseq2/ASW/ASW_dds.rds'
    output:
        power_res = 'output/triple_deseq2/ASW/location/power_res.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/location_power_analysis.log'
    script:
        'src/all_analysis/all_location_power_analysis.R'

###################
## make dds & QC ##
###################

## PC1 is sex - control for in all analyses
rule ASW_PC1_genes:
    input:
        genes_file = 'output/triple_deseq2/ASW/PCA/top500_genes.csv',
        trinotate_file = 'data/asw_mh_transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv'
    output:
        genes_annots = 'output/triple_deseq2/ASW/PCA/top500_genes_annots.csv'
    log:
        'output/logs/ASW_PC1_genes.log'
    singularity:
        bioconductor_container
    script:
        'src/QC/ASW_PC1_genes.R'

rule PCA:
    input:
        dds_file = 'output/triple_deseq2/{species}/{species}_dds.rds'
    output:
        PCA_paraloc_1 = 'output/triple_deseq2/{species}/PCA/parasitism_location_1.pdf',
        PCA_paraloc_2 = 'output/triple_deseq2/{species}/PCA/parasitism_location_2.pdf',
        PCA_paraloc_3 = 'output/triple_deseq2/{species}/PCA/parasitism_location_3.pdf',
        PCA_attloc = 'output/triple_deseq2/{species}/PCA/attack_location.pdf',
        PCA_concentration = 'output/triple_deseq2/{species}/PCA/concentration.pdf',
        PCA_RQN = 'output/triple_deseq2/{species}/PCA/RQN.pdf',
        PCA_weightings = 'output/triple_deseq2/{species}/PCA/gene_weightings.csv',
        PCA_top500 = 'output/triple_deseq2/{species}/PCA/top500_genes.csv'
    log:
        'output/logs/PCA_{species}.log'
    singularity:
        bioconductor_container
    script:
        'src/QC/PCA.R'

rule read_mapping_stats:
    input:
        asw_dds_file = 'output/triple_deseq2/ASW/ASW_dds.rds',
        mh_dds_file = 'output/triple_deseq2/Mh/Mh_dds.rds',
        MhFV_dds_file = 'output/triple_deseq2/MhFV/MhFV_dds.rds'
    output:
        read_mapping_table = 'output/triple_deseq2/read_mapping/table.csv',
        read_mapping_summary = 'output/triple_deseq2/read_mapping/summary.csv',
        all_stacked = 'output/triple_deseq2/read_mapping/stacked_bar.pdf',
        asw_box = 'output/triple_deseq2/read_mapping/asw_box.pdf',
        mh_box = 'output/triple_deseq2/read_mapping/mh_box.pdf',
        MhFV_box = 'output/triple_deseq2/read_mapping/MhFV_box.pdf',
        Mh_PCR_target = 'output/triple_deseq2/read_mapping/Mh_PCR_target.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/read_mapping_stats.log'
    script:
        'src/QC/read_mapping.R'

rule tpm_table:
    input:
        tpm_file = 'output/tpms/sample_TPMs.csv'
    output:
        location_tpms = 'output/tpms/location_tpms.csv',
        exposure_tpms = 'output/tpms/exposure_tpms.csv',
        exposure_location_tpms = 'output/tpms/exposure_location_tpms.csv',
        attack_tpms = 'output/tpms/attack_tpms.csv',
        attack_location_tpms = 'output/tpms/attack_location_tpms.csv',
        parasitism_tpms = 'output/tpms/parasitism_tpms.csv',
        parasitism_location_tpms = 'output/tpms/parasitism_location_tpms.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/tpm_table.log'
    script:
        'src/tpm_table.R'

rule make_dds:
    input:
        gene2tx_file = 'data/asw_mh_transcriptome/output/asw_mh_MhFV/asw_mh_MhFV_tx2gene',
        sample_data_file = 'data/sample_table.csv',
        quant_files = expand('output/asw_mh_MhFV_concat_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        salmon_tpm = 'output/tpms/sample_TPMs.csv',
        #salmon_counts = 'output/tpms/sample_counts.csv',
        mh_asw_MhFV_dds = 'output/triple_deseq2/all_dds.rds',
        asw_dds = 'output/triple_deseq2/ASW/ASW_dds.rds',
        mh_dds = 'output/triple_deseq2/Mh/Mh_dds.rds',
        MhFV_dds = 'output/triple_deseq2/MhFV/MhFV_dds.rds'
    log:
        'output/logs/make_dds.log'
    singularity:
        bioconductor_container
    script:
        'src/make_dds.R'

#################################
## map to asw/mh/MhFV combined ##
#################################

rule multiqc:
    input:
        salmon = expand('output/asw_mh_MhFV_concat_salmon/{sample}_quant/quant.sf', sample=all_samples),
        fastqc = 'output/fastqc'
    output:
        'output/multiqc/multiqc_report.html'
    params:
        outdir = 'output/multiqc',
        indirs = ['output/asw_mh_MhFV_concat_salmon', 'output/fastqc', 'output/logs/bbduk_trim']
    log:
        'output/logs/multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc -f ' ##force to write over old output
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

rule concat_salmon_quant:
    input:
        index_output = 'output/asw_mh_MhFV_concat_salmon/transcripts_index/refseq.bin',
        left = 'output/bbduk_trim/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/asw_mh_MhFV_concat_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/asw_mh_MhFV_concat_salmon/transcripts_index',
        outdir = 'output/asw_mh_MhFV_concat_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon/asw_mh_concat_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '--writeUnmappedNames '
        '-p {threads} '
        '&> {log}'

rule concat_salmon_index:
    input:
        fasta = 'data/ASW_Mh_MhFV_updated.fasta'
    output:
        'output/asw_mh_MhFV_concat_salmon/transcripts_index/refseq.bin'
    params:
        outdir = 'output/asw_mh_MhFV_concat_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_mh_concat_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.fasta} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

######################
## prep for mapping ##
######################

rule fastqc:
    input:
        expand('output/bbduk_trim/{sample}_r{n}.fq.gz',
            sample=all_samples, n=[1,2])
    output:
        directory('output/fastqc')
    shell:
        'mkdir -p {output} ; '
        'fastqc --outdir {output} {input}'

rule bbduk_trim:
    input:
        unpack(get_reads)
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

##OG file with sequencing for this project in different structure - two folders
##folder with small files is a small miseq run
##folder with large files are standard hiseq files to be used for analysis