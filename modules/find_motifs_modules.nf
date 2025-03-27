
process homer_find_motifs {
    
    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/homer_v5.1_rj'

    publishDir "./up_peaks_results/gata1_specific", mode: 'copy', pattern: '*'

    debug true

    input:
    
    tuple val(expr_design), val(type_choosen), val(basename), val(filename), path(filepath)

    path(ref_genome)

    path(tf_motif_target)

    output:
    path("${homer_motifs_output_dir}"), emit: motif_directory
    
    path("${gata1_outfile}"), emit: gata1_peaks_tsv

    script:

    homer_motifs_output_dir = "./homer_motifs_${expr_design}_${type_choosen}"

    gata1_outfile = "gata1_in_peaks_outfile.tsv"

    """
    #!/usr/bin/env bash

    #echo "\$(less \${tf_motif_target})"

    findMotifsGenome.pl ${filename} ${ref_genome} ${homer_motifs_output_dir} \
    -size 200 \
    -find ${tf_motif_target} \
    -p 20 \
    > ${gata1_outfile}



    """
}

process annotate_peaks {

    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/homer_v5.1_rj'

    publishDir './annotated_peaks', mode:'copy', pattern:'*'


    input:

    tuple val(expr_design), val(type_choosen), val(basename), val(filename), path(filepath)

    path(ref_genome)

    path(tf_motif_target)



    output:

    path("${motif_annotate_outfile}"), emit: motifs_annotated
    path("${bed_motif_outname}")


    script:

    motif_annotate_outfile = "${expr_design}_motifs_annotated_${type_choosen}.tsv"

    bed_motif_outname = "gata1_${expr_design}_${type_choosen}.bed"

    """
    #!/usr/bin/env bash

    

    annotatePeaks.pl ${filename} ${ref_genome} \
    -m ${tf_motif_target} \
    -size 200 \
    -mbed ${bed_motif_outname} \
    > ${motif_annotate_outfile}
    

    """
}


// this will be very similar to the homer_find_motifs process but with different metadata

process find_motif_in_promoter {

    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/homer_v5.1_rj'

    publishDir "./promoter_motifs_found/${type_chosen}/gata1_specific", mode: 'copy', pattern: '*'


    input:

    tuple val(type_chosen), val(basename), val(filename), path(bed_file)

    path(genome)

    path(motif_file)



    output:

    path("${tsv_out_gata1_promoters}"), emit: promoter_gata1_motifs_tsv

    path("${motif_output_dir}"), emit: motif_out_dir

    script:

    motif_output_dir = "./promoter_motifs_${type_chosen}"

    tsv_out_gata1_promoters = "gata1_in_promoter_${type_chosen}.tsv"

    """
    #!/usr/bin/env bash

    findMotifsGenome.pl ${filename} ${genome} ${motif_output_dir} \
    -size -3000,1000 \
    -find ${motif_file} \
    -p 20 \
    > ${tsv_out_gata1_promoters}


    """
}


// now using annotate peaks just like above but to get the motifs in the promoter regions of up peaks up genes

process annotate_motif_promoter {

    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/homer_v5.1_rj'

    publishDir "./annotated_promoters/${type_chosen}", mode:'copy', pattern:'*'


    input:
    tuple val(type_chosen), val(basename), val(filename), path(bed_file)

    path(genome)

    path(motif_file)


    output:

    path("${bed_outname_promoter_motifs}"), emit: promoter_motifs_bed

    path("${out_gata1_promoters}"), emit: annotated_promoter_motif_tsv


    script:

    bed_outname_promoter_motifs = "promoter_motifs_${type_chosen}.bed"

    out_gata1_promoters = "annotated_gata1_promoter_${type_chosen}.tsv"


    """
    #!/usr/bin/env bash

    annotatePeaks.pl ${filename} ${genome} \
    -size -3000,1000 \
    -m ${motif_file} \
    -mbed ${bed_outname_promoter_motifs} \
    > ${out_gata1_promoters}




    """
}

process annotate_motif_enhancers {

    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/homer_v5.1_rj'

    publishDir "./annotated_enhancers/${type_chosen}", mode:'copy', pattern:'*'


    input:
    tuple val(type_chosen), val(basename), val(filename), path(bed_file)

    path(genome)

    path(motif_file)


    output:

    path("${bed_outname_enhancers_motifs}"), emit: enhancer_motifs_bed

    path("${out_gata1_enhancers}"), emit: annotated_enhancer_motif_tsv


    script:

    bed_outname_enhancers_motifs = "enhancers_motifs_${type_chosen}.bed"

    out_gata1_enhancers = "annotated_gata1_enhancers_${type_chosen}.tsv"


    """
    #!/usr/bin/env bash

    annotatePeaks.pl ${filename} ${genome} \
    -size -3000,1000 \
    -m ${motif_file} \
    -mbed ${bed_outname_enhancers_motifs} \
    > ${out_gata1_enhancers}




    """
}


process bedtools_intersect {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    
    label 'normal_small_resources'

    publishDir "./intersection_bed_output", mode: 'copy', pattern: '*'

    //debug true

    input:
    
    path(target_genes_bed)
    //path(target_genes_bed)

    tuple val(type_chosen), val(basename), val(filename), path(bed_file)




    output:

    path("${intersect_file_name}"), emit: intersect_target_gene_bed


    script:
    
    // name_file_1 = "${target_genes_bed[0].baseName}"
    // name_file_2 = "${target_genes_bed[1].baseName}"
    // name_file_3 = "${target_genes_bed[2].baseName}"
    // name_file_4 = "${target_genes_bed[3].baseName}"
    // name_file_5 = "${target_genes_bed[4].baseName}"
    // name_file_6 = "${target_genes_bed[5].baseName}"
    // name_file_7 = "${target_genes_bed[6].baseName}"
    // name_file_8 = "${target_genes_bed[7].baseName}"
    // name_file_9 = "${target_genes_bed[8].baseName}"
    // name_file_10 = "${target_genes_bed[9].baseName}"


    intersect_file_name = "${target_genes_bed.baseName}_intersect.bed"

    """
    #!/usr/bin/env bash

    echo "these are the target gene regions as a list no comma: ${target_genes_bed}"


    ######### parameters  ########

    #-a : bed file that is showing coordinates of interest in experiment.
    #-b : enter the databases of coordinates wanted as over lap in a
    #-wa : write the original entry in A for each overlap
    #-wb: write the original entry in B for each overlap
    #-sortout : when using multiple databases (-b), sort the output DB hits for each record. 
    #-names : will let me put the name i want that represents the name where the alignment came from.
    #-filenames : will report the file name automatically
    ##############################
    
    

    #bedtools intersect -a \${filename} \
    -b \${target_genes_bed} \
    -wb \
    > \${intersect_file_name}


    # By default, if an overlap is found, bedtools intersect reports the shared interval between the two overlapping features.
    # so the resulting file when the genes are (a) and the peaks are (b) will show each interval where the gene overlapped with any parts of any peak 
    # this is why we will get different coordinates for the same gene

    bedtools intersect -a ${target_genes_bed} \
    -b ${filename} \
    > ${intersect_file_name}






    """
}


process kenttools_get_regions {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/uscs_utils_rj'

    publishDir "./enhancer_gene_coordinates", mode:'copy', pattern: '*'

    label 'normal_small_resources'


    input:

    path(gene_txt_file_path)


    output:

    path("${out_file_gene_coordinates}"), emit: enhancer_gene_coor_bed



    script:

    out_file_gene_coordinates = "enhancer_gene_coordinates.bed"

    """
    #!/usr/bin/env bash

    wget -O refGene.txt.gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
    gunzip refGene.txt.gz
    cut -f2-11 refGene.txt > refGene.genePred

    genePredToBed refGene.genePred refGene.bed



    awk -F"\t" '{print \$13"\t"\$2}' refGene.txt > gene_symbol_to_refseq.txt
    #cut -f13,2 refGene.txt > gene_symbol_to_refseq.txt


    grep -F -f ${gene_txt_file_path} gene_symbol_to_refseq.txt | cut -f2 > refseq_ids.txt

    

    grep -F -f refseq_ids.txt refGene.bed > "${out_file_gene_coordinates}"


    # now need to change the gene ids back to gene names tomorow

    """
}

process get_gene_ids_regions {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/uscs_utils_rj'

    publishDir "./h1_up_unchanging_gene_coordinates", mode:'copy', pattern: '*'

    label 'normal_small_resources'


    input:
    
    path(gene_ids_txt)



    output:

    path("${output_bed_for_genes}"), emit: h1_genes_coordinates


    script:

    basename = "${gene_ids_txt.baseName}"

    output_bed_for_genes = "${basename}_coordinates.bed"


    """
    #!/usr/bin/env bash

    wget -O Homo_sapiens.GRCh38.gtf.gz "ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.*.gtf.gz"
    #gunzip Homo_sapiens.GRCh38.gtf.gz

    zcat Homo_sapiens.GRCh38.gtf.gz | awk '\$3=="gene" {
    match(\$0, /gene_id "([^"]+)"/, gene_id);
    chr = (\$1 ~ /^chr/) ? \$1 : "chr"\$1;
    print chr "\t" (\$4 - 1) "\t" \$5 "\t" gene_id[1] "\t" \$7
    }' > gene_bed.bed

    # Filter BED file to keep only the gene IDs in the input file
    grep -F -f ${gene_ids_txt} gene_bed.bed > tmp && mv tmp "${output_bed_for_genes}"









    """
}


process plot_motifs_per_gene {


    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/r_language'

    input:

    path(up_genes_motif_tsv)

    path(unchanging_genes_motif_tsv)


    output:


    script:


    """
    #!/usr/bin/env Rscript

    library(dplyr)

    up_genes = read.csv('gata1_in_promoter_up_h1_genes.tsv', sep='\t',header=TRUE)
    unchanging_genes = read.csv('gata1_in_promoter_unchanging_h1_genes.tsv',sep='\t',header=TRUE)
    
    
    
    up_genes\$GeneID <- sapply(strsplit(as.character(up_genes\$PositionID), "-"), `[`, 1)

    unchanging_genes\$GeneID <- sapply(strsplit(as.character(unchanging_genes\$PositionID), "-"), `[`, 1)

    
    
    
    unchanging_gene_counts <- unchanging_genes %>%
    group_by(GeneID) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count))

    up_gene_counts <- up_genes %>%
    group_by(GeneID) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count))

    # the genes in up genes are not found in the unchanging genes
    
    # just checking something
    sum(up_gene_counts\$Count)
    #[1] 5604
    sum(unchanging_gene_counts\$Count)
    #[1] 71766










    """
}