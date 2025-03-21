
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

    publishDir "./promoter_motifs_upgenes_uppeaks/gata1_specific", mode: 'copy', pattern: '*'


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

    publishDir './annotated_promoters', mode:'copy', pattern:'*'


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