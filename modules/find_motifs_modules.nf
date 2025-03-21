
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