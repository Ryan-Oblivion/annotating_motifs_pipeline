nextflow.enable.dsl=2

// the bed files that will be used for finding motifs in homer
// /lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/

params.scrvslow_up = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/rep-peaks-scrvslow-up.bed'

up_scrvslow_ch = Channel.fromPath(params.scrvslow_up)


params.scrvslow_nochange = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/rep-peaks-scrvslow-nochange.bed'

nochange_scrvslow_ch  = Channel.fromPath(params.scrvslow_nochange)

// now to get the genome
params.ref_genome = '/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa'
ref_genome_ch = Channel.fromPath(params.ref_genome)

// get the gata1 motif file. the user would be able to specify the path to any specific tf that they want to find the motif position for
params.query_motif = '/lustre/fs4/home/rjohnson/pipelines/h1_motif_analysis/bin/gata1.motif'
motif_query_tf_ch = Channel.fromPath(params.query_motif)



include {
    
    homer_find_motifs;
    annotate_peaks

}from './modules/find_motifs_modules.nf'


workflow {

    // just checking if the paths are visible
    //up_scrvslow_ch.view()
    //nochange_scrvslow_ch.view()

    // manipulating the channel to get some base names and tokens for meta data
    up_scrvslow_ch
        .map { file -> 

        basename = file.baseName
        filename = file.name
        tokens = basename.tokenize("-")
        tuple( tokens[2],"${tokens[3]}_${tokens[1]}", basename, filename, file)
        
        }
        .set {up_scrvslow_meta}

    //up_scrvslow_meta.view()

    // using this process to get homer motifs
    homer_find_motifs(up_scrvslow_meta, ref_genome_ch, motif_query_tf_ch)

    // the output of gata1 motifs is better here becasue it annotates the actual peaks with that information
    annotate_peaks(up_scrvslow_meta, ref_genome_ch, motif_query_tf_ch)

}
