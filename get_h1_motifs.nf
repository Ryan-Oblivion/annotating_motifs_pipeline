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

// now getting the file that contains the up peaks and up genes so I can find the gata1 motifs that are bound to promoters or enhancer regions of upregulated genes
params.up_peaks_up_genes = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/multi-results/ATACxRNA/newRNA/rep-up-peaks-50kb-upgenes.bed'
up_peaks_up_genes_ch = Channel.fromPath(params.up_peaks_up_genes)


include {
    
    homer_find_motifs;
    annotate_peaks;
    find_motif_in_promoter;
    annotate_motif_promoter;
    bedtools_intersect;
    annotate_motif_promoter as annotate_motif_promoter_2

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




    // finding the gata1 motifs in bed regions that correspond to up peaks and up genes

    // lets tokenize the file name and get metadata
    up_peaks_up_genes_ch
        .map {file ->

        basename = file.baseName
        filename = file.name
        tokens = basename.tokenize("-")
        tuple("${tokens[1]}_${tokens[2]}_${tokens[4]}", basename, filename, file)


        }
        .set{ uppeaks_upgenes_meta_ch}

    find_motif_in_promoter(uppeaks_upgenes_meta_ch, ref_genome_ch, motif_query_tf_ch)

    annotate_motif_promoter(uppeaks_upgenes_meta_ch, ref_genome_ch, motif_query_tf_ch)


    // okay according to the paper 'https://www.science.org/doi/10.1126/science.aag2445', the two enhancers tested, e-GATA1 or e-HDAC6, are responcible for regulating GATA1.
    // so I need to get the bed regions where these two enhancers are, then do an intersection to find which regions are in the bed files we have as data.
    // save that new intersection bed file that represents the regions in our data that has the enhancers we are looking for
    // then annotate it to see if there are any motifs in the regions

    // I went to UCSC genome browser and found the location of the two genes GATA1 and HDAC6 and downloaded the bed regions of those positions.
    // if i find which genomic locations in our uppeaks-upgenes overlapp with those two regions, i can then look upstream, using homer tool, to find if there are any gata1 motifs found in the enhancer regions of the two genes
    // I put the two files from ucsc in the bin directory and labeled them 'GATA1_genomic_location.bed', 'HDAC6_genomci_location.bed'


    // getting the target gene bed files

    params.target_gene_positions = file('/lustre/fs4/home/rjohnson/pipelines/h1_motif_analysis/bin/*genomic_location.bed')

    target_gene_beds_ch = Channel.value(params.target_gene_positions)

    // make a proces that takes the location bed files of each gene and intersect it with uppeaks up genes bedfile.

    //target_gene_beds_ch.collect().toList().view()
    bedtools_intersect(target_gene_beds_ch.collect(), uppeaks_upgenes_meta_ch) 

    // i used about 10 genes that were considered enhancers in the two papers hera sent. gata1 and hdac6 didnt give any intersections with our data
    // many of the target genes didnt give many intersections either

    intersect_bed_file_ch = bedtools_intersect.out.intersect_target_gene_bed

    intersect_bed_file_ch
        .map {file ->

        basename = file.baseName
        filename = file.name
        tokens = basename.tokenize("-,_")
        tuple("${tokens[1]}_${tokens[2]}_${tokens[4]}_${tokens[5]}", basename, filename, file)


        }
        .set{intersect_bed_metadata}

    intersect_bed_metadata.view()
    annotate_motif_promoter_2(intersect_bed_metadata, ref_genome_ch, motif_query_tf_ch)


}


