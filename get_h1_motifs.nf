nextflow.enable.dsl=2

// the bed files that will be used for finding motifs in homer
// /lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/

params.scrvslow_up = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/rep-peaks-scrvslow-up.bed'

up_scrvslow_ch = Channel.fromPath(params.scrvslow_up)


params.scrvslow_nochange = '/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/DEseq/rep-peaks-scrvslow-nochange.bed'

nochange_scrvslow_ch  = Channel.fromPath(params.scrvslow_nochange)

// now to get the genome
params.ref_genome = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa')
ref_genome_ch = Channel.value(params.ref_genome)

// get the gata1 motif file. the user would be able to specify the path to any specific tf that they want to find the motif position for
params.query_motif = file('/lustre/fs4/home/rjohnson/pipelines/h1_motif_analysis/bin/gata1.motif')
motif_query_tf_ch = Channel.value(params.query_motif)

// now getting the file that contains the up peaks and up genes so I can find the gata1 motifs that are bound to promoters or enhancer regions of upregulated genes
params.up_peaks_up_genes = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/multi-results/ATACxRNA/newRNA/rep-up-peaks-50kb-upgenes.bed')
up_peaks_up_genes_ch = Channel.value(params.up_peaks_up_genes)


// getting the supplement 2 list from the nature paper on k562 cell enhancers.
params.supplement_gRNAs =  '/lustre/fs4/home/rjohnson/pipelines/h1_motif_analysis/bin/41467_2024_52490_MOESM5_ESM.csv'
gRNA_info_csv = Channel.fromPath(params.supplement_gRNAs)


////////////////////////////////////////////

// this is the file from hera that contains the Gene_IDs, in the Gene_ID column, for genes that are up (up genes) and not significant (unchanging genes) in the "type" column
params.target_h1_genes_csv = '/lustre/fs4/home/rjohnson/pipelines/h1_motif_analysis/bin/WTvslow_with_gene_names.csv'
h1_target_genes_ch = Channel.fromPath(params.target_h1_genes_csv)


params.wtvs_lowup = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowup-genebody.bed')
wtvslowup_genebody_ch = Channel.value(params.wtvs_lowup)

params.wtvs_lowdown_nochange = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowdown-basemeanmatchnochange-genebody.bed')
wtvslowdown_nochange_ch = Channel.value(params.wtvs_lowdown_nochange)





include {
    
    homer_find_motifs;
    annotate_peaks;
    find_motif_in_promoter;
    annotate_motif_promoter;
    bedtools_intersect;
    annotate_motif_enhancers;
    kenttools_get_regions;
    get_gene_ids_regions;
    get_gene_ids_regions as get_gene_ids_regions_2;
    find_motif_in_promoter as find_motif_in_promoter_2;
    annotate_motif_promoter as annotate_motif_promoter_2;
    plot_motifs_per_gene

}from './modules/find_motifs_modules.nf'

include {

    analyzing_basemean_genes_workflow


}from './workflows/motif_workflow_basemean.nf'

workflow {


    ////// Making a new workflow for organization ////

    analyzing_basemean_genes_workflow(wtvslowup_genebody_ch, wtvslowdown_nochange_ch, ref_genome_ch, motif_query_tf_ch )





    //////////////// starting a new ////////////////

    // this will parse the file and give back only the up gene names
    // h1_target_genes_ch
    //     .splitCsv(header:true)
    //     .filter{ row ->

    //     up_genes = row.type.startsWith("Up")

    //     //unchanging= row.type.startsWith("Not significant")
    //     //tuple(up_genes, unchanging)

    //     }
    //     .map {row ->
            
    //         //up_gene_id = tuple(row.Gene_ID, row.type)
    //         row.Gene_ID

    //     }
    //     .collectFile(name:'up_h1_genes.txt', newLine: true, storeDir:'./bin')
    //     .set{up_h1_genes_file_ch}


    // // this will parse the file and give back the unchanging gene list
    // h1_target_genes_ch
    //     .splitCsv(header:true)
    //     .filter{ row ->

    //     up_genes = row.type.startsWith("Not significant")

    //     //unchanging= row.type.startsWith("Not significant")
    //     //tuple(up_genes, unchanging)

    //     }
    //     .map {row ->
            
    //         //unchanging_gene_id = tuple(row.Gene_ID, row.type)
    //         row.Gene_ID

    //     }
    //     .collectFile(name:'unchanging_h1_genes.txt', newLine: true, storeDir:'./bin')
    //     .set{unchanging_h1_genes_file_ch}
        
        
    // // now that i have the list of genes that are unchanging and the list that were up genes in their respective files, i can make a process to get the coordinates for each gene
    // // make a process to get regions from gene_ids
    // // i could've just combined the two files and flattened them and sent that into one process and nextflow wouldve handeled parallelizing the process
    // get_gene_ids_regions(up_h1_genes_file_ch)  // this will get the up gene regions
    // get_gene_ids_regions_2(unchanging_h1_genes_file_ch)  // this gets the unchanging gene regions 


    // // combine the out files from each channel

    // up_h1_coor_file = get_gene_ids_regions.out.h1_genes_coordinates
    // unchanging_h1_coor_file = get_gene_ids_regions_2.out.h1_genes_coordinates

    // up_and_unchanging_coor_files = up_h1_coor_file.combine(unchanging_h1_coor_file).flatten()
    // //up_and_unchanging_coor_files.view()


    // // now looking for the motifs of the upgenes and motifs of unchanging genes

    // up_h1_coor_file
    //     .map {file ->

    //     basename = file.baseName
    //     filename = file.name
    //     tokens = basename.tokenize("_")
    //     tuple("${tokens[0]}_${tokens[1]}_${tokens[2]}", basename, filename, file)


    //     }
    //     .set{ upgenes_h1_meta_ch} // how the channel looks [up_h1_genes, up_h1_genes_coordinates, up_h1_genes_coordinates.bed, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/h1_motif_analysis/work/65/29edcaa27520db3a4eaac85e689cff/up_h1_genes_coordinates.bed]

    // find_motif_in_promoter(upgenes_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)

    // annotate_motif_promoter(upgenes_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)

    // // and motifs of unchanging genes

    // unchanging_h1_coor_file
    //     .map {file ->

    //     basename = file.baseName
    //     filename = file.name
    //     tokens = basename.tokenize("_")
    //     tuple("${tokens[0]}_${tokens[1]}_${tokens[2]}", basename, filename, file)


    //     }
    //     .set{ unchanging_h1_meta_ch} // [unchanging_h1_genes, unchanging_h1_genes_coordinates, unchanging_h1_genes_coordinates.bed, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/h1_motif_analysis/work/40/0a4dcfac72aebd8cc283899e449b60/unchanging_h1_genes_coordinates.bed]
    
    
    
    // find_motif_in_promoter_2(unchanging_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)

    // annotate_motif_promoter_2(unchanging_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)

    // //getting the gene motif tsv files
    // up_genes_motif_tsv = find_motif_in_promoter.out.promoter_gata1_motifs_tsv  // using this file gata1_in_promoter_up_h1_genes.tsv
    // unchanging_genes_motif_tsv = find_motif_in_promoter_2.out.promoter_gata1_motifs_tsv

    // // now make a process that works on comparing how many motifs show up for each gene
    // plot_motifs_per_gene(up_genes_motif_tsv, unchanging_genes_motif_tsv)

    // // now to find the intersection of the up and unchanging regions with the uppeaks upgenes file
    // // this is the uppeaks-upgenes file up_peaks_up_genes_ch
    // // using this process bedtools_intersect

    // // I already created a meta channel for the uppeaks up genes channel and that is used as input for the bedtools_intersect process

    // up_peaks_up_genes_ch
    //     .map {file ->

    //     basename = file.baseName
    //     filename = file.name
    //     tokens = basename.tokenize("-")
    //     tuple("${tokens[1]}_${tokens[2]}_${tokens[4]}", basename, filename, file)


    //     }
    //     .set{ uppeaks_upgenes_meta_ch}

    

    
    
    
    // not doing intersection yet

    // put this meta channel second in the bedtools_intersect process
    // put the up genes file and unchanging coordinate files first
    //bedtools_intersect(up_and_unchanging_coor_files, uppeaks_upgenes_meta_ch )

        
        
        //.view() 

        /*.multiMap { up, unchanging ->

            up_genes: up
            unchanging_genes: unchanging

        }
        //.set{h1_genes_separated}
        //.view()
    
    //h1_genes_separated.up_genes.view{ it -> "the up genes: $it"}
    //h1_genes_separated.unchanging_genes.view{ it -> "the unchanging genes: $it"}
    */






    ////////////////////////////////////////////////

    // just checking if the paths are visible
    //up_scrvslow_ch.view()
    //nochange_scrvslow_ch.view()

    // manipulating the channel to get some base names and tokens for meta data
    
    /*
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

    //intersect_bed_metadata.view()
    annotate_motif_enhancers(intersect_bed_metadata, ref_genome_ch, motif_query_tf_ch)


    //gRNA_info_csv.view()
    // make a process to grab only the lines that have guides that target enhancers.
    // didnt need to make a process. I just manipulated the file in nextflow
    gRNA_info_csv
        .splitCsv(header:true)
        .filter {row ->

            name_of_guide = row.target_guide
            enhancer_guides = name_of_guide.endsWith("enhancer")

            // now only keep the gene names from the rows that have enhancer as the target guides

            //gene_name = enhancer_guides.target_gene

            
            
            
        }
        .map {row ->

            gene_name = row.target_gene

        }
        .distinct()
        .collectFile(name: 'enhancer_genes.txt', storeDir: './bin', newLine: true )
        .set{enhancer_gene_list_ch}
        //.view{it -> "this is a gene that was targeted from enhancer guide: ${it}"}

    // now making a process that gets the enhancer_genes file and uses kent-tools to query ucsc and get the regions.

    enhancer_gene_list_ch.view()
    kenttools_get_regions(enhancer_gene_list_ch)

*/


}


