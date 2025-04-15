

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

}from '../modules/find_motifs_modules.nf'

workflow analyzing_basemean_genes_workflow {



    take:
    wtvslowup_gene
    wtvslowdown_nochange
    ref_genome_ch
    motif_query_tf_ch




    main:

    //wtvslowup_gene.view()

    // so i need to get these files in the format needed for the process

    wtvslowup_gene
        .splitCsv(header:false, sep:'\t')
        .map{row ->

        // i want to only have the important columns chr name, start, end, geneid
        [row[0],row[1],row[2],row[3]].join('\t')

        }
        .collectFile(name:'wtvslow_up_genebody.bed', newLine:true, storeDir:'./bin')
        .set{wtvslowup_gene_ch}

    // now making a txt file for the no change genes

    wtvslowdown_nochange
        .splitCsv(header:false, sep:'\t')
        .map{row ->

        // i want to only have the important columns chr name, start, end, geneid
        [row[0],row[1],row[2],row[3]].join('\t')

        }
        .collectFile(name:'WTvslow_down-basemeanmatch_nochange-genebody.bed', newLine:true, storeDir:'./bin')
        .set{wtvslowdown_nochange_gene_ch}

    // just copying the steps from the main workflow but using these files above

    // get_gene_ids_regions(wtvslowup_gene_ch)  // this will get the up gene regions
    // get_gene_ids_regions_2(wtvslowdown_nochange_gene_ch)  // this gets the unchanging gene regions

    // // combine the out files from each channel

    // up_h1_coor_file = get_gene_ids_regions.out.h1_genes_coordinates
    // unchanging_h1_coor_file = get_gene_ids_regions_2.out.h1_genes_coordinates

    // up_and_unchanging_coor_files = up_h1_coor_file.combine(unchanging_h1_coor_file).flatten()

    // now looking for the motifs of the upgenes and motifs of unchanging genes


    ///////////////
    // using the file directly from the meta channel because it already has the coordinates and I dont need to do the process that gets coordinates


    wtvslowup_gene_ch
        .map {file ->

        basename = file.baseName
        filename = file.name
        tokens = basename.tokenize("_")
        tuple("${tokens[0]}_${tokens[1]}_${tokens[2]}", basename, filename, file)


        }
        .set{ upgenes_h1_meta_ch} // how the channel looks [up_h1_genes, up_h1_genes_coordinates, up_h1_genes_coordinates.bed, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/h1_motif_analysis/work/65/29edcaa27520db3a4eaac85e689cff/up_h1_genes_coordinates.bed]

    find_motif_in_promoter(upgenes_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)


    // and motifs of unchanging genes

    wtvslowdown_nochange_gene_ch
        .map {file ->

        basename = file.baseName
        filename = file.name
        tokens = basename.tokenize("_")
        tuple("${tokens[0]}_${tokens[1]}_${tokens[2]}", basename, filename, file)


        }
        .set{ unchanging_h1_meta_ch} // [unchanging_h1_genes, unchanging_h1_genes_coordinates, unchanging_h1_genes_coordinates.bed, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/h1_motif_analysis/work/40/0a4dcfac72aebd8cc283899e449b60/unchanging_h1_genes_coordinates.bed]
    
    
    
    find_motif_in_promoter_2(unchanging_h1_meta_ch, ref_genome_ch, motif_query_tf_ch)


    //getting the gene motif tsv files
    up_genes_motif_tsv = find_motif_in_promoter.out.promoter_gata1_motifs_tsv  // using this file gata1_in_promoter_up_h1_genes.tsv
    unchanging_genes_motif_tsv = find_motif_in_promoter_2.out.promoter_gata1_motifs_tsv

    // now make a process that works on comparing how many motifs show up for each gene
    plot_motifs_per_gene(up_genes_motif_tsv, unchanging_genes_motif_tsv)










}


