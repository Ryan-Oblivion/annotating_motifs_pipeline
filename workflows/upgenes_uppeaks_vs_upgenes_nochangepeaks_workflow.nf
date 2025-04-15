
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



workflow upgenes_uppeaks_vs_upgenes_no_changepeaks_workflow {


    take:

    upgenes_uppeaks
    upgenes_nochange_peaks

    ref_genome_ch

    motif_query_tf_ch 




    main:

    

    
}