process CREATE_INPUT_FILES {
    
    tag "$nmr_finemap_region"
    label 'rProcess'

    input: 
    each path(nmr_finemap_master_file)
    tuple val(nmr_finemap_region), path(nmr_finemap_region_dir)
    each path(nmr_finemapping_results_directory)

    output:
    tuple val(nmr_finemap_region), 
          path("finemap_results_region.tsv"), 
          path("ld_file"), 
          path("snplist_file")

    script:
    """
    create_input_files.R \
        $nmr_finemap_master_file \
        $nmr_finemap_region \
        $nmr_finemap_region_dir \
        $nmr_finemapping_results_directory \
        $params.personal_r_library
    """
}
