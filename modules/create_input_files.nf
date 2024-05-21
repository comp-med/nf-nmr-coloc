process CREATE_INPUT_FILES {
    
    cache true
    tag "$olink_id"
    label 'rProcess'

    input: 
    each path(nmr_finemap_master_file)
    tuple val(nmr_finemap_region), path(nmr_finemap_region_dir)
    each path(nmr_finemapping_results_directory)

    output:
    val nmr_finemap_region
    

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
