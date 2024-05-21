process LD_INPUT_FILES {
    
    cache true
    tag "$olink_id"
    label 'rProcess'

    input: 
    tuple val(olink_id), 
          val(chr),
          val(pos_start), 
          val(pos_end), 
          path(olink_file)
    each path(ukb_bgen_directory)
    each path(ukb_sample_inclusion_list)

    output:
    tuple path("snplist.*.lst"),
          path("snpz.*.z"),
          path("filtered.*.bgen"), 
          path("filtered.*.bgen.bgi"),
          path("master.*.z"),
          emit: files
    tuple val(olink_id),
          val(chr),
          val(pos_start),
          val(pos_end),
          path("res_olink.tsv"),
          emit: ids

    script:
    """
    ld_input_files.R \
        $olink_id \
        $chr \
        $pos_start \
        $pos_end \
        $olink_file \
        $ukb_bgen_directory \
        $ukb_sample_inclusion_list \
        $params.personal_r_library \
        $params.bgenix_binary_path
    """
}
