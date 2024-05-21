process CALCULATE_LD {
    
    cache true
    tag "$olink_id"
    
    input:
    tuple path(snplist),
          path(snpz),
          path(bgen_filtered), 
          path(bgen_filtered_index),
          path(master_file)
    tuple val(olink_id),
          val(chr),
          val(pos_start),
          val(pos_end),
          path(olink_file)

    output:
    tuple val(olink_id),
          val(chr),
          val(pos_start),
          val(pos_end),
          path("ld.*.ld"),
          emit: ld

    script:
    """
    $params.ldstore_binary_path \
    --in-files ${master_file} \
    --write-text \
    --n-threads $params.ldstore_threads \
    --read-only-bgen # \
    # 2>&1 /dev/null > /dev/null # redirect all output to /dev/null
    """
}
