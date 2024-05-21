process RUN_COLOC {
    
    tag "$id"
    label 'rProcess'
    publishDir = [
        [
            path: { "${params.outDir}/" },
            mode: 'copy',
            pattern: "graphics/*.png"
        ],
        [
            path: { "${params.outDir}/" },
            mode: 'copy',
            pattern: "tables/*.tsv"
        ]
    ]

    input:
    tuple val(id), path(snplist), path(top_snp), path(ld), path(res_olink)
    each path(outcome_sumstat_files)
    each path(outcome_data_dictionary)
    each path(biomart_gene_annotation)

    output:
    path "graphics/*.png", optional: true
    path "tables/*.tsv", optional: true
    script:
    """
    run_coloc.R \
        $id \
        $snplist \
        $top_snp \
        $ld \
        $res_olink \
        $outcome_sumstat_files \
        $outcome_data_dictionary \
        ${task.cpus} \
        ${params.personal_r_library} \
        ${params.genome_build} \
        $biomart_gene_annotation \
        ${params.recombination_rate_map_dir}
    """
}
