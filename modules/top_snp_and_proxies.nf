process TOP_SNP_AND_PROXIES {

    cache true
    tag "$id"
    label 'rProcess'

    input:
    tuple val(id), path(ld), path(res_olink)

    output:
    tuple val(id), path("${id}.snplist"), path("top_snp.tsv"), path(ld), path("res_olink_with_proxies.tsv")

    script:
    """
    top_snp_and_proxies.R ${id} ${ld} ${res_olink} ${params.personal_r_library}
    """
}
