process TOP_SNP_AND_PROXIES {

    cache true
    tag "$nmr_finemap_region"
    label 'rProcess'

    input:
    tuple val(nmr_finemap_region),
          path(finemap_results_region), 
          path(ld_matrix),
          path(snplist)

    output:
    tuple val(nmr_finemap_region),
          path("proxy_snplist.tsv"),
          path("top_snp.tsv"),
          path(ld_matrix),
          path("res_finemapping_with_proxies.tsv")

    script:
    """
    top_snp_and_proxies.R \
        ${nmr_finemap_region} \
        ${finemap_results_region} \
        ${ld_matrix} \
        ${snplist} \
        ${params.personal_r_library}
    """
}
