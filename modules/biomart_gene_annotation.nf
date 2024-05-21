// This is needed in the plotting function within the colocalization functions
process BIOMART_GENE_ANNOTATION {
    
    cache true
    label 'rProcess'

    output:
    path("Genes.GRCh37.complete.tsv"), emit: gene_annotation

    script:
    """
    #! /usr/bin/env Rscript

    library(biomaRt)
    library(data.table)
      
    # Need to specify the genome build from the settings
    genome_build <- "$params.genome_build"
    genome_build <- as.numeric(gsub("GRCh", "", genome_build))
    gene.ensembl <- useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        GRCh = genome_build
    )
    tmp.genes <- getBM(
      attributes = c(
        'chromosome_name',
        'start_position',
        'end_position',
        'ensembl_gene_id',
        'external_gene_name',
        'gene_biotype'
      ),
      filters = c('chromosome_name'),
      values = list(c(1:22, "X")),
      mart = gene.ensembl
    )
    tmp.genes <- unique(tmp.genes[, c(
      'chromosome_name', 
      'start_position',
      'end_position',
      'ensembl_gene_id',
      'external_gene_name',
      'gene_biotype'
    )])   
    fwrite(
        tmp.genes,
        "Genes.GRCh37.complete.tsv",
        sep = "\t"
        )
    """
}
