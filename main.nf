#!/usr/bin/env nextflow

// https://github.com/AustralianBioCommons/Nextflow_DSL2_template
nextflow.enable.dsl=2

// HEADER ---------------------------------------------------------------------

log.info """\
===============================================================================
Colocalization Analysis Pipeline based on UKB NMR GWAS results
===============================================================================

Created by the Computational Medicine Group | BIH @ Charité

===============================================================================
Workflow run parameters 
===============================================================================
input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}
===============================================================================

"""

// Help function
def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf 

  Required Arguments:

  <TODO>

  Optional Arguments:

  --outDir	Specify path to output directory. Default is `output/`
	
""".stripIndent()
}

// MODULES --------------------------------------------------------------------

include {

    CREATE_INPUT_FILES;

} from './modules/create_input_files.nf'

include {

    TOP_SNP_AND_PROXIES

} from './modules/top_snp_and_proxies.nf'

// include {
// 
//     BIOMART_GENE_ANNOTATION;
// 
// } from './modules/biomart_gene_annotation.nf'
// 
// include {
// 
//     RUN_COLOC
// 
// } from './modules/run_coloc.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {

    nmr_finemapping_master_file_ch = Channel.fromPath(
      "$params.nmr_finemap_master_table",
      checkIfExists: true,
      type: 'file'      
    )

    // Channel for each region
    nmr_credible_set_regions_ch = Channel.fromPath(
      "$params.nmr_ld_directory/*",
      checkIfExists: true,
      type: 'dir'
    ).map({ path -> ["${path.getName()}", "${path}"] })
    
    // Channel for the finemapping results directory
    nmr_finemapping_results_directory_ch = Channel.fromPath(
      "$params.nmr_credible_sets_directory/",
      checkIfExists: true,
      type: 'dir'
    )
    
    // For each region, get all credible set data for each phenotype 
    // from master file
    coloc_input_files = CREATE_INPUT_FILES (
      nmr_finemapping_master_file_ch,
      nmr_credible_set_regions_ch,
      nmr_finemapping_results_directory_ch 
    )

    // Extract the content of the ld matrix and snplist files to values
    coloc_input_files = coloc_input_files
 	.map({
	    region, region_file, ld_file, snplist_file -> 
		[region, region_file, ld_file.text, snplist_file.text]
	 })

    // Use the joint input data
    coloc_input_files_with_top_snps = TOP_SNP_AND_PROXIES(
         coloc_input_files
    )

    // For colocalization, outcome data is needed
    // outcome_sumstat_file_ch = Channel.fromPath(
    //     "${params.outcome_sumstat_directory}/**/${params.genome_build}/*.tsv.bgz",
    //     checkIfExists: true
    // )
    // .map({ path -> "${path.getSimpleName()}, ${path}" })
    // .collectFile( 
    //     name: "all_outcome_files.csv",
    //     newLine: true
    // )
// 
    // // Outcome data dictionary with sample sizes needed for colocalization
    // outcome_data_dictionary_ch = Channel.fromPath(
    //     "$params.outcome_data_dictionary",
    //     checkIfExists: true
    // )
// 
    // // Download the gene annotation from ENSEMBL using biomaRt
    // biomart = BIOMART_GENE_ANNOTATION()
// 
    // // This is the workhorse of the pipeline - the colocalization analysis
    // RUN_COLOC(
    //     snplist_ld_coloc_input,
    //     outcome_sumstat_file_ch,
    //     outcome_data_dictionary_ch,
    //     biomart.gene_annotation
    // )
}

// SUMMARY --------------------------------------------------------------------

workflow.onComplete {
summary = """
===============================================================================
Workflow execution summary
===============================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

===============================================================================
  """
println summary
}
