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

// include {
// 
//     TOP_SNP_AND_PROXIES
// 
// } from './modules/top_snp_and_proxies.nf'
// 
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
    ).view()

    // Channel for each region
    nmr_credible_set_regions_ch = Channel.fromPath(
      "$params.nmr_ld_directory/*",
      checkIfExists: true,
      type: 'dir'
    ).map({ path -> "${path.getName()}, ${path}" })
     .view()
    
    // Channel for the finemapping results directory
    nmr_credible_set_regions_ch = Channel.fromPath(
      "$params.nmr_credible_sets_directory/",
      checkIfExists: true,
      type: 'dir'
    ).view()
    
    // For each region, get all credible set data for each phenotype 
    // from master file
    //CREATE_INPUT_FILES (
    //  
    //)

    // SCALLOP protein input is a list of files from the specified directory
    // scallop_cis_credible_sets_files = Channel.fromPath(
    //     "${params.scallop_credible_sets_directory}/fine.mapped.*.txt",
    //     checkIfExists: true,
    //     type: 'file'
    //     )
// 
    // // Turn directory with input files into list with Protein ID,
    // // chromosome, start, end and file path tuple
    // parsed_scallop_input = scallop_cis_credible_sets_files 
    //     .map { file -> [
    //         file.getName().tokenize('.')[2],
    //         file.getName().tokenize('.')[3],
    //         file.getName().tokenize('.')[4],
    //         file.getName().tokenize('.')[5],
    //         file] 
    //         }
// 
    // // Create input files necessary to calculate LD between SNPs
    // coloc_input = CREATE_INPUT_FILES(
    //     parsed_scallop_input,
    // )
// 
    // // Create a joint channel by key for each ID from the two processes
    // ld_calculation_input_with_key = ld_calculation_input.ids
    //     .map({
    //             olink_id, chr, pos_start, pos_end, olink_file ->
    //             ["${olink_id}.${chr}.${pos_start}.${pos_end}", olink_file]
    //          })
    // ld_matrix_with_key = ld_matrix.ld
    //     .map(
    //     {
    //         olink_id, chr, pos_start, pos_end, ld ->
    //         ["${olink_id}.${chr}.${pos_start}.${pos_end}", ld]
    //         }
    //     )
    // snp_files_and_ld_matrix = ld_matrix_with_key
    //     .join(ld_calculation_input_with_key)
// 
    // // Use the joint channel as input
    // snplist_ld_coloc_input = TOP_SNP_AND_PROXIES(
    //      snp_files_and_ld_matrix
    // )
// 
    // // For colocalization, outcome data is needed
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
