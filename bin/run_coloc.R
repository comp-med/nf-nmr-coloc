#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# These are the arguments passed to the process
nmr_region            <- args[1]
snplist               <- args[2]
top_snp               <- args[3]
ld_file               <- args[4]
res_finemapping_file  <- args[5]
outcome_sumstat_files <- args[6]
outcome_dd_file       <- args[7]

# one core for the master process, the rest for the `mclapply`
n_cores <- as.numeric(args[8])
n_cores <- ifelse(n_cores == 1, n_cores, n_cores - 1)
r_lib   <- args[9]
genome_build <- args[10]

# Files used for plotting
biomart_gene_annotation <- args[11]

# Directory containing recombination rate files for each chromosome
recombination_rate_map_dir <- args[12]

dir.create("./tables/", showWarnings = FALSE)
dir.create("./graphics/", showWarnings = FALSE)

# Libraries needed

suppressPackageStartupMessages(library(data.table, lib.loc = r_lib))
suppressPackageStartupMessages(library(susieR, lib.loc = r_lib))
suppressPackageStartupMessages(library(coloc, lib.loc = r_lib))
suppressPackageStartupMessages(library(doMC, lib.loc = r_lib))
suppressPackageStartupMessages(library(Rfast, lib.loc = r_lib))
suppressPackageStartupMessages(library(glue, lib.loc = r_lib))

# ADDITIONAL SCRIPTS ----

source("../../../modules/templates/plot_locus_compare.R")
source("../../../modules/templates/naive_coloc.R")
source("../../../modules/templates/susie_coloc.R")

# READ DATA ----

outcome_sumstat_files <- fread(
  outcome_sumstat_files,
  sep = ",",
  na.strings = "",
  header = FALSE
)
outcome_sumstat_files[, V1 := gsub(paste0("_", tolower(genome_build)), "", V1)]

# scallop summary statistics
res_finemapping <- fread(res_finemapping_file, na.strings = "")

## import phenotypes to be tested
outcome_dd <- fread(outcome_dd_file, na.strings = "")

nmr_region_parts <- unlist(strsplit(nmr_region, split = "\\."))
# olink <- nmr_region_parts[1] # TODO: this is now several
nmr_phenotypes <- unique(res_finemapping$pheno)
chr.s <- res_finemapping$chrom[1]
pos.s <- res_finemapping$startpos_region[1]
pos.e <- res_finemapping$endpos_region[1]

ld <- fread(
  ld_file,
  data.table = FALSE
)

top.snp <- fread(
  top_snp,
  na.strings = ""
) 

# get recombination rate
recomb_map_chr <- ifelse(chr.s != 23, chr.s, "X")

# TODO: Do this with variable genome build!
rec_file <- paste0(
  recombination_rate_map_dir,
  "/genetic_map_GRCh37_chr",
  recomb_map_chr,
  ".txt"
)

# RUN COLOC ----

# run for each exposure and outcome!
test_plan <- expand.grid(
  exposure = nmr_phenotypes,
  outcome = outcome_sumstat_files$V1, 
  stringsAsFactors = FALSE
)
setDT(test_plan)
test_plan[, test_index := .I]

## run through: outcome_dd$processed_data_name
res_coloc <- mclapply(
  test_plan$test_index,
  function(i) {
    
    # i <-  1
    x <- test_plan[test_index == (i), outcome]
    y <- test_plan[test_index == (i), exposure]
      
    message(glue("[LOG] Outcome: `{x}`"))
    message(glue("[LOG] Exposure: `{y}`"))
    
    # get finemapping info for the single phenotype, only
    res_finemapping_pheno <- copy(res_finemapping[pheno == (y)])
    top_snp_pheno <- copy(top.snp[pheno == (y)])
    phenotype_file <- outcome_sumstat_files[V1 == (x), V2]
    
    # create a snplist for each separate phenotype to get the proxys for
    # Needs to contain `SNP` to grep the first line!
    joint_snplist <- fread(snplist)
    fwrite(
      list(c("SNP", joint_snplist[pheno == (y), id_proxy])), 
      "single_phenotype_snplist.tsv",
      sep = "\t",
      col.names = FALSE
    )
    
    outcome_sumstats_call <- glue("zgrep -wF -f single_phenotype_snplist.tsv {phenotype_file}")
    
    ## get the relevant associations (A2 is the effect allele)
    outcome_stat <- fread(cmd = outcome_sumstats_call)
    
    ## proceed only if at least suggestive evidence
    enough_evidence_to_continue <- nrow(outcome_stat) > 0 & 
      min(outcome_stat$P, na.rm = TRUE) < 1e-4
    
    # Return empty list if FALSE, otherwise, continue
    if (isFALSE(enough_evidence_to_continue)) {
      message("[LOG] Not enough evidence to continue colocalization for this exposure/outcome combination!")
      return(list())
    }
    
    outcome_stats <- glue(
      "zcat {phenotype_file} | ",
      "awk -v chr={chr.s} -v low={pos.s} -v upp={pos.e} ",
      "'{{if(($2 == chr && $3 >= low && $3 <= upp) || NR == 1) print $0}}'"
    )
    
    outcome_stats <- fread(cmd = outcome_stats, sep = "\t")
    
    # Sanity check in case the columns are re-arranged.
    stopifnot(
      "[ERROR] No Outcome Summary Statistics found!" = nrow(outcome_stats) > 0
    )
    
    ## make unique for some reason
    outcome_stats <- unique(outcome_stats)
    
    ## create marker_name to enable mapping to the LD file
    outcome_stats[, marker_name := paste0(
      CHR, ":", BP, ":", pmin(A1, A2), "_", pmax(A1, A2)
    )]
    res_finemapping_pheno[, marker_name := paste0(
      chrom, ":", genpos, ":", pmin(allele0, allele1), "_", pmax(allele0, allele1)
    )]
    
    #-----------------------------#
    ##-- combined set of stats --##
    #-----------------------------#
    
    # check SNP that are not matching
    n_mismatch <- res_finemapping_pheno[
      !(marker_name %in% outcome_stats$marker_name) &
        !(id %in% outcome_stats$SNP), .N
    ]
    message(glue("[LOG] Found {n_mismatch} SNP with rsID mismatches."))
    
    if (is.null(outcome_stats$N)) {
      outcome_stats$N <- outcome_dd[phenotype_id == x, sample_size]
      message(glue("[LOG] Using single sample size value for `runsusie()`"))
    }
    
    ## combine
    res_all <- merge(
      res_finemapping_pheno,
      outcome_stats[, c("marker_name", "A1", "A2", "BETA", "SE", "P", "N")],
      by = "marker_name",
      suffixes = c("_exposure", "_outcome")
    )
    
    # Coloc depends on whether the trais is `quantitative`/`binary`
    outcome_type <- outcome_dd[phenotype_id == (x), trait_type]
    
    # Binary traits need the number of cases, annotate them according to 
    # availability
    if (outcome_type == "binary" & is.null(outcome_stats$N_CAS)) {
      
      res_all$outcome_n_cases <- outcome_dd[phenotype_id == x, n_cases]
      
    } else if (outcome_type == "binary" & !is.null(outcome_stats$N_CAS)) {
      
      res_all <- merge(
        res_all,
        outcome_stats[, .(marker_name, outcome_n_cases = N_CAS)]
      )
      
    }
    
    if (nrow(res_all) == 0) {
      message("[LOG] No SNP overlap between exposure and outcome data!")
      message(glue("[LOG] Skipping Colocalization with `{x}`"))
      return(list(NA))
    }
    
    ## align effect estimates (this currently ignores INDELs --> check!)
    res_all[, Effect_outcome := ifelse(toupper(allele1) == A2, BETA, -BETA)]
    res_all[, StdErr_outcome := SE]
    
    
    # create a regular pvalue for the exposure (is log10)
    res_all[, pval_marginal := 10^(-pval_marginal)]
    res_all[, pval_joint := 10^(-pval_joint)]
    
    #-----------------------------#
    ## --      naive coloc     --##
    #-----------------------------#
    
    ## run coloc
    res_naive <- naive.coloc(
      res_all = res_all,
      ld = ld,
      x = x,
      y = y,
      x_is_binary = ifelse(outcome_type == "binary", TRUE, FALSE),
      nmr_region = nmr_region,
      res_finemapping = res_finemapping_pheno,
      biomart_gene_annotation = biomart_gene_annotation,
      rec_file = rec_file
    )
    ## add type of coloc
    res_naive[, type := "naive"]
    
    #-----------------------------#
    ## --       susie coloc     --##
    #-----------------------------#
    
    ## import function to do so
    
    ## run coloc
    res_susie <- susie.coloc(
      res.all = res_all, 
      res.olink = res_finemapping_pheno,
      ld = ld, 
      top.snp = top_snp_pheno,
      x = x,
      y = y,
      nmr_region = nmr_region,
      x_is_binary = ifelse(outcome_type == "binary", TRUE, FALSE),
      biomart_gene_annotation = biomart_gene_annotation,
      rec_file = rec_file
    )
    
    ## add type of coloc
    if (!is.null(res_susie) > 0) {
      res_susie[, type := "susie"]
    }
    
    #-----------------------------#
    ## --    create return     --##
    #-----------------------------#
    
    ## combine
    res_coloc <- rbindlist(
      list(res_naive, res_susie),
      use.names = TRUE,
      fill = TRUE
    )
    
    ## return results
    return(res_coloc)
    
  }, mc.cores = n_cores)

## combine
res_coloc <- rbindlist(res_coloc, fill = TRUE, use.names = TRUE)

if(nrow(res_coloc) > 0) {
  ## store results
  write.table(
    res_coloc,
    paste0(paste("./tables/coloc_results_region", nmr_region, sep = "_"), ".tsv"),
    sep = "\t",
    row.names = FALSE
  )
}

