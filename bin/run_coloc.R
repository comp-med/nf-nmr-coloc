#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# These are the arguments passed to the process
olink_id              <- args[1]
snplist               <- args[2]
top_snp               <- args[3]
ld_file               <- args[4]
res_olink_file        <- args[5]
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
res_olink <- fread(res_olink_file, na.strings = "")

## import phenotypes to be tested
outcome_dd <- fread(outcome_dd_file, na.strings = "")

olink_id_parts <- unlist(strsplit(olink_id, split = "\\."))
olink <- olink_id_parts[1]
chr.s <- olink_id_parts[2]
pos.s <- olink_id_parts[3]
pos.e <- olink_id_parts[4]

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

## run through: outcome_dd$processed_data_name
res_coloc <- mclapply(
  outcome_sumstat_files$V1,
  # x <- outcome_sumstat_files$V1[42]
  # x <- outcome_sumstat_files$V1[43]
  function(x) {
    
    message(glue("[LOG] Phenotype: `{x}`"))
    
    phenotype_file <- outcome_sumstat_files[V1 == (x), V2]
    outcome_sumstats_call <- glue("zgrep -wF -f {snplist} {phenotype_file}")
    
    ## get the relevant associations (A2 is the effect allele)
    outcome_stat <- fread(cmd = outcome_sumstats_call)
    
    ## proceed only if at least suggestive evidence
    if (nrow(outcome_stat) > 0 & min(outcome_stat$P, na.rm = TRUE) < 1e-4) {
      
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
      
      ## create MarkerName to enable mapping to the LD file
      outcome_stats[, MarkerName := paste0(
        CHR, ":", BP, ":", pmin(A1, A2), "_", pmax(A1, A2)
      )]
      
      #-----------------------------#
      ##-- combined set of stats --##
      #-----------------------------#
      
      # check SNP that are not matching
      n_mismatch <- res_olink[
        (MarkerName %in% outcome_stats$MarkerName) &
          !(id %in% outcome_stats$SNP), .N
      ]
      message(glue("[LOG] Found {n_mismatch} SNP with rsID mismatches."))
      
      if (is.null(outcome_stats$N)) {
        outcome_stats$N <- outcome_dd[phenotype_id == x, sample_size]
        message(glue("[LOG] Using single sample size value for `runsusie()`"))
      }
      
      ## combine
      res_all <- merge(
        res_olink,
        outcome_stats[, c("MarkerName", "A1", "A2", "BETA", "SE", "P", "N")],
        by = "MarkerName",
        suffixes = c(".pQTL", ".outcome")
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
          outcome_stats[, .(MarkerName, outcome_n_cases = N_CAS)]
        )
        
      }
      
      if (nrow(res_all) == 0) {
        message("[LOG] No SNP overlap between Olink and CVD data!")
        message(glue("[LOG] Skipping Colocalization with `{x}`"))
        return(list(NA))
      }
      
      ## align effect estimates (this currently ignores INDELs --> check!)
      res_all[, Effect.outcome := ifelse(toupper(EA) == A2, BETA, -BETA)]
      res_all[, StdErr.outcome := SE]
      
      #-----------------------------#
      ## --      naive coloc     --##
      #-----------------------------#
      
      ## run coloc
      res_naive <- naive.coloc(
        res.all = res_all,
        ld = ld,
        x = x,
        x_is_binary = ifelse(outcome_type == "binary", TRUE, FALSE),
        olink_id = olink_id,
        res.olink = res_olink,
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
        res.olink = res_olink,
        ld = ld, 
        top.snp = top.snp,
        x = x, 
        x_is_binary = ifelse(outcome_type == "binary", TRUE, FALSE),
        olink_id = olink_id,
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
    }
  }, mc.cores = n_cores)

## combine
res_coloc <- rbindlist(res_coloc, fill = T)

## store results
write.table(
  res_coloc,
  paste("./tables/coloc.results", olink_id, "tsv", sep = "."),
  sep = "\t",
  row.names = FALSE
)
