#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# Get arguments from nextflow process
nmr_finemap_master_file <- args[1]
nmr_finemap_region <- args[2]
nmr_finemap_region_dir <- args[3]
nmr_finemapping_results_directory <- args[4]
r.lib <- args[5]

# Required packages
library("data.table", lib.loc = r.lib)
library("glue", lib.loc = r.lib)

# Main Script ----

# Do not export the LD matrix itself, just the link
ld_file <- glue("{Sys.readlink(nmr_finemap_region_dir)}/ldmat.ld")
snplist_file <- glue("{Sys.readlink(nmr_finemap_region_dir)}/snplist.txt")

# This contains the order of the LD matrix!
snplist <- fread(snplist_file, header = FALSE)
setnames(snplist, "snp")

# Get all finemapping result directories for the input region
master_file <- fread(nmr_finemap_master_file)
master_file <- master_file[region_name == nmr_finemap_region]

finemap_results_region <- lapply(master_file$index, function(x) {
  fread(glue("{nmr_finemapping_results_directory}/{x}/finemapped_results.txt"))
})
finemap_results_region <- rbindlist(finemap_results_region)
finemap_results_region[, chrom := ifelse(chrom == "X", 23, chrom)]

## create identifier column in results to keep the mapping to the LD matrix
stopifnot("SNP shoul be in snplist!" = all(finemap_results_region$id %in% snplist$snp))
m1 <- match(finemap_results_region$id, snplist$snp)
nrow(finemap_results_region)
finemap_results_region[, snp_id := m1]

## order, keep in mind that cs -1 means not included in the credible set
finemap_results_region <- finemap_results_region[order(-cs, -pip), .SD, by = pheno]

## create indicator to select
finemap_results_region[, ind := 1:.N, by = .(cs, pheno)]

# Output ----

fwrite(
  finemap_results_region,
  "finemap_results_region.tsv",
  sep = "\t"
)

# write the location of the files to a file to extract later
cat(ld_file, file = "ld_file",sep = "")
cat(snplist_file, file = "snplist_file",sep = "")