#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# Get arguments from nextflow process

# TODO remove!
# nmr_finemap_master_file <- "/sc-projects/sc-proj-computational-medicine/people/Martijn/02_UKBB_NMR/03_GWAS_NMR/02_fine_mapping/input/unique_finemapped_regions.txt"
# nmr_finemap_region <- "10_1"
# nmr_finemap_region_dir <- "/sc-projects/sc-proj-computational-medicine/people/Martijn/02_UKBB_NMR/03_GWAS_NMR/02_fine_mapping/output/region_data/10_1/"
# nmr_finemapping_results_directory <- "/sc-projects/sc-proj-computational-medicine/people/Martijn/02_UKBB_NMR/03_GWAS_NMR/02_fine_mapping/output/finemapping/"

nmr_finemap_master_file <- args[1]
nmr_finemap_region <- args[2]
nmr_finemap_region_dir <- args[3]
nmr_finemapping_results_directory <- args[4]
r.lib <- args[5]

# Required packages
library("data.table", lib.loc = r.lib)
library("glue", lib.loc = r.lib)

# Main Script ----

# Get all finemapping result directories for the input region
master_file <- fread(nmr_finemap_master_file)
master_file <- master_file[region_name == nmr_finemap_region]

finemap_results_region <- lapply(master_file$index, function(x) {
  fread(glue("{nmr_finemapping_results_directory}/{x}/finemapped_results.txt"))
})
finemap_results_region <- rbindlist(finemap_results_region)

# TODO Change col names?

finemap_results_region[, chr := ifelse(chr == "X", 23, chr)]

## create identifier column in results to keep the mapping to the LD matrix
finemap_results_region[, snp.id := 1:nrow(finemap_results_region)]

## order, keep in mind that cs -1 means not included in the credible set
finemap_results_region <- finemap_results_region[order(-cs, -pip)]

## create indicator to select
finemap_results_region[, ind := 1:.N, by = "cs"]

fwrite(
  finemap_results_region,
  "finemap_results_region.tsv",
  sep = "\t"
)

# export the paths of these
ld_file <- glue("{nmr_finemap_region_dir}/ldmat.ld")
snplist_file <- glue("{nmr_finemap_region_dir}/snplist.txt")



#### TODO output
# snplist
# res.olink
# ld
