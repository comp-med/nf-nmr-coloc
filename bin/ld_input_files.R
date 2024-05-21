#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

# Get arguments from nextflow process
olink              <- args[1]
chr.s              <- as.numeric(args[2])
pos.s              <- as.numeric(args[3])
pos.e              <- as.numeric(args[4])
olink.path         <- args[5]
bgen.directory     <- args[6]
ukb.inclusion.list <- args[7]
r.lib              <- args[8]
bgenix.bin         <- args[9]

# Required packages
library("data.table", lib.loc = r.lib)
library("glue", lib.loc = r.lib)

# Main Script ----

res.olink <- fread(olink.path)
res.olink[, chr := ifelse(chr == "X", 23, chr)]

## write list of SNPs to be queried to file
snplist_output_file <- paste(
    "snplist", olink, chr.s, pos.s, pos.e, "lst",
    sep = "."
)
write.table(
  res.olink$rsid,
  snplist_output_file,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

## snp data to be queried
# Positions are in build37, might differ from sentinel variant
# file with positions in build38
tmp.z <- as.data.frame(res.olink[, c("rsid", "chr", "pos", "NEA", "EA")])
names(tmp.z) <- c("rsid", "chromosome", "position", "allele1", "allele2")

## adopt chromosome if needed
if (chr.s < 10) {
  tmp.z$chromosome <- paste0("0", tmp.z$chromosome)
} else if (chr.s == 23) {
  tmp.z$chromosome <- "X"
}

## write to file
write.table(
  tmp.z,
  paste(
    "snpz",
    olink, chr.s, pos.s,
    pos.e, "z", sep = "."
  ),
  row.names = FALSE,
  quote = FALSE
)

# Filter BGEN Files ----

bgen_output_file <- paste(
  "filtered",
  olink, chr.s, pos.s, pos.e,
  "bgen",
  sep = "."
)
bgen_index_output_file <- paste(
  "filtered",
  olink, chr.s, pos.s, pos.e,
  "bgen.bgi",
  sep = "."
)

if (chr.s == 23) {

  bgen_file <- glue("{bgen.directory}/ukb22828_cX_b0_v3.bgen")
  bgenix_call <- glue(
    "{bgenix.bin} ",
    "-g {bgen_file} ",
    "-incl-rsids {snplist_output_file} > {bgen_output_file}")
  system(bgenix_call)

} else {

  bgen_file <- glue("{bgen.directory}/ukb22828_c{chr.s}_b0_v3.bgen")
  bgenix_call <- glue(
    "{bgenix.bin} ",
    "-g {bgen_file} ",
    "-incl-rsids {snplist_output_file} > {bgen_output_file}")
  system(bgenix_call)

}

bgenix_index_call <- glue("{bgenix.bin} -g {bgen_output_file} -index")
system(bgenix_index_call)

## --> create master file for LDstore2 <-- ##

## assign entries
m.file <- data.frame(
  z = paste("snpz", olink, chr.s, pos.s, pos.e, "z", sep = "."),
  bgen = bgen_output_file,
  bgi = bgen_index_output_file,
  ld = paste("ld", olink, chr.s, pos.s, pos.e, "ld", sep = "."),
  incl = ukb.inclusion.list,
  n_samples = ifelse(chr.s != 23, 487409, 486757)
)

## write to file
write.table(
  m.file,
  paste(
    "master",
    olink,
    chr.s,
    pos.s,
    pos.e,
    "z",
    sep = "."
  ),
  sep = ";",
  row.names = FALSE,
  quote = FALSE
)

setDT(res.olink)

## create identifier column in results to keep the mapping to the LD matrix
res.olink[, snp.id := 1:nrow(res.olink)]

## order, keep in mind that cs -1 means not included in the credible set
res.olink <- res.olink[order(-cs, -pip)]

## create indicator to select
res.olink[, ind := 1:.N, by = "cs"]

fwrite(
  res.olink,
  "res_olink.tsv",
  sep = "\t"
)

