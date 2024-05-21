#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

olink          <- args[1]
ld_file        <- args[2]
res_olink_file <- args[3]
r_lib          <- args[4]

library(data.table, lib.loc = r_lib)

## Load data from channel paths
ld <- fread(ld_file, data.table = FALSE, na.strings = "")
res.olink <- fread(res_olink_file, na.strings = "")

## get top SNPs

top.snp <- res.olink[ind == 1 & cs > 0]

## add LD columns: be careful 'cs' does not mean that numbers match,
## but is a legacy from fine-map
for (j in top.snp$cs) {
  res.olink[, paste0("R2.", j)] <- ld[
    res.olink$snp.id,
    top.snp$snp.id[which(top.snp$cs == j)]]^2
}

## get all proxy SNPs
proxy.snps <- paste0("R2.", top.snp$cs)
proxy.snps <- apply(res.olink[, ..proxy.snps], 1, function(x) {
  ifelse(sum(x >= .8) > 0, T, F)
})
## may include SNPs not in the credible set
proxy.snps <- res.olink[proxy.snps]

## store the LD pattern across top SNPs
## (convert to data frame to ease downstream operations)
ld.top.snps <- as.data.frame(res.olink)
ld.top.snps <- lapply(top.snp$cs, function(x) {

  ## get all SNPs and corresponding LD
  tmp <- paste0("R2.", x)
  tmp <- ld.top.snps[
    which(ld.top.snps[, tmp] >= .8),
    c("MarkerName", "id", "rsid", tmp)
  ]

  ## edit names
  names(tmp) <- c("MarkerName.proxy", "id.proxy", "rsid.proxy", "R2")

  ## add top SNP
  tmp <- merge(
    as.data.frame(top.snp[x, c("MarkerName", "id", "rsid", "cs")]),
    tmp,
    suffix = c(".lead", ".proxy")
  )

  ## do some renaming to ease downstream coding
  names(tmp) <- c(
    "MarkerName.lead", "id.lead", "rsid.lead",
    "cs", "MarkerName.proxy", "id.proxy",
    "rsid.proxy", "R2.proxy"
  )

  ## return
  return(tmp)
})

## combine everything
ld.top.snps <- do.call(rbind, ld.top.snps)

fwrite(
  res.olink,
  "res_olink_with_proxies.tsv", 
  sep = "\t"
) 

## write SNPs to file to be queried
fwrite(
  top.snp,
  "top_snp.tsv", 
  sep = "\t"
) 

## write SNPs to file to be queried
write.table(
  unique(c("SNP", ld.top.snps$rsid.proxy)),
  paste0(
    olink,
    ".snplist"),
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
