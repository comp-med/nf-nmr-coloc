#! /usr/bin/env Rscript

## get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

region          <- args[1]
res_finemap_file <- args[2]
ld_file        <- args[3]
snplist_file  <- args[4]
r_lib          <- args[5]

library(data.table, lib.loc = r_lib)

## Load data from channel paths
res_finemap <- fread(res_finemap_file, na.strings = "")
ld <- fread(ld_file, data.table = FALSE, na.strings = "")

# For each phenotype ----
res_finemap_all_phenos <- lapply(unique(res_finemap$pheno), function(p) {
  # p <- unique(res_finemap$pheno)[1]
  
  res_finemap_pheno <- res_finemap[pheno == (p)]
  
  ## get top SNPs
  top_snp_finemap <- res_finemap_pheno[ind == 1 & cs > 0, ]
  
  ## add LD columns: be careful 'cs' does not mean that numbers match,
  ## but is a legacy from fine-map
  for (j in sort(top_snp_finemap$cs)) {
    res_finemap_pheno[, paste0("R2.", j)] <- ld[
      res_finemap_pheno$snp_id, 
      top_snp_finemap$snp_id[which(top_snp_finemap$cs == j)
      ]
    ]^2
  }
  
  ## get all proxy SNPs (LD >0.8)
  proxy_snps <- paste0("R2.", sort(top_snp_finemap$cs))
  proxy_snps <- apply(res_finemap_pheno[, ..proxy_snps], 1, function(x) {
    ifelse(sum(x >= .8) > 0, T, F)
  })
  ## may include SNPs not in the credible set
  proxy_snps <- res_finemap_pheno[proxy_snps]
  
  ## store the LD pattern across top SNPs
  ## (convert to data frame to ease downstream operations)
  ld_top_snp_finemaps <- as.data.frame(res_finemap_pheno)
  ld_top_snp_finemaps <- lapply(top_snp_finemap$cs, function(x) {
    
    ## get all SNPs and corresponding LD
    tmp <- paste0("R2.", x)
    tmp <- ld_top_snp_finemaps[
      which(ld_top_snp_finemaps[, tmp] >= .8),
      c("id", tmp)
    ]
    
    ## edit names
    names(tmp) <- c("id_proxy", "R2")
    
    ## add top SNP
    tmp <- merge(
      as.data.frame(top_snp_finemap[x, c("id", "cs")]),
      tmp,
      suffix = c("_lead", "_proxy")
    )
    
    ## do some renaming to ease downstream coding
    names(tmp) <- c(
      "id_lead", "cs", "id_proxy", "R2_proxy"
    )
    
    ## return
    return(tmp)
  })
  
  ## combine everything
  ld_top_snp_finemaps <- rbindlist(
    ld_top_snp_finemaps, 
    use.names = TRUE, 
    fill = TRUE
  )
  
  ld_top_snp_finemaps[, pheno := p]
  
  return(
    list(
      res_finemap_pheno,
      top_snp_finemap,
      ld_top_snp_finemaps
    )
  )
})

res_finemap <- rbindlist(lapply(res_finemap_all_phenos, `[[`, 1), fill = TRUE)
top_snp_finemap <- rbindlist(lapply(res_finemap_all_phenos, `[[`, 2))
ld_top_snp_finemaps <- rbindlist(lapply(res_finemap_all_phenos, `[[`, 3))
ld_top_snp_finemaps <- ld_top_snp_finemaps[, .(id_proxy, pheno)]
ld_top_snp_finemaps <- ld_top_snp_finemaps[ , .(id_proxy = unique(id_proxy)), pheno]

fwrite(
  res_finemap,
  "res_finemapping_with_proxies.tsv", 
  sep = "\t"
) 

## write SNPs to file to be queried
fwrite(
  top_snp_finemap,
  "top_snp.tsv", 
  sep = "\t"
) 

## write SNPs to file to be queried
fwrite(
  ld_top_snp_finemaps,
  "proxy_snplist.tsv", 
  sep = "\t"
)