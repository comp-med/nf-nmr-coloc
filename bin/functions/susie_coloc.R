#! /usr/bin/env Rscript

####################################################
## function to compute colocalisation assuming more
## than one causal variant with susie coloc

susie.coloc <- function(
  res.all,
  res.olink,
  ld,
  top.snp,
  x,
  y,
  nmr_region,
  x_is_binary,
  olink_id,
  biomart_gene_annotation,
  rec_file
) {
  ## 'res.all' -- data set containing merged and aligned statistics
  ## 'ld'      -- corresponding LD matrix
  ## 'top.snp' -- top snps for pQTL credible sets

  #-----------------------------------------#
  ## --         prepare coloc            --##
  #-----------------------------------------#

  ## order by position
  res.all <- as.data.table(res.all)
  res.all <- res.all[order(genpos)]

  ## edit names
  rownames(ld) <- colnames(ld)

  ## prepare input
  D1 <- list(
    beta = res.all$Effect,
    varbeta = res.all$StdErr^2,
    type = "quant",
    sdY = 1,
    N = 434646,
    snp = paste0("V", res.all$snp_id),
    position = 1:nrow(res.all),
    ## subset the LD matrix to what is needed
    LD = as.matrix(ld[paste0("V", res.all$snp_id), paste0("V", res.all$snp_id)])
  )

  ## binary outcome
  if (isTRUE(x_is_binary)) {
    
    ## binary outcome
    D2 <- list(
      beta = res.all$Effect_outcome, 
      varbeta = res.all$StdErr_outcome^2,
      type = "cc",
      N = max(res.all$N),
      s = max(res.all$outcome_n_cases) / max(res.all$N),
      snp = paste0("V", res.all$snp_id),
      position = 1:nrow(res.all),
      ## subset the LD matrix to what is needed
      LD = as.matrix(ld[paste0("V", res.all$snp_id), paste0("V", res.all$snp_id)])
    )
    
    ## Quantitative Outcome
  } else if (isFALSE(x_is_binary)) {
    
    if (is.null(res_all$FRQ)) {
      return(list())
    }
    
    D2 <- list(
      beta = res.all$Effect_outcome, 
      varbeta = res.all$StdErr_outcome^2,
      type = "quant",
      N = max(res.all$N),
      MAF = res.all$FRQ,
      snp = paste0("V", res.all$snp_id),
      position = 1:nrow(res.all),
      ## subset the LD matrix to what is needed
      LD = as.matrix(ld[paste0("V", res.all$snp_id), paste0("V", res.all$snp_id)])
    )
    
  } else {
    stop("[ERROR] Colocalization outcome must be annotated as `binary` or `quantitative`")
  }
  
  ## Olink protein
  set.seed(42)
  susie.olink <- tryCatch(
    {
      runsusie(
        D1,
        max_iter = 10000,
        repeat_until_convergence = FALSE,
        L = ifelse(nrow(top.snp) == 1, 2, nrow(top.snp))
      )
    },
    error = function(e) {
      message(e)
      }
  )

  ## trait
  susie.trait <- tryCatch(
    {
      runsusie(
        D2,
        # r2.prune = .25,
        max_iter = 10000,
        repeat_until_convergence = FALSE,
        L = 5
      )
    },
    error = function(e) {
      message(e)
    }
  )

  ## check whether both have outcome
  if (!(length(susie.olink) > 1 & 
        length(susie.trait) > 1)) {
    return(list())
  }
  
  ## additional check, whether both traits have at least
  ## one credible set to be tested
  if (!(length(summary(susie.olink)$cs) > 0 &
        length(summary(susie.trait)$cs) > 0
  )) {
    return(list())
  }
  
  #-----------------------------------------#
  ## --             run coloc            --##
  #-----------------------------------------#
  
  ## run coloc with susie input
  res.coloc <- coloc.susie(
    susie.olink, 
    susie.trait, 
    p12 = 5e-6
  )
  
  #-----------------------------------------#
  ## --  cross-check with fine-mapping   --##
  #-----------------------------------------#
  
  ## add LD with fine-mapped variants
  res.coloc <- as.data.table(res.coloc$summary)
  
  ## some processing to ease downstream mapping
  res.coloc[, hit1 := as.numeric(gsub("V", "", hit1))]
  res.coloc[, hit2 := as.numeric(gsub("V", "", hit2))]
  
  ## compute ld between selected lead hits
  res.coloc$ld.top <- apply(
    res.coloc[, c("hit1", "hit2"), drop = FALSE],
    1,
    function(k) ld[k[1], k[2]]^2
  )
  
  #-----------------------------------------#
  ## --     add additional information    --##
  #-----------------------------------------#
  
  ## add effect estimates (restrict to selected pQTLs, which means
  ## no effect estimates for the regional lead for the trait are
  ## taken forward)
  res.coloc <- merge(
    res.coloc,
    res.all[, c(
      "snp_id", "id", "allele1", "allele0",
      "pval_marginal", "beta_marginal", "se_marginal",
      "pip", "cs", "Effect_outcome",
      "StdErr_outcome", "P"
    )],
    by.x = "hit1", by.y = "snp_id"
  )
  
  ## add phecode
  res.coloc$olink <- paste(nmr_region, y, sep = ": ")
  res.coloc$outcome.id <- x
  res.coloc$exposure <- y
  res.coloc$region <- nmr_region
  
  ## add rsids and LD with other variants
  res.coloc <- merge(
    res.coloc,
    res.olink[, c("snp_id", "id", "marker_name")],
    by.x = "hit2",
    by.y = "snp_id",
    suffixes = c("_exposure", "_outcome")
  )
  
  #-----------------------------------------#
  ## --         draw selected            --##
  #-----------------------------------------#
  
  if (
    sum(res.coloc$PP.H4.abf > .7) > 0 |
    sum(res.coloc$ld.top > .8) > 0
  ) {
    png(
      paste0(paste("graphics/susie",nmr_region,y,x, sep = "_"), ".png"),
      width = 16,
      height = 8,
      units = "cm",
      res = 200
    )
    par(
      mar = c(1.5, 1.5, 1, .5),
      mgp = c(.6, 0, 0),
      cex.axis = .5,
      cex.lab = .5,
      tck = .01,
      cex.main = .6,
      font.main = 2
    )
    ## more complex layout for gene assignment
    layout(
      matrix(c(1, 1, 1, 2, 3, 4), 3, 2),
      heights = c(.43, .37, .2)
    )
    plot.locus.compare(
      sum.stat = res.all,
      sum.coloc = res.coloc,
      ld,
      a.vars = unique(
        c(res.coloc$id_exposure,
          res.coloc$id_outcome)
      ),
      biomart_gene_annotation = biomart_gene_annotation,
      rec_file = rec_file
    )
    dev.off()
  }
  
  # this was only needed for the plot
  res.coloc$olink <- NULL
  
  ## do some renaming
  setnames(
    res.coloc, 
    old = c(
      "hit2", "hit1", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", 
      "PP.H3.abf", "PP.H4.abf", "idx1", "idx2", "ld.top", "id_exposure", 
      "allele1", "allele0", "pval_marginal", "beta_marginal", "se_marginal", 
      "pip", "cs", "Effect_outcome", "StdErr_outcome", "P", "outcome.id", 
      "exposure", "region", "id_outcome", "marker_name"),
    new =  c(
      "snp_id_exposure",
      "snp_id_outcome", 
      "n_snps", 
      "pp_h0_abf",
      "pp_h1_abf",
      "pp_h2_abf",
      "pp_h3_abf",
      "pp_h4_abf",
      "idx1",
      "idx2",
      "ld_top",
      "id_exposure", 
      "allele1",
      "allele0",
      "pval_exposure",
      "beta_exposure",
      "se_exposure",
      "pip",
      "cs",
      "beta_outcome",
      "se_outcome",
      "pval_outcome",
      "outcome",
      "exposure",
      "region",
      "id_outcome",
      "marker_name"
    )
  )
  
  ## add R2 with top credible set
  res.coloc <- merge(
    res.coloc,
    res.all[, c(
      "marker_name",
      grep("R2\\.", names(res.all), value = TRUE)
    ), with = FALSE],
    by = "marker_name"
  )
  
  # change names to fit convention
  r2_cols <- grep("R2\\.", names(res.coloc), value = TRUE)
  setnames(
    res.coloc,
    old = r2_cols, 
    new = gsub(".", "_", tolower(r2_cols), fixed = TRUE)
  )
  
  ## write results to file
  return(res.coloc)
}
