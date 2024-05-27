#! /usr/bin/env Rscript

####################################################
## function to compute colocalisation assuming one
## causal variant

naive.coloc <- function(
  res_all,
  ld,
  x,
  y, 
  x_is_binary,
  nmr_region,
  res_finemapping,
  biomart_gene_annotation,
  rec_file
  # other files from plot 1 & 2
) {
  ## 'res_all' -- data set containing merged and aligned statistics
  ## 'ld'      -- corresponding LD matrix

  #-----------------------------------------#
  ## -- 	         sanity check            --##
  #-----------------------------------------#

  ## top signal for the protein in all data
  ts <- res_finemapping$snp_id[which.max(abs(res_finemapping$beta_marginal / res_finemapping$se_marginal))]
  ## top signal for protein in the overlap
  is <- res_all$snp_id[which.max(abs(res_all$beta_marginal / res_all$se_marginal))]
  ## get the top SNP for the outcome
  io <- res_all$snp_id[
    which.max(abs(res_all$Effect_outcome / res_all$StdErr_outcome))
  ]

  ## conserved signal for phecode
  ld.sens <- ld[ts, is]^2
  ## ld between lead signals
  ld.ovl <- ld[is, io]^2

  #-----------------------------------------#
  ## -- 	            run coloc            --##
  #-----------------------------------------#

  ## order by position
  setDT(res_all)
  res_all <- res_all[order(genpos)]

  ## prepare input
  D1 <- list(
    beta = res_all$beta_marginal,
    varbeta = res_all$se_marginal^2,
    type = "quant",
    sdY = 1,
    N = 434646, # This is fixed for each metabolite
    snp = res_all$snp_id,
    position = 1:nrow(res_all)
  )
  
  ## binary outcome
  if (isTRUE(x_is_binary)) {
    
    # needs: N, s
    
    D2 <- list(
      beta = res_all$Effect_outcome,
      varbeta = res_all$StdErr_outcome^2,
      type = "cc",
      N = max(res_all$N),
      s = max(res_all$outcome_n_cases) / max(res_all$N),
      snp = res_all$snp_id,
      position = 1:nrow(res_all)
    )
    
  ## Quantitative Outcome
  } else if (isFALSE(x_is_binary)) {
    
    if (is.null(res_all$FRQ)) {
      return(list())
    }
    
    D2 <- list(
      beta = res_all$Effect_outcome,
      varbeta = res_all$StdErr_outcome^2,
      type = "quant",
      N = max(res_all$N),
      MAF = res_all$FRQ,
      snp = res_all$snp_id,
      position = 1:nrow(res_all)
    )
    
  } else {
    stop("[ERROR] Colocalization outcome must be annotated as `binary` or `quantitative`")
  }

  ## do naive coloc as well
  naive.coloc <- coloc.signals(D1, D2, method = "single", p12 = 5e-6)

  ## add checks to the data
  naive.coloc$summary$R2.1 <- ld.sens
  naive.coloc$summary$ld.top <- ld.ovl

  ## add the trait id and label
  naive.coloc$summary$outcome.id <- x
  naive.coloc$summary$exposure <- y
  naive.coloc$summary$region <- nmr_region
  naive.coloc$summary$olink <- paste(nmr_region, y, sep = ": ")

  
  #-----------------------------------------#
  ## --         draw selected            --##
  #-----------------------------------------#

  if (naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.top > .8) {

    png(
      paste0(paste("graphics/coloc",nmr_region,y,x, sep = "_"), ".png"),
      width = 16,
      height = 8,
      units = "cm",
      res = 300
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
      sum.stat = res_all,
      sum.coloc = naive.coloc$summary,
      ld = ld,
      a.vars = res_all$id[
        which(
          res_all$snp_id == naive.coloc$summary$best4
        )
      ],
      biomart_gene_annotation = biomart_gene_annotation,
      rec_file = rec_file
    )
    dev.off()
  }
  
  ## add effect estimate top coloc snp
  naive.coloc <- as.data.table(naive.coloc$summary)
  naive.coloc[, best4 := as.numeric(best4)]
  
  # this was only needed for the plot
  naive.coloc$olink <- NULL
  
  naive.coloc <- merge(
    naive.coloc, 
    res_all[, c(
      "snp_id", "id", "marker_name", "allele1", "allele0",
      "pval_marginal", "beta_marginal", "se_marginal",
      "pip", "cs", "Effect_outcome",
      "StdErr_outcome", "P"
    )],
    by.x = "best4",
    by.y = "snp_id"
  )
  
  setnames(
    naive.coloc, 
    old = c(
      "best4", "hit2", "hit1", "nsnps", "PP.H0.abf", "PP.H1.abf", 
      "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "best1", "best2", "hit1.margz", 
      "hit2.margz", "R2.1", "ld.top", "outcome.id", "exposure", "region", 
      "id", "marker_name", "allele1", "allele0", "pval_marginal", "beta_marginal", 
      "se_marginal", "pip", "cs", "Effect_outcome", "StdErr_outcome", "P"),
    new = c(
      "snp_id_h4",
      "snp_id_h2",
      "snp_id_h1",
      "n_snps",
      "pp_h0_abf",
      "pp_h1_abf",
      "pp_h2_abf",
      "pp_h3_abf",
      "pp_h4_abf",
      "snp_id_exposure",
      "snp_id_outcome",
      "zscore_exposure",
      "zscore_outcome",
      "r2_1",
      "ld_top",
      "outcome",
      "exposure",
      "region",
      "id_exposure",
      "marker_name",
      "allele1",
      "allele0",
      "pval_exposure",
      "beta_exposure",
      "se_exposure",
      "pip",
      "cs",
      "beta_outcome",
      "se_outcome",
      "pval_outcome"
    )
    )

  ## write results to file
  return(naive.coloc)
}
