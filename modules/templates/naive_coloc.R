#! /usr/bin/env Rscript

####################################################
## function to compute colocalisation assuming one
## causal variant

naive.coloc <- function(
  res.all,
  ld,
  x,
  x_is_binary,
  olink_id,
  res.olink,
  biomart_gene_annotation,
  rec_file
  # other files from plot 1 & 2
) {
  ## 'res.all' -- data set containing merged and aligned statistics
  ## 'ld'      -- corresponding LD matrix

  #-----------------------------------------#
  ## -- 	         sanity check            --##
  #-----------------------------------------#

  ## top signal for the protein in all data
  ts <- res.olink$snp.id[which.max(abs(res.olink$Effect / res.olink$StdErr))]
  ## top signal for protein in the overlap
  is <- res.all$snp.id[which.max(abs(res.all$Effect / res.all$StdErr))]
  ## get the top SNP for the outcome
  io <- res.all$snp.id[
    which.max(abs(res.all$Effect.outcome / res.all$StdErr.outcome))
  ]

  ## conserved signal for phecode
  ld.sens <- ld[ts, is]^2
  ## ld between lead signals
  ld.ovl <- ld[is, io]^2

  #-----------------------------------------#
  ## -- 	            run coloc            --##
  #-----------------------------------------#

  ## order by position
  res.all <- as.data.table(res.all)
  res.all <- res.all[order(pos)]

  ## prepare input
  D1 <- list(
    beta = res.all$Effect,
    varbeta = res.all$StdErr^2,
    type = "quant",
    sdY = 1,
    N = max(res.all$TotalSampleSize),
    MAF = res.all$MAF,
    snp = res.all$snp.id,
    position = 1:nrow(res.all)
  )
  
  ## binary outcome
  if (isTRUE(x_is_binary)) {
    
    # needs: N, s
    
    D2 <- list(
      beta = res.all$Effect.outcome,
      varbeta = res.all$StdErr.outcome^2,
      type = "cc",
      MAF = res.all$MAF,
      N = max(res.all$N),
      s = max(res.all$outcome_n_cases) / max(res.all$N),
      snp = res.all$snp.id,
      position = 1:nrow(res.all)
    )
    
  ## Quantitative Outcome
  } else if (isFALSE(x_is_binary)) {
    
    D2 <- list(
      beta = res.all$Effect.outcome,
      varbeta = res.all$StdErr.outcome^2,
      type = "quant",
      MAF = res.all$MAF,
      N = max(res.all$N),
      snp = res.all$snp.id,
      position = 1:nrow(res.all)
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
  naive.coloc$summary$olink <- olink_id

  #-----------------------------------------#
  ## --         draw selected            --##
  #-----------------------------------------#

  if (naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.top > .8) {

    png(
      paste0(
        "graphics/coloc.",
        x,
        ".",
        olink_id,
        ".png"),
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
      sum.stat = res.all,
      sum.coloc = naive.coloc$summary,
      ld = ld,
      a.vars = res.all$rsid[
        which(
          res.all$snp.id == naive.coloc$summary$best4
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
  naive.coloc <- merge(naive.coloc, res.all[, c(
    "snp.id", "MarkerName", "rsid", "EA", "NEA",
    "MAF", "Effect", "StdErr", "Pvalue",
    "pip", "cs", "Effect.outcome",
    "StdErr.outcome", "P"
  )],
  by.x = "best4", by.y = "snp.id"
  )

  ## do some renaming
  names(naive.coloc) <- c(
    "snp.id.H4", "snp.id.H2", "snp.id.H1", "nsnps",
    "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
    "snp.id.protein", "snp.id.outcome", "zscore.protein", "zscore.outcome",
    "R2.1", "ld.top", "outcome.id", "olink",
    "MarkerName", "rsid.protein", "EA", "NEA", "MAF",
    "Effect.protein", "StdErr.protein", "Pvalue.protein", "pip", "cs",
    "Effect.outcome", "StdErr.outcome", "Pvalue.outcome"
  )

  ## write results to file
  return(naive.coloc)
}
