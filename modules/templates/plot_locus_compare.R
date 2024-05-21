#! /usr/bin/env Rscriptplot

###########################################
## function to plot stacked locus zoom plot
## colouring for different variants

plot.locus.compare <- function(
  sum.stat,
  sum.coloc,
  ld,
  a.vars = NA,
  biomart_gene_annotation,
  rec_file
) {
  ## 'sum.stat'  -- summary statistics for plotting
  ## 'sum.coloc' -- summary from coloc
  ## 'ld'        -- LD matrix (created based on LD-store)
  ## 'a.vars'    -- list of rsids to highlight and colour by

  ## package for faster data reading
  require(data.table)

  #-----------------------#
  ## --  perpare input  --##
  #-----------------------#

  ## convert to data frame if needed
  sum.stat <- as.data.frame(sum.stat)

  ## create log10p
  if ("LOG10P" %in% names(sum.stat)) {
    sum.stat$log10p.olink <- sum.stat$LOG10P
  } else {
    sum.stat$log10p.olink <- -pchisq(
      (sum.stat[, "Effect"] / sum.stat[, "StdErr"])^2,
      df = 1,
      lower.tail = FALSE,
      log.p = TRUE
    ) / log(10)
  }
  ## trait
  sum.stat$log10p.trait <- -pchisq(
    (sum.stat[, "Effect.outcome"] / sum.stat[, "StdErr.outcome"])^2,
    df = 1,
    lower.tail = FALSE,
    log.p = TRUE
  ) / log(10)

  ## delete possible previously computed LD values
  sum.stat[, grep("R2", names(sum.stat))] <- NULL

  ## loop through all SNPs shared
  for (j in 1:length(a.vars)) {
    ## depending on how the matrix is named
    if (paste0("V", sum.stat$snp.rsid[1]) %in% colnames(ld)) {
      ## map variant to rsid
      ii <- sum.stat$snp.rsid[which(sum.stat$rsid == a.vars[j])]
      ## add to the data
      sum.stat[, paste0("R2.", j)] <- ld[sum.stat$snp.rsid, ii]^2
    } else {
      ## map variant to rsid
      ii <- sum.stat$snp.id[which(sum.stat$rsid == a.vars[j])]
      ## add to the data
      sum.stat[, paste0("R2.", j)] <- ld[sum.stat$snp.id, ii]^2
    }
  }

  ## create a couple of colour gradients
  cls <- c("#E58606", "#5D69B1", "#52BCA3", "#99C945",
    "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F",
    "#ED645A", "#CC3A8E", "#A5AA99", "#00A4CC",
    "#F95700", rainbow(10)
  )

  ## create colouts list
  col.list <- lapply(1:20, function(x) {
    colorRampPalette(c("white", cls[x]))(100)
  })

  ## load recombination rates
  rec <- data.frame(
    fread(rec_file, sep = "\t", header = TRUE)
  )

  ## subset to region of interest
  rec <- subset(
    rec, 
    Position.bp. >= min(sum.stat$pos, na.rm = TRUE) &
      Position.bp. <= max(sum.stat$pos, na.rm = TRUE)
  )

  #------------------------------------#
  ## --        oppose p-values       --##
  #------------------------------------#

  plot(log10p.trait ~ log10p.olink,
    xlim = c(0, max(sum.stat$log10p.olink) * 1.05),
    data = sum.stat,
    cex = .4,
    col = "grey90",
    xlab = bquote(
      .(sum.coloc$phenotype.olink) ~ "GWAS" ~ ~ -log[10]("p-value")
    ),
    ylab = bquote(
      .(sum.coloc$phenotype.ieu) ~ "GWAS" ~ ~ -log[10]("p-value")
    )
  )

  ## add variants to be highlighted
  for (j in 1:length(a.vars)) {
    cat("\n plot LD with ", a.vars[j], "\n")

    ## rsidentify all up to R2>.3
    tmp <- sum.stat[which(sum.stat[, paste0("R2.", j)] >= .3), ]
    ## add points
    points(
      tmp$log10p.olink, 
      tmp$log10p.trait,
      cex = .5,
      pch = 21, 
      col = "grey10",
      bg = col.list[[j]][ceiling(tmp[, paste0("R2.", j)] * 100)],
      lwd = .2
    )

    ## add lead variant as diamond
    tmp <- subset(sum.stat, rsid == a.vars[j])
    print(tmp)
    # print(tmp)
    points(
      tmp$log10p.olink,
      tmp$log10p.trait,
      pch = 23, 
      lwd = .8, 
      cex = .7,
      col = col.list[[j]][ceiling(tmp[, paste0("R2.", j)] * 100)],
      bg = "white",
      type = "p"
    )

    ## annotate
    text(tmp$log10p.olink, tmp$log10p.trait * 1.05,
      labels = a.vars[j], cex = .4, pos = 4, xpd = NA, offset = .15
    )
    # print(tmp)
    ## add arrow to combine both
    arrows(
      tmp$log10p.olink, 
      tmp$log10p.trait,
      tmp$log10p.olink,
      tmp$log10p.trait * 1.05,
      lwd = .5, 
      length = 0,
      xpd = NA
    )
  }

  ## add legend for Coloc
  if (nrow(sum.coloc) == 1) {
    legend(
      ifelse(sum.coloc[1, 8] > .5,
        "topleft",
        "topright"
      ),
      lty = 0, 
      pch = NA,
      cex = .5,
      legend = paste(
        paste0("H", 1:4),
        "=",
        sprintf("%.1f", sum.coloc[1, 5:8] * 100)
      ),
      title = "Posterior prob. [%]"
    )
  } else {
    ## show only evrsidence for H4
    legend("topleft",
      lty = 0, pch = NA, cex = .5,
      legend = paste(
        sum.coloc$rsid.olink,
        "-",
        sum.coloc$rsid.ieu,
        "=", 
        sprintf("%.1f", sum.coloc$PP.H4.abf * 100)
      ),
      title = "Posterior prob. [%]: H4"
    )
  }


  #-----------------------#
  ## --    plot stats   --##
  #-----------------------#

  ## plot the SOMAscan
  par(mar = c(.1, 1.5, 1, .5))

  ## define traits to be plotted
  tr <- c("olink", "trait")

  ## now draw the plot in a loop
  for (k in 1:length(tr)) {
    cat("\n plot trait ", tr[k], "\n")

    ## add all variants
    plot(
      sum.stat[, "pos"], 
      sum.stat[, paste("log10p", tr[k], sep = ".")],
      xlab = "", ylab = expression(-log[10]("p-value")),
      ylim = c(
        0, 
        max(sum.stat[, paste("log10p", tr[k], sep = ".")], na.rm = T) * 1.1),
      cex = .3, 
      col = "grey90",
      xaxt = "n", 
      yaxt = "n"
    )
    axis(2, lwd = .5)

    #------------------------------------#
    ## -- colour bars for each variant --##
    #------------------------------------#

    ## add each variant
    for (j in 1:length(a.vars)) {
      cat("\n plot LD with ", a.vars[j], "\n")

      ## rsidentify all up to R2>.3
      tmp <- sum.stat[which(sum.stat[, paste0("R2.", j)] >= .3), ]
      ## add points
      points(tmp[, "pos"], tmp[, paste("log10p", tr[k], sep = ".")],
        cex = .5, pch = 21, col = "grey10",
        bg = col.list[[j]][ceiling(tmp[, paste0("R2.", j)] * 100)],
        lwd = .2
      )

      ## add lead variant as diamond
      tmp <- subset(sum.stat, rsid == a.vars[j])
      print(tmp)
      points(
        tmp[, "pos", drop = F], 
        tmp[, paste("log10p", tr[k], sep = "."), drop = F],
        pch = 23,
        lwd = .5,
        cex = .7,
        col = col.list[[j]][ceiling(tmp[, paste0("R2.", j)] * 100)],
        bg = "white",
        type = "p"
      )

      ## annotate
      text(
        tmp[, "pos", drop = F], 
        tmp[, paste("log10p", tr[k], sep = "."), drop = F] * 1.05,
        labels = a.vars[j],
        cex = .4,
        pos = 4,
        xpd = NA,
        offset = .15
      )
      print(tmp)
      ## add arrow to combine both
      arrows(tmp[, "pos"], tmp[, paste("log10p", tr[k], sep = ".")],
        tmp[, "pos"], tmp[, paste("log10p", tr[k], sep = ".")] * 1.05,
        lwd = .5, length = 0, xpd = NA
      )

      ## plotting coordinates
      pm <- par("usr")

      ## only in upper most plot
      if (k == 1) {
        ## add colour gradient legend
        l <- seq(
          pm[1] + (pm[2] - pm[1]) * .75, 
          pm[1] + (pm[2] - pm[1]) * .95, 
          length.out = 100)
        ## rectangle for the colours
        rect(
          l - (l[2] - l[1]) / 2,
          pm[3] + (pm[4] - pm[3]) * (.93 - (j * .08)),
          l + (l[2] - l[1]) / 2,
          pm[3] + (pm[4] - pm[3]) * (1 - (j * .08)),
          border = NA,
          col = col.list[[j]]
        )
        ## box
        rect(
          l[1] - (l[2] - l[1]) / 2, 
          pm[3] + (pm[4] - pm[3]) * (.93 - (j * .08)), 
          l[100] + (l[2] - l[1]) / 2,
          pm[3] + (pm[4] - pm[3]) * (1 - (j * .08)),
          border = "black",
          col = NA,
          lwd = .3
        )
        ## add header
        text(
          pm[1] + (pm[2] - pm[1]) * .75,
          pm[3] + (pm[4] - pm[3]) * (.97 - (j * .08)),
          cex = .4, 
          labels = a.vars[j],
          pos = 4,
          offset = .2
        )

        if (j == 1) {
          text(
            pm[1] + (pm[2] - pm[1]) * .75,
            pm[3] + (pm[4] - pm[3]) * .96,
            cex = .5,
            labels = "r2 with", 
            pos = 4,
            offset = .2
          )
        }

        if (j == length(a.vars)) {
          ## simple axis
          text(
            l[c(1, 20, 40, 60, 80, 100)],
            pm[3] + (pm[4] - pm[3]) * (.93 - (j * .08)), 
            labels = c(0, .2, .4, .6, .8, 1), 
            pos = 1,
            cex = .4, 
            offset = .1
          )
        }
      }
    }

    ## add trait
    if (tr[k] == "olink") {
      legend(
        "topleft", 
        pch = NA, 
        lty = 0,
        cex = .5,
        legend = sum.coloc$olink[1],
        bty = "n"
      )
    } else {
      legend(
        "topleft",
        pch = NA,
        lty = 0,
        cex = .5,
        legend = sum.coloc$outcome.id[1],
        bty = "n"
      )
    }
  }

  #-----------------------#
  ## --    plot genes   --##
  #-----------------------#

  cat("\n plot genes ", a.vars[j], "\n")

  #------------------------------------#
  ## --        gene assignment       --##
  #------------------------------------#

  ## import entire list, given that biomart acts weirdly
  tmp.genes <- fread(biomart_gene_annotation)
  ## subset to what is needed
  tmp.genes <- subset(
    tmp.genes,
    chromosome_name == ifelse(sum.stat$chr[1] == 23, "X", sum.stat$chr[1]) &
      start_position >= min(sum.stat$pos, na.rm = T) - 3e6 &
      end_position <= max(sum.stat$pos, na.rm = T) + 3e6
  )

  ## restrict to protein encoding genes for now
  if (nrow(tmp.genes) > 20) {
    tmp.genes.subset <- subset(
      tmp.genes,
      gene_biotype %in% c("protein_coding", "processed_transcript")
    )
    if (nrow(tmp.genes.subset) > 0) {
      tmp.genes <- tmp.genes.subset
    }
  }
  ## sort by start
  tmp.genes <- tmp.genes[order(tmp.genes$start_position), ]

  print(tmp.genes)

  ## dummy for the line in the plot
  tmp.genes$line <- NA
  tmp.genes$line[1] <- 1

  ## start sorting
  l <- 0

  ## loop over everything, collect genes row-wise
  while (sum(is.na(tmp.genes$line)) > 0) {
    ## increase line
    l <- l + 1
    e <- 0
    for (j in 1:nrow(tmp.genes)) {
      ## test whether new start is larger than the current end
      if (tmp.genes$start_position[j] > e + 4e4 & is.na(tmp.genes$line[j])) {
        ## assign line to be drawn
        tmp.genes$line[j] <- l
        ## assign new end
        e <- tmp.genes$end_position[j]
      }
    }
  }

  print(head(tmp.genes))

  ## now plot it just below
  par(mar = c(1.5, 1.5, .1, .5), tck = -.02, bty = "o", lwd = .5)
  plot(range(sum.stat$pos), c(0, max(tmp.genes$line) + .5),
    yaxt = "n", xaxt = "n", ylab = "",
    xlab = paste("Genomic position on chrosome", sum.stat$chr[1]), type = "n",
    ylim = rev(c(0, max(tmp.genes$line) + .5))
  )
  axis(1, lwd = .5)
  ## position of the lead variant
  for (j in 1:length(a.vars)) {
    ii <- which(sum.stat$rsid == a.vars[j])
    abline(v = sum.stat$pos[ii], lwd = .3, col = cls[j])
  }
  ## add genes
  arrows(
    tmp.genes$start_position,
    tmp.genes$line, 
    tmp.genes$end_position, 
    tmp.genes$line, lwd = .5, 
    length = 0
  )
  ## add Gene Names on top
  text(tmp.genes$start_position + (
    tmp.genes$end_position - tmp.genes$start_position) / 2, tmp.genes$line,
    cex = .3, font = 3,
    labels = tmp.genes$external_gene_name, pos = 3, offset = .1
  )
}
