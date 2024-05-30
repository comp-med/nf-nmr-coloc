runsusie <- function (d, suffix = 1, maxit = 100, repeat_until_convergence = TRUE, 
          s_init = NULL, ...) 
{
check_dataset(d, suffix, req = c("beta", "varbeta", "LD", 
                                 "snp", "N"))
check_ld(d, d$LD)
if (!("z" %in% names(d))) 
  z = d$beta/sqrt(d$varbeta)
else z = d$z
LD = d$LD[d$snp, d$snp, drop = FALSE]
names(z) = d$snp
snp = d$snp
converged = FALSE
susie_args = list(...)
if ("max_iter" %in% names(susie_args)) {
  maxit = susie_args$max_iter
  susie_args = susie_args[setdiff(names(susie_args), "max_iter")]
}
if (!("n" %in% names(susie_args))) 
  susie_args = c(list(n = d$N), susie_args)
while (!converged) {
  message("running max iterations: ", maxit)
  res = do.call(susie_rss, c(list(z = z, R = LD, max_iter = maxit), 
                             susie_args))
  converged = res$converged
  message("\tconverged: ", converged)
  if (!converged && repeat_until_convergence == FALSE) 
    stop("susie_rss() did not converge in ", maxit, 
         " iterations. Try running with run_until_convergence=TRUE")
  if (!converged) 
    maxit = maxit * 100
}
res = annotate_susie(res, snp, LD)
res
}