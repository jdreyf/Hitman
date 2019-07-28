#' Apply Fisher test to simulated data analyzed by \code{hitman}
#'
#' Apply Fisher test to simulated data analyzed by \code{hitman}.
#'
#' @param effect.v Numeric vector of log fold-changes or percent of phentotypes to add.
#' @param all.feats All feature names. Must overlap with \code{G}.
#' @param prop.other.sig Proportion of genes outside set to have signal added.
#' @inheritParams ezlimma::roast_contrasts
#' @inheritParams hitman
#' @inheritParams sim_barfield

# can't test in same fcn as sim_pants_hitman b/c need to add ev to some subset of genes
# competitive, so need to decide how much signal to give other genes
sim_fisher_hitman <- function(G, E, all.feats, effect.v=c(0, 0.2), alpha=0.05, nsim=99, prop.other.sig=0.01,
                              seed=1, verbose=TRUE){
  g1 <- G[[1]]$genes
  stopifnot(any(g1 %in% all.feats))

  sig.mat <- matrix(NA, nrow=nsim, ncol=length(effect.v),
                         dimnames=list(paste0("sim", 1:nsim), paste0("eff_", effect.v)))
  set.seed(seed)
  err <- stats::rnorm(n=length(E), sd=stats::sd(E)/2)
  phenotype <- E + err

  other.g <- setdiff(all.feats, g1)

  for (sim in 1:nsim){
    for (ev in effect.v){
      obj.test <- matrix(stats::rnorm(n=length(all.feats)*length(E)), ncol=length(E), nrow=length(all.feats),
                         dimnames=list(all.feats, names(E)))
      if (ev > 0){
        # to geneset & 1% of others
        other.sig <- sample(x=other.g, size=prop.other.sig*length(other.g))
        obj.test[c(g1, other.sig),] <- t(t(obj.test[c(g1, other.sig),]) + phenotype + ev*err)
      }

      hm.tab <- hitman(M=obj.test, E=E, Y=phenotype)[, 1, drop=FALSE]
      sig.genes <- rownames(hm.tab)[hm.tab[,1] < alpha]
      ftp <- stats::fisher.test(all.feats %in% sig.genes, all.feats %in% g1, alternative = "greater")$p.value
      sig.mat[sim, paste0("eff_", ev)] <- ftp < alpha
    }
    if (verbose & sim %% 100 == 0) cat("sim: ", sim, "\n")
  }
  (sig.mat <- rbind(avg=colMeans(sig.mat), sig.mat))
}
