#' Simulate & test mediation methods on omics dataset
#'
#' Simulate & test mediation methods on omics dataset.
#'
#' @param b1t2 Numeric value that both theta2 (\code{"t2"}) and beta1 (\code{"b1"}) take for non-null mediators.
#' @param t1 Numeric value of theta2 i.e. the effect of the exposure on the outcome.
#' @param nsamp Number of samples.
#' @param ngene Number of genes other than that of primary interest to simulate.
#' @param FDR FDR threshold to apply.
#' @param rho Inter-gene correlation coefficients.
#' @param prop.consistent Proportion of genes that are consistent mediators. Must be at least \code{1/ngene}.
#' @param prop.inconsistent Proportion of genes that are inconsistent mediators.
#' @inheritParams ezlimma:::sim_fisher
#' @return Matrix with proportion of significant calls for every method for true and null mediators.

sim_omics <- function(b1t2=1, t1=5, nsamp=15, ngene=100, FDR=0.25, rho=0, prop.consistent=1/ngene,
                      prop.inconsistent=0, nsim=10**3, seed=0, verbose=TRUE){

  stopifnot(prop.consistent >= 1/ngene, prop.inconsistent >= 0, prop.consistent <= 1, prop.inconsistent <= 1,
            prop.consistent + prop.consistent <= 1)
  set.seed(seed)
  #t = theta; b = beta
  t0 <- t3 <- b0 <- b2 <- 0.14

  prop.sig.arr <- array(NA, dim=c(3, 3, nsim),
                        dimnames=list(c("hitman", "lotman", "js"), c("n.hits", "power", "false_disc"), paste0("sim_", 1:nsim)))

  n.consistent <- round(prop.consistent*ngene)
  n.inconsistent <- round(prop.inconsistent*ngene)

  for (sim in 1:nsim){
    g.nms <- paste0("gene", 1:ngene)
    consistent_genes <- sample(g.nms, size=n.consistent)
    # returns character(0) if n.inconsistent=0
    inconsistent_genes <- sample(setdiff(g.nms, consistent_genes), size=n.inconsistent)
    med_genes <- union(consistent_genes, inconsistent_genes)

    x <- stats::rnorm(n=nsamp)
    a <- stats::rnorm(n=nsamp)

    Sigma <- matrix(data=rho, nrow=ngene, ncol=ngene)
    diag(Sigma) <- 1
    error_m <- t(MASS::mvrnorm(n=nsamp, mu = rep(0, ngene), Sigma = Sigma))
    rm(Sigma)

    m <- b0 + b2*x + error_m
    dimnames(m) <- list(g.nms, paste0("s", 1:nsamp))
    m[med_genes,] <- m[med_genes,, drop=FALSE] + b1t2*a

    # eq 3; E(Y)
    y <- t0 + t1*a + t3*x + stats::rnorm(n=nsamp)
    y <- y + colSums(b1t2 * m[consistent_genes,, drop=FALSE])
    if (n.inconsistent > 0) y <- y - colSums(b1t2 %*% m[inconsistent_genes,, drop=FALSE])

    names(a) <- names(x) <- names(y) <- paste0("s", 1:nsamp)

    hm.res <- hitman(E=a, M=m, Y=y, covariates = x, verbose=TRUE)
    lm.res <- lotman(E=a, M=m, Y=y, covariates = x, verbose = TRUE)

    js.v <- apply(m, 1, FUN=function(m.v){
      joint_signif_mediation(E=a, M=m.v, Y=y, covariates = x)
    })
    js.res <- data.frame(p=js.v, FDR=stats::p.adjust(js.v, method="fdr"))

    prop.sig.arr["hitman", "n.hits", sim] <- sum(hm.res$EMY.FDR < FDR)
    prop.sig.arr["hitman", "power", sim] <- mean(hm.res[consistent_genes, "EMY.FDR"] < FDR)
    prop.sig.arr["hitman", "false_disc", sim] <-
      sum(hm.res[setdiff(g.nms, consistent_genes), "EMY.FDR"] < FDR)

    prop.sig.arr["lotman", "n.hits", sim] <- sum(lm.res[, "EMY.FDR"] < FDR)
    prop.sig.arr["lotman", "power", sim] <- mean(lm.res[consistent_genes, "EMY.FDR"] < FDR)
    prop.sig.arr["lotman", "false_disc", sim] <-
      sum(lm.res[setdiff(g.nms, consistent_genes), "EMY.FDR"] < FDR)

    prop.sig.arr["js", "n.hits", sim] <- sum(js.res[, "FDR"] < FDR)
    prop.sig.arr["js", "power", sim] <- mean(js.res[consistent_genes, "FDR"] < FDR)
    prop.sig.arr["js", "false_disc", sim] <-
      sum(js.res[setdiff(g.nms, consistent_genes), "FDR"] < FDR)

    if (verbose && sim %% 10 == 0) message("sim: ", sim)
  }
  ret <- data.frame(apply(prop.sig.arr, 1:2, FUN=mean))
  ret$fdr <- as.numeric(apply(prop.sig.arr[,"false_disc",, drop=FALSE], c(1, 2), FUN=sum)/
    apply(prop.sig.arr[,"n.hits",, drop=FALSE], 1:2, FUN=sum))
  ret
}
