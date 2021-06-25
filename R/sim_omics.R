#' Simulate & test mediation methods on omics dataset
#'
#' Simulate & test mediation methods on omics dataset.
#'
#' @param b1t2 Numeric value that both theta2 (\code{"t2"}) and beta1 (\code{"b1"}) take for non-null mediators.
#' @param t1 Numeric value of theta2 i.e. the effect of the exposure on the outcome.
#' @param nsamp Number of samples.
#' @param ngene Number of genes other than that of primary interest to simulate.
#' @param FDR FDR threshold to apply.
#' @param Sigma Matrix of inter-gene correlation coefficients.
#' @param prop.consistent Proportion of genes that are consistent mediators. Must be at least \code{1/ngene}.
#' @param prop.inconsistent Proportion of genes that are inconsistent mediators.
#' @param prop.1c Proportion of genes where exactly one component of the test (E --> M or M --> Y given E) holds.
#' @inheritParams ezlimma:::sim_fisher
#' @inheritParams hitman
#' @return Matrix with proportion of significant calls for every method for true and null mediators.
#' @export

sim_omics <- function(b1t2=1, t1=5, nsamp=15, ngene=100, FDR=0.25, Sigma=diag(ngene), prop.consistent=1/ngene,
                      prop.inconsistent=0, prop.1c=0, nsim=10**3, fdr.method=c("BH", "BY"), seed=0, verbose=TRUE){

  fdr.method <- match.arg(fdr.method, c("BH", "BY"))

  stopifnot(prop.consistent >= 1/ngene, prop.inconsistent >= 0, prop.consistent <= 1, prop.inconsistent <= 1,
            prop.consistent + prop.consistent + prop.1c <= 1, prop.1c >= 0)

  set.seed(seed)
  #t = theta; b = beta
  t0 <- t3 <- b0 <- b2 <- 0.14

  n.consistent <- round(prop.consistent*ngene)
  n.inconsistent <- round(prop.inconsistent*ngene)

  cnms <- c("n.hits", "power", "false_disc", "fdr")

  if (n.inconsistent > 0){
    prop.sig.arr <- array(NA, dim=c(4, length(cnms), nsim),
                          dimnames=list(c("hitman", "lotman", "js", "js_both"), cnms, paste0("sim_", 1:nsim)))
  } else {
    prop.sig.arr <- array(NA, dim=c(3, length(cnms), nsim),
                          dimnames=list(c("hitman", "lotman", "js"), cnms, paste0("sim_", 1:nsim)))
  }

  g.nms <- paste0("gene", 1:ngene)

  for (sim in 1:nsim){
    consistent_genes <- sample(g.nms, size=n.consistent)
    # returns character(0) if n.inconsistent=0
    inconsistent_genes <- sample(setdiff(g.nms, consistent_genes), size=n.inconsistent)
    med_genes <- union(consistent_genes, inconsistent_genes)
    # one component genes
    one_comp_genes_em <- sample(x=setdiff(g.nms, med_genes), size=ceiling(prop.1c*ngene/2))
    one_comp_genes_my <- sample(setdiff(g.nms, union(med_genes, one_comp_genes_em)), size=floor(prop.1c*ngene/2))

    # x is covariate
    x <- stats::rnorm(n=nsamp)
    # a is exposure
    a <- stats::rnorm(n=nsamp)

    error_m <- t(MASS::mvrnorm(n=nsamp, mu = rep(0, ngene), Sigma = Sigma))

    m <- b0 + matrix(b2*x, nrow=ngene, ncol=nsamp, byrow = TRUE) + error_m
    dimnames(m) <- list(g.nms, paste0("s", 1:nsamp))
    # contribution from exposure
    m[med_genes,] <- m[med_genes,, drop=FALSE] + b1t2*a
    if (length(one_comp_genes_em) > 0) m[one_comp_genes_em,] <- m[one_comp_genes_em,, drop=FALSE] + b1t2*a

    # outcome: modified by mediator genes
    y <- t0 + t1*a + t3*x + stats::rnorm(n=nsamp)
    y <- y + colSums(b1t2 * m[consistent_genes,, drop=FALSE])
    if (length(one_comp_genes_my) > 0) y <- y + colSums(b1t2 * m[one_comp_genes_my,, drop=FALSE])
    if (n.inconsistent > 0) y <- y - colSums(b1t2 * m[inconsistent_genes,, drop=FALSE])

    names(a) <- names(x) <- names(y) <- paste0("s", 1:nsamp)

    hm.res <- hitman(E=a, M=m, Y=y, covariates = x, fdr.method = fdr.method, verbose=TRUE)
    lm.res <- lotman(E=a, M=m, Y=y, covariates = x, fdr.method = fdr.method, verbose = TRUE)

    js.v <- apply(m, 1, FUN=function(m.v){
      joint_signif_mediation(E=a, M=m.v, Y=y, covariates = x)
    })
    js.res <- data.frame(p=js.v, FDR=stats::p.adjust(js.v, method=fdr.method))

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

    if (n.inconsistent > 0){
      prop.sig.arr["js_both", "n.hits", sim] <- sum(js.res[, "FDR"] < FDR)
      prop.sig.arr["js_both", "power", sim] <- mean(js.res[med_genes, "FDR"] < FDR)
      prop.sig.arr["js_both", "false_disc", sim] <-
        sum(js.res[setdiff(g.nms, med_genes), "FDR"] < FDR)
    }

    if (verbose && sim %% 100 == 0) message("sim: ", sim)
  }
  fdr.mat <- prop.sig.arr[,"false_disc",, drop=TRUE]/prop.sig.arr[,"n.hits",, drop=TRUE]
  fdr.mat[is.na(fdr.mat)] <- 0
  prop.sig.arr[, "fdr",] <- fdr.mat
  ret <- data.frame(apply(prop.sig.arr, 1:2, FUN=mean))
  ret
}
