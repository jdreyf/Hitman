#' Simulate & test mediation methods on omics dataset
#'
#' Simulate & test mediation methods on omics dataset.
#'
#' @param b1 Numeric value that beta1 (\code{"b1"}), the EM coefficient, takes for non-null mediators.
#' @param t2 Numeric value that theta2 (\code{"t2"}), the MY coefficient, takes for non-null mediators.
#' @param t1 Numeric value of theta1 i.e. the effect of the exposure on the outcome.
#' @param nsamp Number of samples.
#' @param ngene Number of genes.
#' @param FDR FDR threshold to apply.
#' @param Sigma Gene covariance matrix.
#' @param prop.consistent Proportion of genes that are consistent mediators. Must be at least \code{1/ngene}.
#' @param prop.inconsistent Proportion of genes that are inconsistent mediators.
#' @param prop.em Proportion of genes where E --> M but not M --> Y given E.
#' @param prop.my Proportion of genes where M --> Y given E but not E --> M.
#' @param sd.mn Numeric standard deviation of measurement noise.
#' @inheritParams ezlimma:::sim_fisher
#' @inheritParams hitman
#' @return Matrix with proportion of significant calls for every method for true and null mediators.
#' @export

sim_omics <- function(b1=1, t2=b1, t1=5, nsamp=50, ngene=100, FDR=0.25, Sigma=diag(ngene), prop.consistent=1/ngene,
                      prop.inconsistent=0, prop.em=0, prop.my=0, sd.mn=0, nsim=10**3, fdr.method=c("BH", "BY"), seed=0,
                      verbose=TRUE){

  fdr.method <- match.arg(fdr.method, c("BH", "BY"))

  stopifnot(prop.consistent >= 1/ngene, prop.inconsistent >= 0, prop.consistent <= 1, prop.inconsistent <= 1,
            prop.consistent + prop.consistent + prop.em + prop.my <= 1, prop.em >= 0, prop.my >= 0)

  set.seed(seed)
  #t = theta; b = beta
  t0 <- b0 <- b2 <- t3 <- 0.14

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
    one_comp_genes_em <- sample(x=setdiff(g.nms, med_genes), size=floor(prop.em*ngene))
    one_comp_genes_my <- sample(setdiff(g.nms, union(med_genes, one_comp_genes_em)), size=floor(prop.my*ngene))

    # a is exposure
    a <- stats::rnorm(n=nsamp)
    # x is covariate
    x <- stats::rnorm(n=nsamp)

    # signal variance in M
    m_signal_var <- t(MASS::mvrnorm(n=nsamp, mu = rep(0, ngene), Sigma = Sigma))
    m_signal <- b0 + m_signal_var + matrix(b2*x, nrow=ngene, ncol=nsamp, byrow = TRUE)
    dimnames(m_signal) <- dimnames(m_signal_var) <- list(g.nms, paste0("s", 1:nsamp))
    # contribution from exposure
    m_signal[med_genes,] <- m_signal[med_genes,, drop=FALSE] + matrix(b1*a, nrow=length(med_genes), ncol=nsamp, byrow = TRUE)
    if (length(one_comp_genes_em) > 0){
      em.sgns <- sample(x=c(-1, 1), size=length(one_comp_genes_em), replace = TRUE, prob=rep(0.5, 2))
      if (any(em.sgns == -1)){
        g.tmp <- one_comp_genes_em[em.sgns == -1]
        m_signal[g.tmp,] <- m_signal[g.tmp,, drop=FALSE] - matrix(b1*a, nrow=length(g.tmp), ncol=nsamp, byrow = TRUE)
      }
      if (any(em.sgns == 1)){
        g.tmp <- one_comp_genes_em[em.sgns == 1]
        m_signal[g.tmp,] <- m_signal[g.tmp,, drop=FALSE] - matrix(b1*a, nrow=length(g.tmp), ncol=nsamp, byrow = TRUE)
      }
    }
    m <- m_signal + t(MASS::mvrnorm(n=nsamp, mu = rep(0, ngene), Sigma=sd.mn*diag(ngene)))

    # outcome: modified by mediator genes
    y <- t0 + t1*a + t3*x + stats::rnorm(n=nsamp) + colSums(t2 * m_signal[consistent_genes,, drop=FALSE])

    if (length(one_comp_genes_my) > 0){
      my.sgns <- sample(x=c(-1, 1), size=length(one_comp_genes_my), replace = TRUE, prob=rep(0.5, 2))
      if (any(my.sgns == -1)){
        g.tmp <- one_comp_genes_my[my.sgns == -1]
        y <- y - colSums(t2 * m[g.tmp,, drop=FALSE])
      }
      if (any(my.sgns == 1)){
        g.tmp <- one_comp_genes_my[my.sgns == 1]
        y <- y + colSums(t2 * m_signal[g.tmp,, drop=FALSE])
      }
    }

    if (n.inconsistent > 0) y <- y - colSums(t2 * m[inconsistent_genes,, drop=FALSE])

    names(a) <- names(x) <- names(y) <- paste0("s", 1:nsamp)

    hm.res <- hitman(E=a, M=m, Y=y, covariates = x, fdr.method = fdr.method, verbose=TRUE)
    lm.res <- lotman(E=a, M=m, Y=y, covariates = x, fdr.method = fdr.method, verbose = TRUE)

    js.v <- apply(m, MARGIN=1, FUN=function(m.v){
      joint_signif_mediation(E=a, M=m.v, Y=y, covariates = x)[1, "EMY.p"]
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
  } # end for sim
  fdr.mat <- prop.sig.arr[,"false_disc",, drop=TRUE]/prop.sig.arr[,"n.hits",, drop=TRUE]
  fdr.mat[is.na(fdr.mat)] <- 0
  prop.sig.arr[, "fdr",] <- fdr.mat
  ret <- data.frame(apply(prop.sig.arr, 1:2, FUN=mean))
  ret
}
