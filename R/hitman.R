#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}. Before applying Hitman, you should know a priori the direction that \code{E} changes \code{Y} and verify
#' that \code{E} changes \code{Y} significantly and in the same direction here. See examples in vignette.
#'
#' @param E A numeric vector of exposures.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample of mediators.
#' Must have more than one feature and have row names and column names.
#' @param Y A numeric vector of \code{length(E)} of outcomes. Only continuous, normally distributed outcomes
#' currently supported.
#' @param covariates Numeric vector with one element per sample or matrix-like object with rows corresponding
#' to samples and columns to covariates to be adjusted for.
#' @param fam Character string of family to use in generalized linear model of \code{Y}. The default
#' \code{"gaussian"} reduces to the usual linear regression model. See stats::family.
#' @param fdr.method Character string; either "BH" for Benjamini-Hochberg or "BY" for Benjamini-Yekutieli.
#' See stats::p.adjust.
#' @param verbose Logical; should messages be given for lack of association between \code{E} & \code{Y} and filtering?
#' @param check.names Logical; should \code{names(E)==colnames(M) & colnames(M)==names(Y)} be checked?
#' @inheritParams ezlimma::ezcor
#' @return Data frame with columns
#' \describe{
#' \item{EMY.chisq}{Overall chi-square for mediation on 1 degreee of freedom.}
#' \item{EMY.p}{Overall p-value for mediation}
#' \item{EMY.FDR}{Overall FDR for mediation}
#' \item{EM.z}{z-score for E-->M, not accounting for direction}
#' \item{EM.p}{p-value for E-->M, not accounting for direction}
#' \item{MY.z}{z-score for M-->Y, not accounting for direction}
#' \item{MY.p}{p-value for M-->Y, not accounting for direction}
#' }
#' @details \code{E} and \code{Y} cannot have \code{NA}s. \code{M} may have some \code{NA}s, but rows that have
#' less non-missing values than \code{5 + ncol(covariates)} will be filtered out, and if \code{verbose=TRUE},
#' a message will be written with the number of rows filtered out.
#'
#' Larger chi-square values are more significant.
#' @export

hitman <- function(E, M, Y, covariates=NULL, fam= "gaussian", reorder.rows=TRUE, fdr.method=c("BH", "BY"),
                   verbose=TRUE, check.names=TRUE){

  fdr.method <- match.arg(fdr.method, c("BH", "BY"))
  stopifnot(is.numeric(E), is.numeric(Y), !is.na(E), !is.na(Y), is.null(dim(E)), is.null(dim(Y)), stats::var(E) > 0,
            nrow(M) > 1, length(E)==ncol(M), length(Y)==ncol(M), !is.null(rownames(M)), !is.null(colnames(M)),
            length(unique(Y)) >= 3 || fam == "binomial", limma::isNumeric(M) || class(M)=="EList")
  if (check.names) stopifnot(names(E)==colnames(M), colnames(M)==names(Y))

  # ok if covariates is NULL
  my.covar <- cbind(E=E, covariates=covariates)

  # filter
  # regression df = N-k-1, then lose 2 for removeBatchEffect, so N-k-3 >= 1, i.e. need N >= k+4
  min.df <- ncol(my.covar) + 4
  keep.rows <- apply(M, MARGIN=1, FUN=function(vv){
    sum(!is.na(vv)) >= min.df
  })
  if (sum(keep.rows) == 0){
    stop("No rows had ", min.df, " non-missing values.")
  } else if (sum(keep.rows) < nrow(M)){
    if (verbose) message(nrow(M) - sum(keep.rows), " row(s) were filtered for not having ", min.df, " non-missing values.")
    M <- M[keep.rows,, drop=FALSE]
  }
  # if all rows should be kept, do nothing

  # test EY; return ey.sign & weak assoc warning
  if (fam == "gaussian"){
    fm.ey <- stats::lm(Y ~ ., data=data.frame(Y, my.covar))
    tt.ey <- c(EY.t=summary(fm.ey)$coefficients["E", "t value"], EY.p=summary(fm.ey)$coefficients["E", "Pr(>|t|)"])
  } else {
    fm.ey <- stats::glm(Y ~ ., data=data.frame(Y, my.covar), family=fam)
    tt.ey <- c(EY.t=summary(fm.ey)$coefficients["E", "z value"], EY.p=summary(fm.ey)$coefficients["E", "Pr(>|z|)"])
  }

  if (tt.ey["EY.p"] > 0.05 && verbose){
    message("E and Y are not associated at p<0.05 (p=", signif(tt.ey["EY.p"], digits = 2),
            "), so mediation may not be meaningful.")
  }
  ey.sign <- sign(tt.ey["EY.t"])

  # change order of columns so it's consistent with c("MY.p", "MY.slope")
  # include intercept in the design matrix
  design <- stats::model.matrix(~., data=data.frame(my.covar))
  tt.em <- ezlimma::limma_cor(object=M, design=design, coef=2, prefix="EM", cols=c("z", "P.Value"))

  # don't need to recheck names
  tt.my <- limma_pcor(object=M, phenotype=Y, covariates=my.covar, prefix="MY", check.names=FALSE, cols=c("z", "P.Value"))
  tt.my <- tt.my[,setdiff(colnames(tt.my), "MY.FDR")]
  ret <- cbind(tt.em[rownames(tt.my),], tt.my)

  # modify separate columns, to keep stats of two-sided tests for inspection.
  # ret <- cbind(EM_dir.p=ret$EM.p, MY_dir.p=ret$MY.p, ret)
  # # p.cols <- c("EM_dir.p", "MY_dir.p")
  # ret <- modify_hitman_pvalues(tab=ret, overall.sign = ey.sign, p.cols=p.cols)
  # p.cols=c("EM.p", "MY.p")

  stat.cols=c("EM.z", "MY.z")
  p.cols=c("EM.p", "MY.p")
  # use chi-sq so p=1 --> chisq=0
  EMY.chisq <- EMY.p <- rep(NA, nrow(ret))
  sgn <- apply(ret[, stat.cols], MARGIN=1, FUN=function(vv) sign(prod(vv)))
  eq.sgn <- sgn == ey.sign
  neq.sgn <- !eq.sgn
  # order columns
  p.tab.o <- t(apply(data.matrix(ret[, p.cols, drop=FALSE]), MARGIN = 1, FUN=sort, na.last=TRUE))
  colnames(p.tab.o) <- c("minp", "maxp")

  if (any(neq.sgn)){
    EMY.p[which(neq.sgn)] <- 1
  }

  if (any(eq.sgn)){
    EMY.p[which(eq.sgn)] <- 0.5*p.tab.o[which(eq.sgn), "maxp"]
  }

  EMY.chisq <- stats::qchisq(p=EMY.p, df=1, lower.tail = FALSE)
  EMY.FDR <- stats::p.adjust(EMY.p, method = fdr.method)

  ret <- cbind(EMY.chisq, EMY.p, EMY.FDR, ret)
  if (reorder.rows) ret <- ret[order(ret$EMY.p),]
  return(ret)
}
