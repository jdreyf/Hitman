#' High-throughput mediation analysis
#'
#' High-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}. See examples in vignette.
#'
#' @param E A numeric vector of exposures.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample of mediators.
#' Must have more than one feature.
#' @param Y A numeric vector of \code{length(E)} of outcomes.
#' @param covariates Numeric vector with one element per sample or matrix-like object with rows corresponding
#' to samples and columns to covariates to be adjusted for.
#' @param verbose Logical; should messages & warnings about lack of association between \code{E} & \code{Y} be printed?
#' @param check.names Logical; should \code{names(E)==colnames(M) & colnames(M)==names(Y)} be checked?
#' @return Data frame with columns
#' \describe{
#' \item{EMY.p}{Overall p-value for mediation}
#' \item{EMY.FDR}{Overall FDR for mediation}
#' \item{EM_dir.p}{p-value for E-->M accounting for direction of mediation}
#' \item{MY_dir.p}{p-value for M-->Y accounting for direction of mediation}
#' \item{EM.t}{t-statistic for E-->M, not accounting for direction}
#' \item{EM.p}{p-value for E-->M, not accounting for direction}
#' \item{MY.t}{t-statistic for M-->Y, not accounting for direction}
#' \item{MY.p}{p-value for M-->Y, not accounting for direction}
#' }
#' @details \code{E} and \code{Y} cannot have \code{NA}s.
#' @export

hitman <- function(E, M, Y, covariates=NULL, verbose=TRUE, check.names=TRUE){
  stopifnot(is.numeric(E), limma::isNumeric(M), is.numeric(Y), !is.na(E), !is.na(Y), is.null(dim(E)), is.null(dim(Y)),
            stats::var(E) > 0, stats::var(Y) > 0, nrow(M) > 1, length(E)==ncol(M), length(Y)==ncol(M))
  if (check.names){
    stopifnot(names(E)==colnames(M), colnames(M)==names(Y))
  }

  # ok if covariates is NULL
  my.covar <- cbind(E=E, covariates=covariates)

  # test EY; return ey.sign & weak assoc warning
  fm.ey <- stats::lm(Y ~ ., data=data.frame(Y, my.covar))
  tt.ey <- c(EY.t=summary(fm.ey)$coefficients["E", "t value"], EY.p=summary(fm.ey)$coefficients["E", "Pr(>|t|)"])
  if (tt.ey["EY.p"] > 0.1 && verbose){
    warning("E and Y are not associated at p<0.1, so mediation may not be meaningful.")
  }
  ey.sign <- sign(tt.ey["EY.t"])

  # change order of columns so it's consistent with c("MY.p", "MY.slope")
  # include intercept in the design matrix
  design <- stats::model.matrix(~., data=data.frame(my.covar))
  tt.em <- ezlimma::limma_cor(object=M, design=design, coef=2, prefix="EM", cols=c("t", "P.Value"))

  # don't need to recheck names
  tt.my <- limma_pcor(object=M, phenotype=Y, covariates=my.covar, prefix="MY", check.names=FALSE, cols=c("t", "P.Value"))
  tt.my <- tt.my[,setdiff(colnames(tt.my), "MY.FDR")]
  ret <- cbind(tt.em[rownames(tt.my),], tt.my)

  # modify separate columns, to keep stats of two-sided tests for inspection.
  ret <- cbind(EM_dir.p=ret$EM.p, MY_dir.p=ret$MY.p, ret)
  p.cols <- c("EM_dir.p", "MY_dir.p")
  ret <- modify_hitman_pvalues(tab=ret, overall.sign = ey.sign, p.cols=p.cols)

  EMY.p <- apply(ret[,p.cols], MARGIN=1, FUN=function(v){
    max(v)^1.1
  })
  EMY.FDR <- stats::p.adjust(EMY.p, method="BH")
  EMY.z <- stats::qnorm(p=EMY.p, lower.tail = FALSE)
  ret <- cbind(EMY.z, EMY.p, EMY.FDR, ret)
  ret <- ret[order(ret$EMY.p),]
  return(ret)
}
