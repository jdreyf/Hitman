#' Low-throughput mediation analysis (Lotman)
#'
#' Low-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}.
#' @param M A numeric matrix-like data object with one row per feature and one column per sample of mediators.
#' @inherit hitman
#' @export

# need to modify limma_cor, since ezcor does not handle design
lotman <- function(E, M, Y, covariates=NULL, check.names=TRUE){
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
  if (tt.ey["EY.p"] > 0.1){
    warning("E and Y are not associated, so mediation may not be meaningful.")
  }
  ey.sign <- sign(tt.ey["EY.t"])

  # change order of columns so it's consistent with c("MY.p", "MY.slope")
  # include intercept in the design matrix
  des.em <- stats::model.matrix(~., data=data.frame(my.covar))
  tt.em <- ezlimma::limma_cor(object=M, design=des.em, coef=2, prefix="EM", cols=c("t", "P.Value"), moderated=FALSE)

  # don't need to recheck names
  des.my <- stats::model.matrix(~1+Y+my.covar) #test this
  colnames(des.my) <- sub(pattern="^my.covar", "", x=colnames(des.my))
  tt.my <- ezlimma::limma_cor(object=M, design = des.my, prefix="MY", check.names=FALSE,
                      cols=c("t", "P.Value"), moderated=FALSE)
  tt.my <- tt.my[,setdiff(colnames(tt.my), "MY.FDR")]
  ret <- cbind(tt.em[rownames(tt.my),], tt.my)

  # modify separate columns, to keep stats of two-sided tests for inspection.
  ret <- cbind(EM_dir.p=ret$EM.p, MY_dir.p=ret$MY.p, ret)
  p.cols <- c("EM_dir.p", "MY_dir.p")
  ret <- modify_hitman_pvalues(tab=ret, overall.sign = ey.sign, p.cols=p.cols)

  EMY.p <- apply(ret[,p.cols], MARGIN=1, FUN=function(v){
    max(v)^1.2
  })
  EMY.FDR <- stats::p.adjust(EMY.p, method="BH")
  ret <- cbind(EMY.p, EMY.FDR, ret)
  ret <- ret[order(ret$EMY.p),]
  return(ret)
}
