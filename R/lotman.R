#' Low-throughput mediation analysis (Lotman)
#'
#' Low-throughput mediation analysis to test if rows of \code{M} mediate the effect of exposure \code{E} on outcome
#' \code{Y}. Before applying Lotman, you should know a priori the direction that \code{E} changes \code{Y} and verify
#' that \code{E} changes \code{Y} significantly and in the same direction here.
#'
#' @param M A numeric matrix-like data object with one row per feature and one column per sample of mediators.
#' @inherit hitman
#' @inheritParams ezlimma::ezcor
#' @export

# need to modify limma_cor, since ezcor does not handle design
lotman <- function(E, M, Y, covariates=NULL, reorder.rows=TRUE, verbose=TRUE, check.names=TRUE){

  if (is.null(dim(M))){
    M <- matrix(M, nrow=1, dimnames=list("analyte", names(M)))
    reorder.rows <- FALSE
  }
  stopifnot(is.numeric(E), limma::isNumeric(M), is.numeric(Y), !is.na(E), !is.na(Y), is.null(dim(E)),
            is.null(dim(Y)), stats::var(E) > 0, length(unique(Y)) >= 3, length(E)==ncol(M), length(Y)==ncol(M))
  if (check.names) stopifnot(names(E)==colnames(M), colnames(M)==names(Y))

  # ok if covariates is NULL
  my.covar <- cbind(E=E, covariates=covariates)

  # test EY; return ey.sign & weak assoc warning
  fm.ey <- stats::lm(Y ~ ., data=data.frame(Y, my.covar))
  tt.ey <- c(EY.t=summary(fm.ey)$coefficients["E", "t value"], EY.p=summary(fm.ey)$coefficients["E", "Pr(>|t|)"])
  if (tt.ey["EY.p"] > 0.1 && verbose){
    message("E and Y are not associated at p<0.1, so mediation may not be meaningful.")
  }
  ey.sign <- sign(tt.ey["EY.t"])

  # change order of columns so it's consistent with c("MY.p", "MY.slope")
  # include intercept in the design matrix
  # des.em <- stats::model.matrix(~., data=data.frame(my.covar))
  tt.em <- t(apply(as.matrix(M), MARGIN = 1, FUN=function(mm){
    fm <- stats::lm(mm ~ ., data=data.frame(my.covar))
    summary(fm)$coefficients["E", c("t value", "Pr(>|t|)")]
  }))
  colnames(tt.em) <- paste0("EM.", sub("t value", "t",
                                       sub("Pr(>|t|)", "p", colnames(tt.em), fixed=TRUE)))

  # don't need to recheck names
  tt.my <- t(apply(as.matrix(M), MARGIN = 1, FUN=function(vv){
    fm <- stats::lm(Y ~ 1 + vv + ., data=data.frame(my.covar))
    summary(fm)$coefficients[2, c("t value", "Pr(>|t|)")]
  }))
  colnames(tt.my) <- paste0("MY.", sub("t value", "t",
                                       sub("Pr(>|t|)", "p", colnames(tt.my), fixed = TRUE)))

  ret <- cbind(tt.em, tt.my)

  stat.cols=c("EM.t", "MY.t")
  p.cols=c("EM.p", "MY.p")
  # use chi-sq so p=1 --> chisq=0
  EMY.chisq <- EMY.p <- rep(NA, nrow(ret))
  sgn <- apply(ret[, stat.cols, drop=FALSE], MARGIN=1, FUN=function(vv) sign(prod(vv)))
  eq.sgn <- sgn == ey.sign
  neq.sgn <- !eq.sgn
  # order columns
  p.tab.o <- t(apply(data.matrix(ret[, p.cols, drop=FALSE]), MARGIN = 1, FUN=sort))
  colnames(p.tab.o) <- c("minp", "maxp")

  if (any(neq.sgn)){
    EMY.p[which(neq.sgn)] <- 1
  }

  if (any(eq.sgn)){
    EMY.p[which(eq.sgn)] <- 0.5*p.tab.o[which(eq.sgn), "maxp"]
  }

  EMY.chisq <- stats::qchisq(p=EMY.p, df=1, lower.tail = FALSE)
  EMY.FDR <- stats::p.adjust(EMY.p, method = "BH")

  ret <- cbind(EMY.chisq, EMY.p, EMY.FDR, ret)
  if (reorder.rows) ret <- ret[order(ret[, "EMY.p", drop=FALSE]),]
  return(ret)
}
