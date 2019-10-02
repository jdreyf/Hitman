#' ScreenMin multiple hypthesis testing for union hypothesis tests
#'
#' ScreenMin multiple hypthesis testing for union hypothesis tests from Djordjilovic et al. (2018).
#'
#' @param p.tab Numberic matrix with two columns of p-values.
#' Each row is an analyte. No \code{NA}'s allowed.
#' @param thresh Numeric threshold, termed \code{c} in Djordjilovic et al. (2019). Min p-values > \code{thresh}
#' are assigned adjusted p-values of one.
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams ezlimma::ezcor
#' @references Djordjilovic V, Page CM, Gran JM, Nost TH, Sandanger TM, Veierød MB, Thoresen M.
#' Global test for high-dimensional mediation: Testing groups of potential mediators. Stat Med. 2019 Aug 15;38(18):3346-3360.

screen_min <- function(p.tab, thresh=0.05, reorder.rows = TRUE, prefix = NULL, adjust.method = "BH"){
  stopifnot(ncol(p.tab) == 2, is.matrix(p.tab), is.numeric(p.tab), is.finite(p.tab), p.tab >= 0, p.tab <= 1,
            !is.null(rownames(p.tab)))
  p.tab.o <- t(apply(p.tab, MARGIN = 1, FUN=sort))
  colnames(p.tab.o) <- c("minp", "maxp")
  qv <- rep(1, nrow(p.tab.o))

  S.ind <- which(p.tab.o[, "minp"] <= thresh)
  if (length(S.ind) >= 1) qv[S.ind] <- stats::p.adjust(p.tab.o[S.ind, "maxp"], method = adjust.method)
  # modified from limma::eztoptab
  adjp.colnm <- ifelse(adjust.method %in% c("BH", "fdr"), yes="FDR", no=adjust.method)
  ret <- cbind(p.tab.o,  qv)
  dimnames(ret) <- list(rownames(p.tab.o), c(colnames(p.tab.o), adjp.colnm))
  if (reorder.rows) ret <- ret[order(ret[, "minp"] > thresh, ret[, "maxp"]), ]
  if (!is.null(prefix)) colnames(ret) <- paste(prefix, colnames(ret), sep=".")
  ret
}
