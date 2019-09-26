#' ScreenMin multiple hypthesis testing for union hypothesis tests
#'
#' ScreenMin multiple hypthesis testing for union hypothesis tests from Djordjilovic et al. (2018).
#'
#' @param p.tab Table with two columns. Column 1: min of p-values & column 2: max of p-values.
#' Each row is an analyte.
#' @param thresh Numeric threshold, termed \code{c} in Djordjilovic et al. (2019). Min p-values > \code{thresh}
#' are assigned adjusted p-values of one.
#' @inheritParams ezlimma::limma_contrasts
#' @references Djordjilovic V, Page CM, Gran JM, Nost TH, Sandanger TM, Veierød MB, Thoresen M.
#' Global test for high-dimensional mediation: Testing groups of potential mediators. Stat Med. 2019 Aug 15;38(18):3346-3360.

screen_min <- function(p.tab, thresh=0.05, adjust.method = "BH"){
  stopifnot(ncol(p.tab) == 2, nrow(p.tab) > 1, p.tab[,1] <= p.tab[,2])
  qv <- rep(1, nrow(p.tab))
  if (!is.null(rownames(p.tab))) stats::setNames(qv, nm=rownames(p.tab))
  S.ind <- which(p.tab[,1] <= thresh)
  if (length(S.ind) >= 1) qv[S.ind] <- stats::p.adjust(p.tab[,2], method = adjust.method)
  qv
}
