#' High-throughput approach for testing directional replication
#'
#' High-throughput approach for testing replication of two base studies, where each study applies a two-sided test whose
#' signed statistic and p-value are supplied via `tab`, and its desired to find rows where there is replication in a common direction.
#'
#' @param tab Matrix-like object with statistical and p-value columns. Only the signs of the statistics columns are used.
#' `tab` should have non-duplicated row names and should not have missing values.
#' @param cols Vector of column indices or names in the order of `c(stat1, p1, stat2, p2)`.
#' @param prefix Character string of length one with prefix of returned columns, e.g. if `prefix="repl"`, returned columns might
#' be `c("repl.chisq", "repl.p", "repl.FDR")`.
#' Prefix is not added if it is `NA`.
#' @inheritParams hitman
#' @inheritParams ezlimma::limma_cor
#' @return Data frame whose rows correspond to the rows of `tab` with the same row names and whose columns are
#' \describe{
#' \item{chisq}{Chi-square for replication on 1 degreee of freedom.}
#' \item{p}{P-value for replication}
#' \item{FDR}{FDR for replication}
#' }
#' @details Larger chi-square values are more significant.
#' @md
#' @export

hitman_replication <- function(tab, cols=1:4, reorder.rows=FALSE, fdr.method=c("BH", "BY"), prefix=NA){
  fdr.method <- match.arg(fdr.method, c("BH", "BY"))
  stopifnot(nrow(tab) > 0, cols %in% c(1:ncol(tab), colnames(tab)), length(cols)==4, limma::isNumeric(tab))
  stat.cols <- cols[c(1, 3)]
  p.cols <- cols[c(2, 4)]
  # require rownames for consistency with hitman2_replication, which needs them st can reorder
  stopifnot(0 <= tab[, p.cols], tab[, p.cols] <= 1, !is.na(tab), is.logical(reorder.rows), !is.null(prefix), !is.null(rownames(tab)))

  if (all(tab[, stat.cols] >= 0) | all(tab[, stat.cols] <= 0)){
    warning("All stats are the same sign, which is possible but unlikely for two-sided stats.")
  }

  res <- matrix(NA, nrow=nrow(tab), ncol=3, dimnames=list(rownames(tab), c("chisq", "p", "FDR")))
  sgn <- apply(tab[, stat.cols], MARGIN=1, FUN=function(vv) sign(prod(vv)))
  # order columns per row
  p.tab.o <- t(apply(data.matrix(tab[, p.cols]), MARGIN = 1, FUN=sort, na.last=TRUE))
  colnames(p.tab.o) <- c("minp", "maxp")

  if (any(sgn < 1)) res[which(sgn < 1), "p"] <- 1
  if (any(sgn == 1)) res[which(sgn == 1), "p"] <- 0.5*p.tab.o[which(sgn == 1), "maxp"]

  res[, "chisq"] <- stats::qchisq(p=res[, "p"], df=1, lower.tail = FALSE)
  res[, "FDR"] <- stats::p.adjust(res[, "p"], method = fdr.method)

  if (reorder.rows) res <- res[order(res[, "p"]),]
  if (!is.na(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  data.frame(res)
}
