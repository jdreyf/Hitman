#' High-throughput approach for testing directional replication that uses filtering to improve adjusted p-values
#'
#' High-throughput approach for testing replication of two base studies, where each study applies a two-sided test whose
#' signed statistic and p-value are supplied via `tab`, and its desired to find rows where there is replication in a common direction.
#' Hitman2 improves on the adjusted p-values using filtering.
#'
#' @param p.adj.rate Either "FDR" for false discovery rate or "FWER" for family-wise error rate, the rate controlled by the Bonferroni procedure.
#' @inheritParams hitman_replication
#' @inheritParams hitman
#' @inheritParams ezlimma::limma_cor
#' @return Data frame whose rows correspond to the rows of `tab` with the same row names and whose columns are
#' \describe{
#' \item{chisq}{Chi-square for replication on 1 degreee of freedom.}
#' \item{p}{P-value for replication}
#' \item{FDR or FWER}{FDR or FWER for replication}
#' }
#' @details Larger chi-square values are more significant.
#' @md

hitman2_replication <- function(tab, cols=1:4, reorder.rows=FALSE, p.adj.rate=c("FDR", "FWER"), prefix=NA){
  p.adj.rate <- match.arg(p.adj.rate, c("FDR", "FWER"))
  stopifnot(nrow(tab) > 0, cols %in% c(1:ncol(tab), colnames(tab)), length(cols)==4, limma::isNumeric(tab))
  stat.cols <- cols[c(1, 3)]
  p.cols <- cols[c(2, 4)]
  # don't reorder rows if don't have row names to identify the rows
  stopifnot(0 <= tab[, p.cols], tab[, p.cols] <= 1, !is.na(tab), is.logical(reorder.rows), !is.null(prefix), !is.null(rownames(tab)))

  hm <- hitman_replication(tab=tab, cols=cols, reorder.rows=FALSE, fdr.method="BH", prefix=NA) |>
    dplyr::select(!FDR)
  M <- nrow(tab)

  tab2 <- data.frame(tab[, p.cols], hm) |>
    dplyr::mutate(minp = apply(as.matrix(tab[, p.cols]), MARGIN=1, FUN=min), max2 = pmax(p, minp)) |>
    dplyr::arrange(max2) |>
    dplyr::mutate(rnk = 1:M)
  tab2$adj.num <- apply(as.matrix(tab2$max2), 1, FUN=function(xx) sum(tab2$minp <= xx))

  if (p.adj.rate == "FDR"){
    tab2 <- tab2 |> dplyr::mutate(bh.point = adj.num*max2/rnk, FDR = cummin(bh.point[M:1])[M:1]) |>
      dplyr::select(!(minp:bh.point))
  } else {
    tab2 <- tab2 |> dplyr::mutate(holm.point = adj.num*max2, FWER = pmin(cummin(holm.point[M:1])[M:1], 1)) |>
      dplyr::select(!(minp:holm.point))
  }

  if (reorder.rows){
    tab2 <- tab2 |> dplyr::arrange("p")
  } else {
    tab2 <- tab2[rownames(tab),]
  }
  if (!is.na(prefix)) colnames(tab2) <- paste(prefix, colnames(tab2), sep=".")
  tab2 |> dplyr::select(-(1:2))
}
