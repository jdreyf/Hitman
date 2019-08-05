#' Modify p-values of hitman or lotman
#'
#' Modify p-values of hitman or lotman based on sidedness of both tests.
#'
#' @param tab Numeric matrix. Must have "EM" and "MY" stat and ".p" columns.
#' @param overall.sign Sign of overall effect from exposure to outcome, one of 1 or -1.
#' @param stat.cols Vector of length 2 with column names or indices of signed statistics.
#' @param p.cols Vector of length 2 with column names or indices of p-values.
#' @return Matrix with p-value columns modified.

modify_hitman_pvalues <- function(tab, overall.sign, stat.cols=c("EM.z", "MY.z"), p.cols=c("EM.p", "MY.p")){
  stopifnot(overall.sign %in% c(1, -1), nrow(tab) > 0, stat.cols %in% colnames(tab), p.cols %in% colnames(tab),
            length(stat.cols)==2, length(p.cols)==2)

  p1.larger.rows <- which(tab[, p.cols[1]] > tab[, p.cols[2]])
  if (length(p1.larger.rows) > 0){
    # look at other sign & overall sign --> alternative
    alt.p1 <- sign(tab[p1.larger.rows, stat.cols[2]] * overall.sign)

    if (any(alt.p1 == -1)){
      tab[p1.larger.rows[alt.p1 == -1], p.cols[1]] <- ezlimma:::two2one_tailed(tab=tab[p1.larger.rows[alt.p1 == -1], c(stat.cols[1], p.cols[1])],
                                                                               stat.cols=1, p.cols=2, alternative="less")
    }

    if (any(alt.p1 == 1)){
      tab[p1.larger.rows[alt.p1 == 1], p.cols[1]] <- ezlimma:::two2one_tailed(tab=tab[p1.larger.rows[alt.p1 == 1], c(stat.cols[1], p.cols[1])],
                                                                             stat.cols=1, p.cols=2, alternative="greater")
    }
  }

  # MY has larger
  p2.larger.rows <- which(tab[, p.cols[1]] < tab[, p.cols[2]])
  if (length(p2.larger.rows) > 0){
    # look at other sign & overall sign --> alternative
    alt.p2 <- sign(tab[p2.larger.rows, stat.cols[1]] * overall.sign)

    if (any(alt.p2 == -1)){
      tab[p2.larger.rows[alt.p2 == -1], p.cols[2]] <- ezlimma:::two2one_tailed(tab=tab[p2.larger.rows[alt.p2 == -1], c(stat.cols[2], p.cols[2])],
                                                                               stat.cols=1, p.cols=2, alternative="less")
    }

    if (any(alt.p2 == 1)){
      tab[p2.larger.rows[alt.p2 == 1], p.cols[2]] <- ezlimma:::two2one_tailed(tab=tab[p2.larger.rows[alt.p2 == 1], c(stat.cols[2], p.cols[2])],
                                                                            stat.cols=1, p.cols=2, alternative="greater")
    }
  }

  return(tab)
}
