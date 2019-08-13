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

  gr.rows <- list(p1=which(tab[, p.cols[1]] > tab[, p.cols[2]]))
  gr.rows$p2 <- setdiff(1:nrow(tab), gr.rows$p1)

  for (gr.ind in seq_along(gr.rows)){
    rows.tmp <- gr.rows[[gr.ind]]
    if (length(rows.tmp) > 0){
      # look at other sign * overall sign --> alternative
      gr.other.ind <- setdiff(1:2, gr.ind)
      alt.v <- sign(tab[rows.tmp, stat.cols[gr.other.ind]] * overall.sign)

      for (alt.sgn in c(-1, 1)){
        if (any(alt.v == alt.sgn)){
          alt.chr <- ifelse(alt.sgn == 1, "greater", "less")
          tab[rows.tmp[alt.v == alt.sgn], p.cols[gr.ind]] <-
            ezlimma:::two2one_tailed(tab=tab[rows.tmp[alt.v == alt.sgn], c(stat.cols[gr.ind], p.cols[gr.ind]), drop=FALSE],
                                     stat.cols=1, p.cols=2, alternative=alt.chr)
        }
      }# for alt.sgn
    }
  }# for gr.ind
  return(tab)
}
