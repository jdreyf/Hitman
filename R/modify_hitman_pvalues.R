#' Modify p-values of hitman or lotman
#'
#' Modify p-values of hitman or lotman based on sidedness of both tests.
#'
#' @param tab Numeric matrix. Must have "EM" and "MY" stat and ".p" columns.
#' @param overall.sign Sign of overall effect from exposure to outcome, one of 1 or -1.
#' @param stat.cols Vector of length 2 with column names or indices of signed statistics.
#' @param p.cols Vector of length 2 with column names or indices of p-values.
#' @return Matrix with p-value columns modified.

modify_hitman_pvalues <- function(tab, overall.sign, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM.p", "MY.p")){
  stopifnot(overall.sign %in% c(1, -1), nrow(tab) > 0, stat.cols %in% colnames(tab), p.cols %in% colnames(tab),
            length(stat.cols)==2, length(p.cols)==2)

  # account for sign product
  prod.sign <- sign(tab[,stat.cols[1]])*sign(tab[,stat.cols[2]])
  wrong.sign <- which(prod.sign != overall.sign)
  if (length(wrong.sign) > 0){
    tab.ss <- tab[wrong.sign,]

    p1.larger.rows <- which(tab.ss[, p.cols[1]] > tab.ss[, p.cols[2]])
    if (length(p1.larger.rows) > 0){
      p1.pv <- apply(X=tab.ss[p1.larger.rows, c(stat.cols[1], p.cols[1])], MARGIN=1, FUN=function(v){
        # want to flip alternative to enlarge p-value
        alt.tmp <- ifelse (v[1] > 0, "less", "greater")
        ezlimma:::two2one_tailed(tab=t(as.matrix(v)), stat.cols=1, p.cols=2, alternative=alt.tmp)
      })
      # replace
      tab[wrong.sign, p.cols[1]][p1.larger.rows] <- p1.pv
    }
    p2.larger.rows <- setdiff(1:nrow(tab.ss), p1.larger.rows)
    if (length(p2.larger.rows) > 0){
      p2.pv <- apply(X=tab.ss[p2.larger.rows, c(stat.cols[2], p.cols[2])], MARGIN=1, FUN=function(v){
        # want to flip alternative to enlarge p-value
        alt.tmp <- ifelse (v[1] > 0, "less", "greater")
        ezlimma:::two2one_tailed(tab=t(as.matrix(v)), stat.cols=1, p.cols=2, alternative=alt.tmp)
      })
      tab[wrong.sign, p.cols[2]][p2.larger.rows] <- p2.pv
    }
  }# wrong sign
  return(tab)
}
