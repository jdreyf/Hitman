#' Joint significant mediation test
#'
#' Joint significant mediation test.
#'
#' @param M Numeric vector of mediators with one element per sample.
#' @inherit hitman
#' @inheritParams ezlimma::ezcor
#' @return Data frame with columns \code{EMY.p}, \code{EM.p}, and \code{MY.p}.
#' @seealso \url{https://github.com/cedricbatailler/JSmediation}.

# p-value is max of both; since both need to be significant for test to be significant
joint_signif_mediation <- function(E, M, Y, covariates=NULL){
  stopifnot(is.numeric(E), is.numeric(M), is.numeric(Y), !is.na(E), !is.na(M), !is.na(Y), length(E) > 0,
            length(Y)==length(E), length(Y)==length(M))
  if (is.null(covariates)){
    dat <- data.frame(E, M, Y)
    em.fm <- summary(stats::lm(M~E, data=dat))
    my.fm <- summary(stats::lm(Y~M+E, data=dat))
  } else {
    stopifnot(is.numeric(covariates), length(covariates)==length(E) || nrow(covariates)==length(E))
    dat <- data.frame(E, M, Y, covariates)
    em.fm <- summary(stats::lm(M~E+covariates, data=dat))
    my.fm <- summary(stats::lm(Y~M+E+covariates, data=dat))
  }
  em.p <- em.fm$coefficients["E", "Pr(>|t|)"]
  my.p <- my.fm$coefficients["M", "Pr(>|t|)"]
  # ret.p <- matrix(max(em.p, my.p), nrow=1, ncol=1, dimnames = list("row1", "EMY.p"))
  ret.p <- matrix(c(max(em.p, my.p), em.p, my.p), nrow=1, dimnames=list("row1", c("EMY.p", "EM.p", "MY.p")))
  return(ret.p)
}
