#' Call \code{mediate} function in \pkg{mediation} similarly to \code{hitman}
#' 
#' Call \code{mediate} function in \pkg{mediation} similarly to \code{hitman}.
#' 
#' @param M Vector of a potential mediator.
#' @param sims Number of simulations per analyte to calculate mediation.
#' @inheritParams hitman

# gives p=0!
ezmediate <- function(E, M, Y, covariates = NULL, sims=10**4){
  if (!requireNamespace("mediation", quietly = TRUE)){
    stop("Install 'mediation' package.", call. = FALSE)
  }
  stopifnot(is.null(covariates) || ncol(as.matrix(covariates)) == 1, is.null(dim(M)))
  
  if (is.null(covariates)){
    med.fit <- stats::lm(M ~ E)
    out.fit <- stats::lm(Y ~ M + E)
  } else {
    med.fit <- stats::lm(M ~ E + covariates)
    out.fit <- stats::lm(Y ~ M + E + covariates)
  }
  
  med.out <- mediation::mediate(med.fit, out.fit, treat = "E", mediator = "M", sims = sims)
  
  # use d.avg.p =  average of average causal mediation effects per treatment & control
  # return matrix for consistency with hitman
  res <- matrix(med.out$d.avg.p, nrow=1, ncol=1)
  colnames(res) <- "EMY.p"
  res
}