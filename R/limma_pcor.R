#' Test partial correlation of each row of an object to a phenotype vector
#'
#' Test partial correlation of each row of an object to a phenotype vector given covariates. The
#' approach here is to estimate the residuals of regressing each of \code{object} and \code{phenotype}
#' on \code{cbind(1, covariates)}, and then to test the correlation of the residuals, with the appropriate
#' degrees of freedom.
#'
#' @inheritParams ezlimma::limma_cor
#' @inheritParams hitman
#' @inheritParams ezlimma::limma_contrasts
#' @details \code{covariates} should not include the regression intercept, but when called from \link{hitman},
#' it should include the exposure.
#'
#' Estimating the residuals corresponding to \code{object} is analogous to
#' \link[limma]{removeBatchEffect}, and several lines of code are adapted from that function.
#'
#' @return Data frame.

limma_pcor <- function(object, phenotype, covariates, fam="gaussian", reorder.rows=TRUE, prefix=NULL,
                       adjust.method="BH", check.names=TRUE, cols=c("t", "P.Value")){
  stopifnot(length(phenotype)==ncol(object), limma::isNumeric(covariates), nrow(as.matrix(covariates))==ncol(object))
  if (check.names){ stopifnot(names(phenotype)==colnames(object)) }

  # pheno residuals
  # data.frame(cbind()) can handle NULLs but not factors
  dat <- data.frame(cbind(phenotype, covariates))
  # phenotype isn't used in RHS
  pheno.fm <- stats::glm(formula=phenotype ~ ., data=dat, family = fam)
  pheno.res <- stats::residuals(pheno.fm)
  # names get corrupted
  names(pheno.res) <- names(phenotype)

  # object residuals; removed as per supplemental text
  # i can't use limma::removeBatchEffects with design=NULL, so using some of its code here
  covar.mat <- cbind(int=rep(1, ncol(object)), covariates)
  fit <- limma::lmFit(object, design = covar.mat)
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  object.res <- as.matrix(object) - beta %*% t(covar.mat)

  # how many df to remove in limma_cor? number of covariates in covar.mat
  reduce.df <- ncol(covar.mat)

  # need intercept b/c object.res not centered at 0
  des <- stats::model.matrix(~0+pheno.res)
  lc <- ezlimma::limma_cor(object=object.res, design = des, reduce.df=reduce.df, coef = 1, reorder.rows=reorder.rows,
                           prefix=prefix, adjust.method=adjust.method, cols=cols)
  return(lc)
}
