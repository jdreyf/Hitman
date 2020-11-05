#' Test partial correlation of each row of an object to a phenotype vector
#'
#' Test partial correlation of each row of an object to a phenotype vector given covariates.
#'
#' @inheritParams ezlimma::limma_cor
#' @inheritParams hitman
#' @inheritParams ezlimma::limma_contrasts
#' @details \code{covariates} should not include the regression intercept.
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
  # removeBatchEffects doesn't have an intercept
  object.res <- limma::removeBatchEffect(x=object, covariates = covariates)

  # how many df to remove in limma_cor, inc. intercept
  reduce.df <- ncol(as.matrix(covariates)) + 1

  # need intercept b/c object.res not centered at 0
  lc <- ezlimma::limma_cor(object=object.res, phenotype = pheno.res, reduce.df=reduce.df, reorder.rows=reorder.rows,
                           prefix=prefix, adjust.method=adjust.method, cols=cols)
  return(lc)
}
