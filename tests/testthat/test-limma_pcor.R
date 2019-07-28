context("ezpcor")

test_that("batch --> covar", {
  res0 <- limma_pcor(object=M, phenotype=pheno.v, covariates = grp2)
  res1 <- limma_pcor(object=M, phenotype=pheno.v, covariates = as.numeric(grp2==1))
  #parametrization changes only AveExpr
  expect_equal(res0$p, res1$p)

  m2 <- rbind(M, phenotype)
  res2 <- limma_pcor(object=m2, phenotype=pheno.minus.one, covariates = covar)
  expect_equal(rownames(res2)[1], "phenotype")
})

test_that("with covar", {
  m2 <- rbind(M, phenotype, covariates = covar)
  res2 <- limma_pcor(object=m2, phenotype=pheno.minus.one, covariates = covar)
  #possible that phenotype wouldn't be first, but it is most likely
  expect_equal(rownames(res2)[1], "phenotype")

  expect_error(res2 <- limma_pcor(object=m2, phenotype=pheno.minus.one, covariates = grp))
})
