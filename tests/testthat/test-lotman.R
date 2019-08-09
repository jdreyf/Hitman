context("lotman")

test_that("E numeric", {
  hm <- lotman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)

  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)

  #no variance
  expect_error(lotman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))

  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_message(hm3 <- lotman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("E binary", {
  hm <- lotman(E=grp2, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)

  #EMY.p indep of parametrization
  grp3 <- as.numeric(grp==grp[1])
  hm2 <- lotman(E=grp3, M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2$EMY.p)

  covar2 <- rnorm(length(pheno.v))
  expect_message(hm3 <- lotman(E=grp2, M=M[1,, drop=FALSE], Y=pheno.v, covariates=covar2)) # no assoc warning
  expect_true(hm["gene1", "EMY.p"] != hm3[1, "EMY.p"])

  y <- rep(1:3, times=2)
  expect_message(hm4 <- lotman(E=grp2, M=M, Y=rep(1:3, times=2)))
})

test_that("E nominal --> design", {
  grp.tmp <- ezlimma:::batch2design(rep(letters[1:2], each=3))[,1]
  names(grp.tmp) <- colnames(M)

  hm <- lotman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)

  expect_message(hm3 <- lotman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)

  expect_error(lotman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(lotman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(lotman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("M df", {
  mdf <- data.frame(M)
  hm <- lotman(E=ee, M=mdf, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
})

test_that("gene1", {
  # warning: essentially perfect fit: summary may be unreliable
  expect_warning( hm <- lotman(E=grp2, M=M, Y=M[1,]) )
  expect_equal(rownames(hm)[1], "gene1")

  expect_equal(hm["gene1", "MY.p"], hm["gene1", "MY_dir.p"])
  expect_gte(hm["gene1", "EM.p"], hm["gene1", "EM_dir.p"])
  expect_equal(hm["gene1", "EMY.p"], max(hm["gene1", "EM_dir.p"], hm["gene1", "MY_dir.p"]))
})

test_that("NAs", {
  expect_error(lotman(E=grp2, M=M, Y=pheno2))

  grp2[1] <- NA
  expect_error(lotman(E=grp2, M=M, Y=pheno.v))

  M[1, 1] <- NA
  expect_silent(lotman(E=ee, M=M, Y=pheno.v))
})

test_that("testing a gene is independent of other genes", {
  hm1 <- lotman(E=ee, M=M[1,, drop=FALSE], Y=pheno.v, covariates=covar.tmp)
  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_equal(hm1[1, "EMY.p"], hm2["gene1", "EMY.p"])
})
