context("lotman")

test_that("E numeric", {
  hm <- lotman(E=ee, M=M, Y=pheno.v, reorder.rows = FALSE)
  expect_lt(mean(hm[, "EMY.p"] < 0.05), 0.1)

  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp, reorder.rows = FALSE)
  both1 <- rowMeans(cbind(hm[,"EMY.p"], hm2[,"EMY.p"]) == 1) == 1
  expect_lte(mean(hm[!both1, "EMY.p"]==hm2[!both1, "EMY.p"]), 0.01)

  #no variance
  expect_error(lotman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))

  set.seed(0)
  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_message(hm3 <- lotman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3[, "EMY.p"] < 0.05), 0.1)
})

test_that("E binary", {
  hm <- lotman(E=grp2, M=M, Y=pheno.v)
  expect_lt(mean(hm[, "EMY.p"] < 0.05), 0.2)

  #EMY.p indep of parametrization
  grp3 <- as.numeric(grp==grp[1])
  hm2 <- lotman(E=grp3, M=M, Y=pheno.v)
  expect_equal(hm[, "EMY.p"], hm2[, "EMY.p"])

  covar2 <- rnorm(length(pheno.v))
  expect_message(hm3 <- lotman(E=grp2, M=M[1,, drop=FALSE], Y=pheno.v, covariates=covar2)) # no assoc warning
  expect_true(hm["gene1", "EMY.p"] != hm3["EMY.p"])

  y <- rep(1:3, times=2)
  expect_message(hm4 <- lotman(E=grp2, M=M, Y=rep(1:3, times=2)))
})

test_that("E nominal --> design", {
  grp.tmp <- ezlimma:::batch2design(rep(letters[1:2], each=3))[,1]
  names(grp.tmp) <- colnames(M)

  hm <- lotman(E=grp.tmp, M=M, Y=pheno.v, reorder.rows = FALSE)
  expect_lt(mean(hm[, "EMY.p"] < 0.05), 0.2)

  expect_message(hm3 <- lotman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp, reorder.rows = FALSE))
  both1 <- rowMeans(cbind(hm[,"EMY.p"], hm3[,"EMY.p"]) == 1) == 1
  expect_lte(mean(hm[!both1, "EMY.p"] == hm3[!both1, "EMY.p"]), 0.01)

  expect_error(lotman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(lotman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(lotman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("M df", {
  mdf <- data.frame(M)
  hm <- lotman(E=ee, M=mdf, Y=pheno.v)
  expect_lt(mean(hm[, "EMY.p"] < 0.05), 0.1)
})

test_that("gene1", {
  # warning: essentially perfect fit: summary may be unreliable
  expect_warning( hm <- lotman(E=grp2, M=M, Y=M[1,]) )
  expect_equal(rownames(hm)[1], "gene1")
})

test_that("NAs", {
  expect_error(lotman(E=grp2, M=M, Y=pheno2))

  grp2[1] <- NA
  expect_error(lotman(E=grp2, M=M, Y=pheno.v))

  M[1, 1] <- NA
  expect_silent(lotman(E=ee, M=M, Y=pheno.v))
})

test_that("single gene is independent of other genes", {
  hm1 <- lotman(E=ee, M=M[1,, drop=FALSE], Y=pheno.v, covariates=covar.tmp)
  hm2 <- lotman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_equal(signif(hm1["EMY.p"], 3), setNames(signif(hm2["gene1", "EMY.p"], 3), nm="EMY.p"))
  hm3 <- lotman(E=ee, M=M[1,], Y=pheno.v, covariates=covar.tmp)
  expect_true(all(hm1 == hm3))
})

test_that("barfield t2=0.39; b1=0", {
  # MY.t >> 0
  set.seed(1)
  t0 <- t1 <- t3 <- b0 <- b2 <- 0.14
  t2 <- 0.39; b1 <- 0
  nsamp <- 50
  x <- stats::rnorm(n=nsamp)
  a <- stats::rnorm(n=nsamp)
  # eq 1; E(M1)
  em1 <- b0 + b1*a + b2*x
  # normal here ~ log-normal or inv chi sq on counts
  m1 <- stats::rnorm(n=nsamp, mean=em1)
  # eq 3; E(Y)
  ey <- t0 + t1*a + t2*m1 + t3*x
  y <- stats::rnorm(n=nsamp, mean=ey)
  med.res <- lotman(E=a, M=m1, Y=y, covariates = x, verbose = FALSE)
  expect_gt(med.res[1, "EMY.p"], 0.05)
})

test_that("barfield t2=0; b1=0.39", {
  # EM.t >> 0
  set.seed(0)
  t0 <- t1 <- t3 <- b0 <- b2 <- 0.14
  t2 <- 0; b1 <- 0.39
  nsamp <- 10
  x <- stats::rnorm(n=nsamp)
  a <- stats::rnorm(n=nsamp)
  # eq 1; E(M1)
  em1 <- b0 + b1*a + b2*x
  # normal here ~ log-normal or inv chi sq on counts
  m1 <- stats::rnorm(n=nsamp, mean=em1)
  # eq 3; E(Y)
  ey <- t0 + t1*a + t2*m1 + t3*x
  y <- stats::rnorm(n=nsamp, mean=ey)
  med.res <- lotman(E=a, M=m1, Y=y, covariates = x, verbose = FALSE)
  expect_gt(med.res[1, "EMY.p"], 0.05)
})

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.mat <- sim_barfield(med.fnm = "lotman", b1t2.v=c(0, 0.39), nsim = 50, ngene = 0, verbose = FALSE)
  expect_lte(prop.sig.mat[1, 1], 0.05)
  expect_gte(prop.sig.mat[2, 2], 0.1)
})
