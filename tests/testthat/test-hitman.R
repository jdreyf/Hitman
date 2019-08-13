context("hitman")

test_that("E numeric", {
  hm <- hitman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)

  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)

  #no variance
  expect_error(hitman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))

  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_message(hm3 <- hitman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("one-sided for correct sign", {
  hm <- hitman(E=ee, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)

  hm2 <- hitman(E=ee, M=M, Y=pheno.v, covariates=covar.tmp)
  expect_lte(mean(hm$EMY.p==hm2[rownames(hm), "EMY.p"]), 0.01)

  #no variance
  expect_error(hitman(E=numeric(length(pheno.v)), M=M, Y=pheno.v))

  ee2 <- rnorm(length(pheno.v), sd=0.1)
  expect_message(hm3 <- hitman(E=ee2, M=M, Y=pheno.v))
  expect_lt(mean(hm3$EMY.p < 0.05), 0.1)
})

test_that("E binary", {
  hm <- hitman(E=grp2, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)

  #EMY.p indep of parametrization
  grp3 <- as.numeric(grp==grp[1])
  hm2 <- hitman(E=grp3, M=M, Y=pheno.v)
  expect_equal(hm$EMY.p, hm2$EMY.p)

  expect_message(hm3 <- hitman(E=grp2, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)

  y <- rep(1:3, times=2)
  expect_message(hm4 <- hitman(E=grp2, M=M, Y=rep(1:3, times=2)))
})

test_that("E nominal --> design", {
  grp.tmp <- ezlimma:::batch2design(rep(letters[1:2], each=3))[,1]
  names(grp.tmp) <- colnames(M)

  hm <- hitman(E=grp.tmp, M=M, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.2)

  set.seed(0)
  covar.tmp <- rnorm(length(pheno.v))
  # warning: essentially perfect fit: summary may be unreliable
  expect_warning(hm3 <- hitman(E=grp.tmp, M=M, Y=pheno.v, covariates=covar.tmp))
  expect_lte(mean(hm$EMY.p == hm3[rownames(hm), "EMY.p"]), 0.01)

  expect_error(hitman(E=rep("a", length(pheno.v)), M=M, Y=pheno.v))
  expect_error(hitman(E=c(rep("a", length(pheno.v)-1), NA), M=M, Y=pheno.v))
  expect_error(hitman(E=rep(NA, length(pheno.v)-1), M=M, Y=pheno.v))
})

test_that("M df", {
  mdf <- data.frame(M)
  hm <- hitman(E=ee, M=mdf, Y=pheno.v)
  expect_lt(mean(hm$EMY.p < 0.05), 0.1)
})

test_that("gene1", {
  hm <- hitman(E=grp2, M=M, Y=M[1,])
  expect_equal(rownames(hm)[1], "gene1")

  expect_lt(hm["gene1", "MY_dir.p"], hm["gene1", "MY.p"])
  expect_equal(hm["gene1", "EM.p"], hm["gene1", "EM_dir.p"])
  # need to substract for rounding of both values
  expect_lte(hm["gene1", "EMY.p"], max(hm["gene1", "EM.p"], hm["gene1", "MY.p"]))
})

test_that("NAs", {
  expect_error(hitman(E=grp2, M=M, Y=pheno2))

  grp2[1] <- NA
  expect_error(hitman(E=grp2, M=M, Y=pheno.v))

  M[1, 1] <- NA
  expect_silent(hitman(E=ee, M=M, Y=pheno.v))
})

test_that("consistent & inconsistent", {
  n <- 10
  sigma <- 0.25
  E <- rep(0:1, each=n)
  ey <- rnorm(n=2*n, sd=sigma)
  Y <- E+ey
  eps <- rnorm(n=2*n, sd=sigma)
  em.ic <- -ey+eps
  em.c <- ey+eps
  M <- rbind(ics=E+em.ic, cs=E+em.c)
  hm <- hitman(E=E, M=M, Y=Y)
  expect_lt(hm["cs", "EMY.p"], 0.01)
  expect_gt(hm["ics", "EMY.p"], 0.9)
})

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.mat <- sim_barfield(med.fnm = "hitman", b1t2.v=c(0, 0.39), nsim = 50, ngene = 9, verbose = FALSE)
  expect_lte(prop.sig.mat[1, 1], 0.05)
  expect_gte(prop.sig.mat[2, 2], 0.1)
})
