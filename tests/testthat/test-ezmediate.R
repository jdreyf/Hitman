context("ezmediate")

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.size <- sim_barfield(med.fcn = ezmediate, b1t2.v=0, nsim = 10, sims=30)
  expect_lte(prop.sig.size[1, 1], 0.2)

  prop.sig.size2 <- sim_barfield(med.fcn = ezmediate, b1t2.v=0, nsim = 10, sims=30)
  expect_lte(prop.sig.size[1, 1], 0.2)

  prop.sig.power <- sim_barfield(med.fcn = ezmediate, b1t2.v=0.39, nsim = 10, sims=30)
  expect_gte(prop.sig.power[1, 1], 0.2)
})

test_that("covars", {
  ezm.res <- ezmediate(E=grp.mat[,1], M=M[1,], Y=phenotype, sims = 50)
  set.seed(123)
  ezm.res2 <- ezmediate(E=grp.mat[,1], M=M[1,], Y=phenotype, covariates = rnorm(n=ncol(M)), sims = 50)
  expect_true(ezm.res[1,1] != ezm.res2[1,1])
})
