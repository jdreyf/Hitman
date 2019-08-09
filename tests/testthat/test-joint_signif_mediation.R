context("joint_signif_mediation")

# add test with covariates

test_that("compare to JSmediation", {
  set.seed(0)
  n <- 20
  E <- rnorm(n); M <- rnorm(n); Y <- rnorm(n)
  jsp <- joint_signif_mediation(E, M, Y)[1, 1]

  # compare to JSmediation
  # dat <- data.frame(E, M, Y)
  # max p-value from paths b & c is 0.745
  # mdt_simple(data = dat, E, M, Y)

  expect_lte(abs(jsp - 0.745), 0.001)
})

# takes a few sec -- worth it.
test_that("barfield", {
  prop.sig.mat <- sim_barfield(med.fnm = "joint_signif_mediation", b1t2.v=c(0, 0.39), nsim = 50)
  expect_lte(prop.sig.mat[1, 1], 0.03)
  expect_gte(prop.sig.mat[2, 2], max(prop.sig.mat[2, 1], prop.sig.mat[1, 2]))
  expect_lte(prop.sig.mat[2, 2], 0.65)
})
