context("sim-omics")

test_that("sims match", {
  prop.sig.hit <- Hitman::sim_barfield(med.fnm = "hitman", b1t2=0.39, nsim = 10, ngene=9, nsamp = 50, t1=0.59)
  # FDR=0.5 if no other genes affect top FDR; otherwise, FDR will be lower
  sc1 <- Hitman::sim_omics(nsamp=50, ngene=10, FDR=0.5, prop.consistent=1/10, b1t2 = 0.39, t1=0.59, nsim=10,
                            Sigma = diag(10), prop.1c = 0, seed=0)
  expect_gte(sc1["hitman", "power"], prop.sig.hit)
})
