context("sim-omics")

test_that("compare sims sc1", {
  prop.sig.hit <- Hitman::sim_barfield(med.fnm = "hitman", b1=0.39, nsim = 10, ngene=9, nsamp = 50, t1=0.59)
  # FDR=0.5 if no other genes affect top FDR; otherwise, FDR will be lower, so have more power
  ngene <- 10
  sc1 <- Hitman::sim_omics(nsamp=50, ngene=ngene, FDR=0.5, prop.consistent=1/ngene, b1 = 0.39, t1=0.59, nsim=10,
                            Sigma = diag(ngene), prop.1c = 0, seed=0)
  expect_gte(sc1["hitman", "power"], prop.sig.hit)

  ngene <- 20

  # t1 < 0 --> mediators should not be significant
  sc1.negt1 <- Hitman::sim_omics(nsamp=50, ngene=ngene, FDR=0.5, prop.consistent=1/ngene, b1 = 0.39, t1 = -5, nsim=10,
                           Sigma = diag(ngene), prop.1c = 0, seed=0)
  expect_lte(sc1.negt1["hitman", "power"], 0)
  expect_lte(sc1.negt1["lotman", "power"], 0)
  expect_gt(sc1.negt1["js", "power"], 0)

  # FDR effect
  sc1.fdr15 <- Hitman::sim_omics(nsamp=50, ngene=ngene, FDR=0.15, prop.consistent=1/ngene, b1 = 0.39, t1=0.59, nsim=10,
                                      Sigma = diag(ngene), prop.1c = 0, seed=0)
  expect_true(all(sc1.fdr15$n.hits <= sc1$n.hits))
  expect_true(all(sc1.fdr15$power <= sc1$power))

  # b1
  sc1.lob1 <- Hitman::sim_omics(nsamp=50, ngene=ngene, FDR=0.5, prop.consistent=1/ngene, b1 = 0.39, t1=2, nsim=10,
                           Sigma = diag(ngene), prop.1c = 0, seed=0)
  sc1.hib1 <- Hitman::sim_omics(nsamp=50, ngene=ngene, FDR=0.5, prop.consistent=1/ngene, b1 = 1, t1=2, nsim=10,
                                  Sigma = diag(ngene), prop.1c = 0, seed=0)
  expect_true(all(sc1.hib1$power >= sc1.lob1$power))

  # prop but don't know effect
})
