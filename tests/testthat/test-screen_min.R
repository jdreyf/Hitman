context("screen min")

test_that("screen min proof", {
  set.seed(0)
  thresh <- 0.05
  p_tab <- matrix(runif(n=10**4), ncol = 2)
  p_tab <- t(apply(p_tab, 1, FUN=sort))
  colnames(p_tab) <- c("minp", "maxp")

  # check proof that maxp|minp < thresh ~ U(0, 1)
  pp <- p_tab[p_tab[, "minp"] < thresh, "maxp"]
  ks <- ks.test(x=pp, y="punif", alternative = "greater")
  expect_gte(ks$p.value, 0.2)
})
