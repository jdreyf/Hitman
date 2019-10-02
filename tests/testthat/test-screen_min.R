context("screen min")

test_that("screen min proof", {
  set.seed(0)
  thresh <- 0.05
  p_tab0 <- matrix(runif(n=10**4), ncol = 2)
  rownames(p_tab0) <- paste0("g", 1:nrow(p_tab0))
  p_tab <- t(apply(p_tab0, MARGIN=1, FUN=sort))
  colnames(p_tab) <- c("minp", "maxp")

  # check proof that maxp|minp<thresh ~ U(0, 1)
  pp <- p_tab[p_tab[, "minp"] <= thresh, "maxp"]
  ks <- ks.test(x=pp, y="punif", alternative = "greater")
  expect_gte(ks$p.value, 0.2)

  # screen min reproduces above
  sc.min <- screen_min(p.tab=p_tab0, thresh = thresh)
  expect_true(!is.null(rownames(sc.min)))
  expect_true(!is.null(colnames(sc.min)))
  # ordered by minp > thresh
  expect_true(all(diff(sc.min[, "minp"] > thresh) >= 0))
  expect_true(all(sc.min[, 1] <= sc.min[, 2]))

  ks2 <- ks.test(x=sc.min[sc.min[, "minp"] <= thresh, "maxp"], y="punif", alternative = "greater")
  expect_equal(ks$p.value, ks2$p.value)
})
