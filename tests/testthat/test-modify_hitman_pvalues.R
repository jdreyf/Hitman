context("modify_hitman_pvalues")

#gene at bottom of list
#my.p > em.p & wrong sign --> flip my.p sign
test_that("gene88", {
  v <- ret["gene88",] #a matrix
  prod.sign <- sign(v[,"EM.z"])*sign(v[,"MY.z"])
  expect_false(prod.sign == ey.sign)

  expect_gt(v[,"MY.z"], 0)
  expect_gt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"EM.p"], v[,"EM_dir.p"])
  alt.tmp <- ifelse (v[,"EM.z"] > 0, "greater", "less")
  expect_equal(as.numeric(v["MY_dir.p"]),
               as.numeric(ezlimma:::two2one_tailed(v, stat.col = "MY.z", p.col="MY.p", alternative = alt.tmp)))
})

#my.p < em.p, & wrong sign
test_that("gene72", {
  v <- ret["gene72",] #a matrix
  prod.sign <- sign(v[,"EM.z"])*sign(v[,"MY.z"])
  expect_false(prod.sign == ey.sign)

  expect_lt(v[,"MY.z"], 0)
  expect_lt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"MY.p"], v[,"MY_dir.p"])

  alt.tmp <- ifelse (v[,"EM.z"] > 0, "less", "greater")
  expect_equal(as.numeric(v["EM_dir.p"]),
               as.numeric(ezlimma:::two2one_tailed(v, stat.col = "EM.z", p.col="EM.p", alternative = alt.tmp)))
})

test_that("one gene", {
  tt <- matrix(c(0.00578, 0.849, 7.09, 0.00578, -0.215, 0.849), nrow=1,
               dimnames=list("gene1", c('EM_dir.p', 'MY_dir.p', 'EM.t', 'EM.p', 'MY.t', 'MY.p')))
  res1 <- modify_hitman_pvalues(tab=tt, overall.sign = -1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_equal(nrow(res1), 1)
  expect_equal(res1[, "MY_dir.p"]*2, res1[, "MY.p"])

  res2 <- modify_hitman_pvalues(tab=tt, overall.sign = 1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_gte(res2[, "MY_dir.p"], 0.5)

  tt2 <- tt
  tt2[, "EM.t"] <- -1*tt[, "EM.t"]
  res3 <- modify_hitman_pvalues(tab=tt2, overall.sign = -1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_gte(res3[, "MY_dir.p"], 0.5)

  res4 <- modify_hitman_pvalues(tab=tt2, overall.sign = 1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_equal(res4[, "MY_dir.p"]*2, res4[, "MY.p"])
})

test_that("another gene", {
  tt <- matrix(c(0.149, 0.00578, -0.215, 0.149, 7.09, 0.00578), nrow=1,
               dimnames=list("gene1", c('EM_dir.p', 'MY_dir.p', 'EM.t', 'EM.p', 'MY.t', 'MY.p')))
  res1 <- modify_hitman_pvalues(tab=tt, overall.sign = -1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_equal(nrow(res1), 1)
  expect_equal(res1[, "EM_dir.p"]*2, res1[, "EM.p"])

  res2 <- modify_hitman_pvalues(tab=tt, overall.sign = 1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_gte(res2[, "EM_dir.p"], 0.5)

  tt2 <- tt
  tt2[, "EM.t"] <- -1*tt[, "EM.t"]
  res3 <- modify_hitman_pvalues(tab=tt2, overall.sign = -1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_gte(res3[, "EM_dir.p"], 0.5)

  res4 <- modify_hitman_pvalues(tab=tt2, overall.sign = 1, stat.cols=c("EM.t", "MY.t"), p.cols=c("EM_dir.p", "MY_dir.p"))
  expect_equal(res4[, "EM_dir.p"]*2, res4[, "EM.p"])
})
