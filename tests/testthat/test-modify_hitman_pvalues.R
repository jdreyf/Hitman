context("modify_hitman_pvalues")

#gene at bottom of list
#my.p > em.p & wrong sign --> flip my.p sign
test_that("gene88", {
  v <- ret["gene88",] #a matrix
  prod.sign <- sign(v[,"EM.t"])*sign(v[,"MY.t"])
  expect_false(prod.sign == ey.sign)

  expect_gt(v[,"MY.t"], 0)
  expect_gt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"EM.p"], v[,"EM_dir.p"])
  alt.tmp <- ifelse (v[,"EM.t"] > 0, "greater", "less")
  expect_equal(as.numeric(v["MY_dir.p"]),
               as.numeric(ezlimma:::two2one_tailed(v, stat.col = "MY.t", p.col="MY.p", alternative = alt.tmp)))
})

#my.p < em.p, & wrong sign
test_that("gene72", {
  v <- ret["gene72",] #a matrix
  prod.sign <- sign(v[,"EM.t"])*sign(v[,"MY.t"])
  expect_false(prod.sign == ey.sign)

  expect_lt(v[,"MY.t"], 0)
  expect_lt(v[,"MY.p"], v[,"EM.p"])
  expect_equal(v[,"MY.p"], v[,"MY_dir.p"])

  alt.tmp <- ifelse (v[,"EM.t"] > 0, "less", "greater")
  expect_equal(as.numeric(v["EM_dir.p"]),
               as.numeric(ezlimma:::two2one_tailed(v, stat.col = "EM.t", p.col="EM.p", alternative = alt.tmp)))
})
