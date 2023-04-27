context("hitman2_replication")

test_that("correct significances", {
  hmr <- hitman_replication(tab=tab.tmp, reorder.rows = FALSE, fdr.method = "BH")
  hm2r <- hitman2_replication(tab=tab.tmp, reorder.rows = FALSE, p.adj.rate = "FDR")
  hm2r.fwer <- hitman2_replication(tab=tab.tmp, reorder.rows = FALSE, p.adj.rate = "FWER")
  expect_true(all(rownames(hm2r) == rownames(tab.tmp)))
  expect_true(all.equal(hm2r[, 1:2], hmr[, 1:2]))
  expect_true(all.equal(hm2r.fwer[, 1:2], hmr[, 1:2]))
  # r10 has lowest p & FDR
  expect_lt(hm2r[10, "FDR"], hmr[10, "FDR"])
  expect_true(all(hm2r$FDR <= hm2r.fwer$FWER))
})

