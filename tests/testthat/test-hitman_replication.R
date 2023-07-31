context("hitman_replication")

test_that("correct significances", {
  hmr <- hitman_replication(tab=tab.tmp)
  expect_equal(hmr["r10", grep("p$", colnames(hmr))], 0.5*max(tab.tmp["r10", c(2, 4)]))
  # r3 is inconsistent
  expect_equal(hmr["r3", grep("p$", colnames(hmr))], 1)

  hmby <- hitman_replication(tab=tab.tmp, fdr.method = "BY")
  expect_true(all(hmby$FDR >= hmr$FDR))
})


test_that("invalid input", {
  expect_error(tab.tmp |> dplyr::mutate(stat1 = as.character(stat1)) |>
    hitman_replication())
  expect_error(tab.tmp |> hitman_replication(cols=2:5))
  expect_error(tab.tmp |> hitman_replication(cols=colnames(tab.tmp)[-1]))

  tab.na <- tab.tmp
  tab.na[1,1] <- NA
  expect_error(tab.na |> hitman_replication())

  # only matrices can have NULL row names
  tab2 <- tab.tmp |> as.matrix()
  rownames(tab2) <- NULL
  expect_error(tab2 |> hitman_replication(reorder.rows = TRUE))
})

test_that("matches DMT", {
  dmtr <- DirectionalMaxPTest::dmt(tab=tab.tmp)
  hmr <- hitman_replication(tab=tab.tmp)
  expect_true(all(dmtr == hmr))
})
