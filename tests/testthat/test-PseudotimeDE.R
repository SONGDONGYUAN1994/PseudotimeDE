context("Run DE")
library(PseudotimeDE)

test_that("PseudotimeDE works", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")

  res1 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                    ori.tbl = LPS_ori_tbl,
                                    sub.tbl = LPS_sub_tbl[1:100],
                                    sce = LPS_sce,
                                    model = "nb")
  expect_equal(length(res1), 10)

  res2 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                    ori.tbl = LPS_ori_tbl,
                                    sub.tbl = NULL,
                                    sce = LPS_sce,
                                    model = "auto")
  expect_equal(length(res2), 10)

  res3 <- PseudotimeDE::runPseudotimeDE(gene = c("CCL5", "CXCL10"),
                                     ori.tbl = LPS_ori_tbl,
                                     sub.tbl = NULL,
                                     sce = LPS_sce,
                                     model = "auto")
  expect_equal(length(res3), 2)


})
