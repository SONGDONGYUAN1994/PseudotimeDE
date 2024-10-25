context("Run DE")
library(PseudotimeDE)

test_that("PseudotimeDE works", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")

  res1 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                    ori.tbl = LPS_ori_tbl,
                                    sub.tbl = LPS_sub_tbl[1:100],
                                    mat = LPS_sce,
                                    model = "nb")
  expect_equal(length(res1), 12)
  expect_false(is.na(res1$test.statistics))

  res2 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                    ori.tbl = LPS_ori_tbl,
                                    sub.tbl = NULL,
                                    mat = LPS_sce,
                                    model = "auto")
  expect_equal(length(res2), 12)
  expect_false(is.na(res2$test.statistic))

  res3 <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                     ori.tbl = LPS_ori_tbl,
                                     sub.tbl = LPS_sub_tbl[1:100],
                                     mat = LPS_sce,
                                     model = "auto",
                                     mc.cores = 1)
  
  expect_equal(dim(res3)[1], 3)
  expect_false( any(is.na(res3$test.statistics[1:2])) )
  expect_true( is.na(res3$test.statistics[[3]]) )
  expect_contains(class( res3$notes[[3]] ),
                  "error")
  

  res4 <- PseudotimeDE::plotCurve(gene.vec = c("CCL5", "CXCL10"),
                                        ori.tbl = LPS_ori_tbl,
                                        mat = LPS_sce,
                                        model.fit = res3$gam.fit[1:2])
  expect_equal(class(res4)[1], "gg")

  res5 <- PseudotimeDE::plotUncertainty(ori.tbl = LPS_ori_tbl,
                                        sub.tbl = LPS_sub_tbl[1:100])
  expect_equal(class(res5)[1], "gg")
})

