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


test_that("runPseudotimeDE works with Seurat object", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")
  
  requireNamespace("Seurat")
  
  suppressWarnings(
    LPS_seurat <- SeuratObject::as.Seurat(LPS_sce) |>
      SeuratObject::RenameAssays(assay.name = "originalexp",
                                 new.assay.name = "RNA")
  )
  
  
  res_sce <- runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                             ori.tbl = LPS_ori_tbl,
                             sub.tbl = LPS_sub_tbl[1:100],
                             mat = LPS_sce,
                             model = "nb",
                             mc.cores = 1)
  
  
  res_seurat <- runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                ori.tbl = LPS_ori_tbl,
                                sub.tbl = LPS_sub_tbl[1:100],
                                mat = LPS_seurat,
                                model = "nb",
                                mc.cores = 1)
  
  stable_colnames <- c("fix.pv", "emp.pv", "rank", "test.statistics", "aic", "expv.mean", "expv.zero")
  
  expect_equal(res_sce[1:2, stable_colnames],
               res_seurat[1:2, stable_colnames])
  expect_contains(class(res_seurat$notes[[3]]),
                  "error")
})


test_that("runPseudotimeDE works with matrix input", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")
  
  LPS_count_mat <- SingleCellExperiment::counts(LPS_sce)
  
  res_sce <- runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                             ori.tbl = LPS_ori_tbl,
                             sub.tbl = LPS_sub_tbl[1:100],
                             mat = LPS_sce,
                             model = "nb",
                             mc.cores = 1)
  
  
  res_count_mat <- runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                   ori.tbl = LPS_ori_tbl,
                                   sub.tbl = LPS_sub_tbl[1:100],
                                   mat = LPS_count_mat,
                                   model = "nb",
                                   mc.cores = 1)
  
  expect_equal(res_sce[1:2, ], res_count_mat[1:2, ])
  expect_contains(class(res_count_mat$notes[[3]]),
                  "error")
})


test_that("We can change parameters", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")
  
  res_k6 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb")
  
  expect_identical(as.character(formula(res_k6$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')
  
  res_k7 <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb",
                                       k = 7,
                                       knots = c(0:6/6))
  
  expect_identical(as.character(formula(res_k7$gam.fit))[[3]],
                   's(pseudotime, k = 7, bs = "cr")')
  
  res_run_k6 <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "auto",
                                              mc.cores = 1)
  
  expect_identical(as.character(formula(res_run_k6$gam.fit[[1]]))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')
  
  
  res_run_k7 <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "auto",
                                              k = 7,
                                              knots = c(0:6/6),
                                              mc.cores = 1)
  
  expect_identical(as.character(formula(res_run_k7$gam.fit[[1]]))[[3]],
                   's(pseudotime, k = 7, bs = "cr")')
  
  
})


test_that("qgam works",{
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")
  
  res_nb <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb")
  expect_equal(length(res_nb), 12)
  expect_false(is.na(res_nb$test.statistics))
  
  capture_output(
    res_qgam <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                           ori.tbl = LPS_ori_tbl,
                                           sub.tbl = LPS_sub_tbl[1:100],
                                           mat = LPS_sce,
                                           model = "qgam")
  )
  
  expect_equal(length(res_qgam), 12)
  expect_false(is.na(res_qgam$test.statistic))
  
  expect_equal(res_nb$expv.quantile, res_qgam$expv.quantile)
  
  
  
  
  res_run_nb <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "nb",
                                              mc.cores = 1)
  
  expect_equal(dim(res_run_nb)[1], 3)
  expect_false( any(is.na(res_run_nb$test.statistics[1:2])) )
  expect_true( is.na(res_run_nb$test.statistics[[3]]) )
  expect_contains(class( res_run_nb$notes[[3]] ),
                  "error")
  
  capture_output(
    res_run_qgam <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                                  ori.tbl = LPS_ori_tbl,
                                                  sub.tbl = LPS_sub_tbl[1:100],
                                                  mat = LPS_sce,
                                                  model = "qgam",
                                                  mc.cores = 1)
  )
  
  expect_equal(dim(res_run_qgam)[1], 3)
  expect_false( any(is.na(res_run_qgam$test.statistic[1:2])) )
  expect_true( is.na(res_run_qgam$test.statistic[[3]]) )
  expect_contains(class( res_run_qgam$notes[[3]] ),
                  "error")
  
  expect_equal(res_run_nb$expv.quantile,
               res_run_qgam$expv.quantile)
  
})



test_that("We can change spline basis", {
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")


  # just pseudotimeDE

  res_cr <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb")

  expect_identical(as.character(formula(res_cr$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')



  res_cr_formula <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb",
                                       formula = expv ~ s(pseudotime, k = 6,bs = 'cr') )

  expect_identical(as.character(formula(res_cr_formula$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')
  
  expect_equal(res_cr, res_cr_formula)
  
  
  res_cs <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                               ori.tbl = LPS_ori_tbl,
                                               sub.tbl = LPS_sub_tbl[1:100],
                                               mat = LPS_sce,
                                               model = "nb",
                                               formula = expv ~ s(pseudotime, k = 6, bs = 'cs') )
  
  expect_identical(as.character(formula(res_cs$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cs")')

  expect_equal(res_cr$expv.quantile,
               res_cs$expv.quantile)

  
  res_cc <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       formula = expv ~ s(pseudotime, k = 6, bs = 'cc') ,
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb")

  expect_identical(as.character(formula(res_cc$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cc")')

  expect_equal(res_cr$expv.quantile,
               res_cc$expv.quantile)

  
  
  
  # test in runPseudotimeDE
  
  res_run_nb <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "nb",
                                              mc.cores = 1)
  
  expect_equal(dim(res_run_nb)[1], 3)
  expect_false( any(is.na(res_run_nb$test.statistics[1:2])) )
  expect_true( is.na(res_run_nb$test.statistics[[3]]) )
  expect_contains(class( res_run_nb$notes[[3]] ),
                  "error")
  
  
  res_run_cs <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              formula = expv ~ s(pseudotime, k = 6, bs = 'cs'),
                                              mat = LPS_sce,
                                              model = "nb",
                                              mc.cores = 1)
  
  
  expect_equal(dim(res_run_cs)[1], 3)
  expect_false( any(is.na(res_run_cs$test.statistics[1:2])) )
  expect_true( is.na(res_run_cs$test.statistics[[3]]) )
  expect_contains(class( res_run_cs$notes[[3]] ),
                  "error")
  
  expect_equal(res_run_nb$expv.quantile,
               res_run_cs$expv.quantile)

})



test_that("We can use adaptive smoothing in PseudotimeDE", {
  
  skip(message = "While possible, too slow/memory-hungry")
  
  data("LPS_sce")
  data("LPS_ori_tbl")
  data("LPS_sub_tbl")
  
  
  # just pseudotimeDE
  
  res_cr <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb")
  
  expect_identical(as.character(formula(res_cr$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')
  
  
  
  
  res_ad_cr <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                          ori.tbl = LPS_ori_tbl,
                                          sub.tbl = LPS_sub_tbl[1:100],
                                          mat = LPS_sce,
                                          model = "nb",
                                          formula = expv ~ s(pseudotime, k = 12, bs = 'ad', xt = list(bs = 'cr')),
                                          k = 12,
                                          knots = c(0:11/11))
  
  expect_identical(as.character(formula(res_ad_cr$gam.fit))[[3]],
                   's(pseudotime, k = 12, bs = "ad", xt = list(bs = "cr"))')
  
  expect_equal(res_cr$expv.quantile,
               res_ad_cr$expv.quantile)
  
  
  
})


