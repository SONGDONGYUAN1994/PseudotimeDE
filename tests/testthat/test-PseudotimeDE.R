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
  
  res_cs <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb",
                                       bs = "cs")
  
  expect_identical(as.character(formula(res_cs$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cs")')
  
  
  
  res_cc <- PseudotimeDE::pseudotimeDE(gene = "CCL5",
                                       ori.tbl = LPS_ori_tbl,
                                       sub.tbl = LPS_sub_tbl[1:100],
                                       mat = LPS_sce,
                                       model = "nb",
                                       bs = "cc")
  
  expect_identical(as.character(formula(res_cc$gam.fit))[[3]],
                   's(pseudotime, k = 6, bs = "cc")')
  
  
  
  # test runPseudotimeDE
  
  res_run_cr <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "auto",
                                              mc.cores = 1)
  
  expect_identical(as.character(formula(res_run_cr$gam.fit[[1]]))[[3]],
                   's(pseudotime, k = 6, bs = "cr")')
  
  
  res_run_cs <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "auto",
                                              bs = "cs",
                                              mc.cores = 1)
  
  expect_identical(as.character(formula(res_run_cs$gam.fit[[1]]))[[3]],
                   's(pseudotime, k = 6, bs = "cs")')
  
  expect_true(is.na( res_run_cs$test.statistics[[3]] ))
  expect_identical(res_run_cr$expv.mean, res_run_cs$expv.mean)
  
  relative_diff <- (res_run_cr$test.statistics - res_run_cs$test.statistics)/res_run_cr$test.statistics
  expect_true(all( relative_diff[1:2] < .1 ))
  
  
  res_run_cc <- PseudotimeDE::runPseudotimeDE(gene.vec = c("CCL5", "CXCL10", "JustAJoke"),
                                              ori.tbl = LPS_ori_tbl,
                                              sub.tbl = LPS_sub_tbl[1:100],
                                              mat = LPS_sce,
                                              model = "auto",
                                              bs = "cc",
                                              mc.cores = 1)
  
  expect_identical(as.character(formula(res_run_cc$gam.fit[[1]]))[[3]],
                   's(pseudotime, k = 6, bs = "cc")')
  
  expect_true(is.na( res_run_cc$test.statistics[[3]] ))
  
  expect_identical(res_run_cr$expv.mean, res_run_cc$expv.mean)
  
})


