#' Computing t*
#'
#' Computes the t* U-statistic for input data pairs
#' (x_1,y_1), (x_2,y_2), ..., (x_n,y_n)
#' using the improved algorithm based on Heller and Heller (2016) <arXiv:1605.08732>
#' building off of the work of Weihs, Drton, and Leung (2015)
#' <DOI:10.1007/s00180-015-0639-x>.
#'
#' @export
#'
#' @param x A numeric vector of x values (length >= 4).
#' @param y A numeric vector of y values, should be of the same length as x.
#'
#' @return The numeric value of the t* statistic.
#'
#' @references
#' Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based
#' on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014),
#' no. 2, 1006--1028.
#' \cr\cr
#' Heller, Yair and Heller, Ruth. "Computing the Bergsma Dassios
#' sign-covariance." arXiv preprint arXiv:1605.08732 (2016).
#' \cr\cr
#' Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
#' Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
#'
#' @examples
#' \dontrun{
#' library(TauStar)
#'
#' # Compute t* for a concordant quadruple
#' tStar(c(1,2,3,4), c(1,2,3,4)) # == 2/3
#'
#' # Compute t* for a discordant quadruple
#' tStar(c(1,2,3,4), c(1,-1,1,-1)) # == -1/3
#'
#' # Compute t* on random normal iid normal data
#' set.seed(23421)
#' tStar(rnorm(4000), rnorm(4000)) # near 0
#'
#' # Compute t* as a v-statistic
#' set.seed(923)
#' tStar(rnorm(100), rnorm(100), vStatistic = TRUE)
#'
#' # Compute an approximation of tau* via resampling
#' set.seed(9492)
#' tStar(rnorm(10000), rnorm(10000), resample = TRUE, sampleSize = 30,
#'       numResamples = 5000)
#' }
tStar <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Input x and y to tStar must be numeric.")
  }
  if (length(x) != length(y) || length(x) < 4) {
    stop("Input x and y to tStar are of the wrong length, they must both have equal length < 4.")
  }
  return(TStarHellerAndHellerRCPP(x, y))
}
#' Check if Vector of Probabilities
#'
#' Checks if the input vector has entries that sum to 1 and are non-negative
#'
#' @param probs the probability vector to check
#'
#' @return TRUE if conditions are met, FALSE if otherwise
isProbVector <- function(probs) {
  probSum = sum(probs)
  return(abs(probSum - 1) < 10^-7 && all(probs >= 0))
}

#' Null asymptotic distribution of t* in the mixed case
#'
#' Density, distribution function, quantile function and random generation for
#' the asymptotic null distribution of t* in the mixed case. That is, in the
#' case that t* is generated a sample from an independent bivariate distribution
#' where one coordinate is marginally discrete and the other marginally
#' continuous.
#'
#' @export
#'
#' @param x the value (or vector of values) at which to evaluate the function.
#' @param error a tolerated error in the result. This should be considered as a
#'        guide rather than an exact upper bound to the amount of error.
#' @param lower.tail a logical value, if TRUE (default), probabilities are
#'        \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.
#' @param probs a vector of probabilities corresponding to the (ordered)
#'        support the marginally discrete random variable. That is, if the
#'        marginally discrete distribution has support \eqn{u_1,...,u_n}
#'        then the ith entry of probs should be the probability of seeing
#'        \eqn{u_i}.
#'
#' @name pMixHoeffInd
#' @rdname pMixHoeffInd
#'
#' @return pMixHoeffInd gives the distribution
#'         function.
pMixHoeffInd <- function(x, probs, lower.tail = T, error = 0.01) { # 10^-6
  if (!isProbVector(probs)) {
    stop("probs in pMixHoeffInd is not a probability vector.")
  }
  if (any(probs == 1)) {
    lowerTailProb = if (x >= 0) 1 else 0
  } else {
    eigenP = eigenForDiscreteProbs(probs)
    lowerTailProb = HoeffIndMixedCdfRCPP(x + 2 * sum(eigenP), eigenP, error)
  }
  if (lower.tail) {
    return(lowerTailProb)
  }
  return(1 - lowerTailProb)
}

#' Determine if input data is discrete
#'
#' Attempts to determine if the input data is from a discrete distribution. Will
#' return true if the data type is of type integer or there are non-unique
#' values.
#'
#' @param x a vector which should be determined if discrete or not.
#'
#' @return the best judgement of whether or not the data was discrete
isDiscrete = function(x) {
  if (is.integer(x) || (length(unique(x)) != length(x))) {
    return(T)
  }
  return(F)
}

#' Is Vector Valid Data?
#'
#' Determines if input vector is a valid vector of real valued observations
#'
#' @param x the vector to be tested
#'
#' @return TRUE or FALSE
isValidDataVector <- function(x) {
  return(is.numeric(x) || is.integer(x))
}

#' Test of Independence Using the Tau* Measure
#'
#' Performs a (consistent) test of independence between two input vectors using
#' the asymptotic (or permutation based) distribution of the test statistic t*.
#' The asymptotic results hold in the case that x is generated from either a
#' discrete or continous distribution and similarly for y (in particular it is
#' allowed for one to be continuous while the other is discrete). The asymptotic
#' distributions were computed in Nandy, Weihs, and Drton (2016)
#' <http://arxiv.org/abs/1602.04387>.
#'
#' @export
#'
#' @param x a vector of sampled values.
#' @param y a vector of sampled values corresponding to x, y must be the same
#'        length as x.
#'
#' @return a list with class "tstest" recording the outcome of the test.
#'
#' @references
#' Preetam Nandy, Luca Weihs, and Mathias Drton. Large-Sample Theory for the
#' Bergsma-Dassios Sign Covariance. arXiv preprint arXiv:1602.04387. 2016.
#'
#' @examples
#' set.seed(123)
#' x = rnorm(100)
#' y = rbinom(100, 1, 0.5)
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # big p-value
#'
#' y = y + x # make x and y correlated
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # small p-value
tauStarTest <- function(x, y, error = 0.01, mode = "auto") {
  if (!isValidDataVector(x) || !isValidDataVector(y) ||
      length(x) != length(y)) {
    stop(paste("vectors inputted to tauStarTest must be of type numeric or",
               "integer and must be the same length"))
  }
  if (mode != "auto" && mode != "permute" && mode != "asymptotic"){
    stop("Mode must be auto, permute or asymptotic")
  }  

  xIsDis = isDiscrete(x)
  yIsDis = isDiscrete(y)
  if (xIsDis && !yIsDis) {
    z = x
    x = y
    y = z
  }
  xIsDis = isDiscrete(x)
  yIsDis = isDiscrete(y)
  n = length(x)
  
  toReturn = list()
  class(toReturn) = "tstest"
  toReturn$x = x
  toReturn$y = y
  toReturn$tStar = tStar(x, y)
  
  if (yIsDis){
    if (length(unique(y)) == 1){
      toReturn$pVal = 1
      return(toReturn)
    }  
    if (length(unique(y)) <= 50) 
      p = as.numeric(table(y))/n
    else {
      lbr = quantile(y, 0:50/50)
      p = hist(y, breaks = unique(lbr), plot = F)$counts/n
    }

    if (mode == "permute" || (mode == "auto" && any(p > 0.99))){
      taustar_p = NULL
      for(i in c(1:100)){
        temp = tStar(sample(x), y)
        taustar_p = c(taustar_p, temp)
      }
      toReturn$pVal = 1 - sum(toReturn$tStar >= taustar_p) / 100
    }
    else{
      toReturn$pVal = 1 - pMixHoeffInd(n * toReturn$tStar, probs = p, error = error)
    }  
      
  } else {
    toReturn$pVal = NA
  }
 
  return(toReturn)
}

#' @title Perform Differential Expression Test by tauStar on One Gene
#'
#' @description
#' Test if one gene is differentially expressed along pseudotime.
#' @param gene A vector of pseudotime.
#' @param count.v Rhe expression data.
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{taustar}}{The estimated Bergsma-Dassios Correlation Coefficient}
#'   \item{\code{taustar.pv}}{The p-value}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom stats fitted
#'
#' @export pseudotimeDE
#'
#' @author Yuheng Lai
tauStarDE <- function(pseudotime,
                      count.v, mode = "auto") {
  res = tauStarTest(pseudotime, count.v, mode)
  tau = res$tStar
  tau_pv = res$pVal

  return(list(taustar = tau,
              taustar.pv = tau_pv
  ))
}

#' @title Perform Differential Expression Test by tauStar on a List of Genes
#'
#' @description
#' Test if each gene in a list of genes is differentially expressed along pseudotime.
#' A wrapper of \code{tauStarDE::tauStarDE}.
#'
#' @param gene.vec A vector of genes. It should be a subset of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param mat The input expression data. It can be:
#' (1) A SingleCellExperiment object which contain the expression data;
#' (2) An matrix;
#' (3) A Seurat object which contain the expression data.
#' Its row names should be genes and col names should be cells.
#' @param assay.use The \code{assay} used in SingleCellExperiment or \code{slot} used in Seurat. Default is \code{counts}.
#' @param seurat.assay The \code{assay} used in Seurat. Default is \code{'RNA'}.
#' @param mc.cores Number of cores for computing.
#' @param mc.preschedule See \code{mclapply}. Default is TRUE.
#' @param SIMPLIFY A logic variable whether to return a tibble (TRUE) or a list of lists (FALSE). Default is TRUE.
#' @return A tibble of summary results of genes
#'
#' @examples
#' data("LPS_sce")
#' data("LPS_ori_tbl")
#' data("LPS_sub_tbl")
#' res <- PseudotimeDE::runTauStarDE(gene.vec = c("CCL5", "CXCL10"),
#' ori.tbl = LPS_ori_tbl, sub.tbl = LPS_sub_tbl[1:10], mat = LPS_sce, model = "nb")
#'
#' @export runTauStarDE
#' @author Yuheng Lai

runTauStarDE <- function(gene.vec,
                         ori.tbl,
                         mat,
                         assay.use = "counts",
                         seurat.assay = 'RNA',
                         mode = 'auto',
                         mc.cores = 2,
                         mc.preschedule = TRUE,
                         SIMPLIFY = TRUE) {
  ## ori.tbl is a dataframe containing cell and pseudotime. "cell" must be the same order as in mat
  num_cell <- length(ori.tbl$cell)
  pseudotime <- ori.tbl$pseudotime
  
  num_total_cell <- dim(mat)[2]
  
  if(class(mat)[1] == "SingleCellExperiment") {
    expv <- SummarizedExperiment::assay(mat, assay.use)[gene.vec, ori.tbl$cell]
    count.v <- expv
  }
  else if(class(mat)[1] == "Seurat") {
    if(assay.use == 'logcounts'){
      assay_alter <- 'data'
    }else{
      assay_alter <- assay.use
    }
    expv <- Seurat::GetAssayData(object = mat, slot = assay_alter, assay = seurat.assay)[gene.vec, ori.tbl$cell]
    count.v <- expv
  }
  else {
    expv <- mat[gene.vec, ori.tbl$cell]
    count.v <- expv
  }
  
  
  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- mc.cores
  
  res <- BiocParallel::bplapply(gene.vec, function(x, ...) {
    cur_res <- tryCatch(expr = tauStarDE(pseudotime = pseudotime,
                                         count.v = count.v[x, ], mode = mode), #input only the target gene
                        error = function(e) {
                          list(taustar = NA,
                               taustar.pv = NA)
                        })
    cur_res
  },
  BPPARAM = BPPARAM)
  
  
  if(SIMPLIFY) {
    res <- simplify2array(res)
    res <- t(res)
    rownames(res) <- gene.vec
    res <- tibble::as_tibble(res, rownames = "gene")
    res <- tidyr::unnest(res, cols = c("taustar", "taustar.pv"))
  }
  
  res
}
