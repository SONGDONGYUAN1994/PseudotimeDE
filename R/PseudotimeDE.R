#' @title Perform Differential Expression Test on One Gene
#'
#' @description
#' Test if one gene is differentially expressed along pseudotime.
#' @param gene A string of gene name. It should be one of the row names in sce.
#' @param ori.tbl A tibble or dataframe which contains the original cells and pseudotime as two columns.
#' @param sub.tbl A list of tibbles or dataframes where each is the fit of a subsample. Each element is the same format as ori.tbl.
#' @param sce A SingleCellExperment object which contain the count data. Its row names should be genes and col names should be cells.
#' @param model A string of the model name. One of \code{nb}, \code{zinb} and \code{auto}.
#' @param k A integer of the basis dimension. Default is 6. The reults are usually robust to different k; we recommend to use k from 5 to 10.
#' @param knots A numeric vector of the location of knots.
#' @param fix.weight A logic variable indicating if the ZINB-GAM will use the zero weights from the original model.
#' @param aicdiff A numeric variable of the threshold of modle selection. Only works when \code{model = `auto`}.
#' @param seed A numeric variable of the random seed. It mainly affects the fitting of null distribution.
#'
#'
#' @return A list with the components:
#' \describe{
#'   \item{\code{fix.pv}}{The p-value assuming the pseudotime is fixed}
#'   \item{\code{emp.pv}}{The permutation p-value}
#'   \item{\code{para.pv}}{The permutation p-value by fitting a parametric null distribution}
#'   \item{\code{rank}}{The estimated effect degree of freedom of the original model}
#'   \item{\code{gam.fit}}{The fitted gam model on original data}
#'   \item{\code{zinf}}{Whether the model is zero inflated}
#'   \item{\code{aic}}{The AIC of the orignial model}
#'   \item{\code{expv.quantile}}{Quantiles of the log counts plus 1}
#'   \item{\code{expv.mean}}{Mean of the log counts plus 1}
#'   \item{\code{expv.zero}}{Zero proportions of counts}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assays colData
#' @importFrom methods is
#' @importFrom stats AIC binomial dgamma dnbinom fitted logLik model.frame model.matrix model.response pchisq pf pgamma predict qchisq quantile rbinom rgamma rmultinom update
#' @importFrom tibble as_tibble
#' @importFrom mgcv gam nb
#' @importFrom fitdistrplus fitdist
#' @importFrom mixtools gammamixEM
#'
#' @export pseudotimeDE
#'
#' @author Dongyuan Song


pseudotimeDE <- function(gene,
                         ori.tbl,
                         sub.tbl,
                         sce,
                         model = c("nb", "zinb", "auto"),
                         k = 6,
                         knots = c(0:5/5),
                         fix.weight = TRUE,
                         aicdiff = 10,
                         seed = 123) {

  ## ori.tbl is a dataframe containing cell and pseudotime. "cell" must be the same order as in sce
  stopifnot(is(sce, "SingleCellExperiment"))
  ## Set seed
  set.seed(seed)

  fit.formula <- stats::as.formula(paste0("expv ~ s(pseudotime, k = ", k, ",bs = 'cr')"))

  num_cell <- length(ori.tbl$cell)
  pseudotime <- ori.tbl$pseudotime

  expv <- assays(sce)$counts[gene, ori.tbl$cell]
  count.v <- expv

  dat <- cbind(pseudotime, expv) %>% as_tibble()

  expv.quantile <- quantile(log(count.v + 1))
  expv.mean <- mean(log(count.v + 1))
  expv.zero <- sum(count.v == 0)/length(count.v)

  fit.nb <- fit_gam(dat, zinf = FALSE, use_weights = FALSE, k = k, knots = knots)
  fit.zinb <- zinbgam(fit.formula, ~ logmu,
                      data = dat, knots = knots, k = k)
  aic.zinb <- fit.zinb$aic
  aic.nb <- AIC(fit.nb)

  ## zinb model only uses the mu part.
  fit.zinb <- fit.zinb$fit.mu

  if(model == "nb") {
    fit <- fit.nb
    zinf <- FALSE
    aic <- aic.nb
  }
  else if (model == "zinb"){
    fit <- fit.zinb
    zinf <- TRUE
    aic <- aic.zinb
  }
  else if (model == "auto") {
    if(aic.nb - aic.zinb > aicdiff) {
      fit <- fit.zinb
      zinf <- TRUE
      aic <- aic.zinb
    }
    else {
      fit <- fit.nb
      zinf <- FALSE
      aic <- aic.nb
    }
  }
  else {stop("Specified 'method=' parameter is invalid. Must be one of 'nb', 'zinb', 'auto'.")}



  ## ori p-value
  edf1 <- sum(fit$edf1)
  s.pv <- summary(fit)$s.pv

  ## Simon 2012 test statistic for sample. Copy from mgcv, but drop the subsampling part.
  res.df <- -1
  p <- fit$coefficients
  V <- fit$Vp
  rank <- sum(fit$edf1) - 1
  Xp <- model.matrix(fit)

  Tr <- testStat(p = p[2:k], X = Xp[, 2:k,drop=FALSE], V = V[2:k, 2:k,drop=FALSE], rank = rank, res.df = -1, type = 0)
  Tr <- Tr$stat

  if(is.null(sub.tbl)) {
    return(list(fix.pv = s.pv,
                emp.pv = NA,
                para.pv = NA,
                rank = rank,
                gam.fit = fit,
                zinf = zinf,
                aic = aic,
                expv.quantile = expv.quantile,
                expv.mean = expv.mean,
                expv.zero = expv.zero
    ))
  }

  ##  Subsample Tr
  n.boot <- length(sub.tbl)

  ## Redefine count.v since the original fit not necessarily includes all cell
  expv <- assays(sce)$counts[gene, ]
  count.v <- expv

  ## Use weights

  if(fix.weight && zinf) {
    ## Use the weights estimated from ori pseudotime

    cell_weights <- rep(1, num_cell)
    cell_weights <- fit.zinb$prior.weights
    names(cell_weights) <- colnames(sce)
  }
  else {
    ## All weights are 1
    cell_weights <- rep(1, num_cell)
    names(cell_weights) <- colnames(sce)}

  ## Start permutation
  boot_tbl <- tibble::tibble(id = seq_len(length(sub.tbl)), time.res = sub.tbl)

  boot_models <- boot_tbl %>%
    dplyr::mutate(splits = purrr::map(time.res, function(x){
      x <- cbind(expv = count.v[x$cell], pseudotime = base::sample(x$pseudotime), cellWeights = cell_weights[x$cell]) %>% as.data.frame(); x
    })) %>%
    dplyr::mutate(model = lapply(X = splits, FUN = fit_gam, zinf = zinf, use_weights = fix.weight, k = k, knots = knots),
                  stat = sapply(model, function(fit) {
                    if(is.logical(fit)) {Tr <- NA}
                    else {
                      edf1 <- sum(fit$edf1)
                      res.df <- -1
                      p <- fit$coefficients
                      V <- fit$Vp
                      rank <- sum(fit$edf1) - 1
                      Xp <- model.matrix(fit)

                      Tr <- testStat(p = p[2:k], X = Xp[, 2:k,drop=FALSE], V = V[2:k, 2:k,drop=FALSE], rank = rank, res.df = -1, type = 0)
                      Tr <- as.numeric(Tr$stat)}

                    return(Tr)
                  }))

  ## empirical permutation p-value
  s.pv.p <- (sum(Tr <= boot_models$stat)+1)/(n.boot+1)

  ## parametric permutation p-value
  smooth_pv <- suppressMessages(cal_pvalue(y = Tr, x = boot_models$stat, plot.fit = FALSE))

  return(list(fix.pv = s.pv,
              emp.pv = s.pv.p,
              para.pv = smooth_pv,
              rank = rank,
              gam.fit = fit,
              zinf = zinf,
              aic = aic,
              expv.quantile = expv.quantile,
              expv.mean = expv.mean,
              expv.zero = expv.zero
  ))
}


## Define the fit_gam function

fit_gam <- function(dat, nthreads = 1, zinf, use_weights, knots, k) {
  if(use_weights) {cellWeights <- dat$cellWeights}
  else cellWeights <- rep(1, dim(dat)[1])

  fit.formula <- stats::as.formula(paste0("expv ~ s(pseudotime, k = ", k, ",bs = 'cr')"))

  fit.gam <- tryCatch(expr = {
    if(!zinf) {
      fit <- gam(fit.formula, family = nb(link = "log"),
                 data = dat,
                 knots = list(pseudotime = knots), control = list(nthreads = nthreads))
    }
    else {
      if(!use_weights) {
        fit <- zinbgam(fit.formula, ~ logmu,
                       data = dat, knots = knots, k = k)
        fit <- fit$fit.mu
      }
      else {
        fit <- gam(fit.formula, family = nb(link = "log"),
                   data = dat,
                   knots = list(pseudotime = knots), control = list(nthreads = nthreads), weights = cellWeights)
      }
    }
    return(fit)
  },
  error = function(e){return(FALSE)})
  fit.gam
}




## Zero Inflated NB GAM. Refer to https://github.com/AustralianAntarcticDataCentre/zigam.
## Log ZINB density
dzinb.log <- function(x,mu,pi,shape) {
  logp <- log(1 - pi)+dnbinom(x,size=shape,mu=mu,log=T)
  logp[x==0] <- log(exp(logp[x==0])+(pi[x==0]))
  logp
}

## Main function

zinbgam <- function(mu.formula,
                    pi.formula,
                    data,
                    mu=NULL,
                    pi=NULL,
                    min.em=5,
                    max.em=50,
                    tol=1.0e-4,
                    k,
                    knots) {
  ## As matrix
  #data <- data.frame(data)
  ## Extract the response y
  mf <- model.frame(update(mu.formula,.~1),data=data)
  y <- model.response(mf)
  N <- length(y)

  ## Response for pi component is the weights
  pi.formula <- update(pi.formula, z ~ .)

  ## Get inital NB fit
  fit.gam <- gam(mu.formula,
                 data = data, family = nb(link = "log"), knots = list(pseudotime = knots))
  ## Set initial pi, mu
  if(is.null(mu)) mu <- fitted(fit.gam)
  if(is.null(pi)) pi <- mean(y>0)

  logL <- double(max.em)
  theta <- fit.gam$family$getTheta(TRUE)
  for(k in 1:max.em) {
    ## Evaluate weights for current iteration
    z <- ifelse(y==0,pi/(pi+ (1-pi)*dnbinom(0,size=theta,mu=mu)),0)
    ## Update the data (with mu and z)
    data$z <- z
    data$mu <- mu
    data$logmu <- log(mu)
    ## Update models for current iteration
    fit.pi <- suppressWarnings(gam(pi.formula,family=binomial(),data=data))
    fit.mu <- suppressWarnings(gam(mu.formula,weights=1-z,family=nb(link = log), data=data, knots = list(pseudotime = knots)))
    pi <- predict(fit.pi,type="response")
    mu <- predict(fit.mu,type="response")
    theta <- fit.mu$family$getTheta(TRUE)
    ## Evaluate likelihood
    logL[k] <- sum(dzinb.log(y,mu,pi,theta))
    if(k >= 2) {
      if(abs(logL[k]-logL[k-1]) < tol) {
        logL <- logL[1:k]
        break
      }
    }
  }

  ## Calculate degrees of freedom and aic
  df <- attr(logLik(fit.pi),"df")+attr(logLik(fit.mu),"df")
  aic <- 2*(df-logL[length(logL)])
  ## Return results
  fit <- list(fit.mu=fit.mu,
              fit.pi=fit.pi,
              mu=mu,
              pi=pi,
              z=z,
              aic=aic,
              logL=logL,
              theta=theta)
  class(fit) <- "zinbgam"
  fit
}




###### mgcv function #######
#### These functions copied from Simon Woods's mgcv package.
##  testStat
testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
  ## Implements Wood (2013) Biometrika 100(1), 221-228
  ## The type argument specifies the type of truncation to use.
  ## on entry `rank' should be an edf estimate
  ## 0. Default using the fractionally truncated pinv.
  ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
  ## res.df is residual dof used to estimate scale. <=0 implies
  ## fixed scale.

  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ## remove possible ambiguity from statistic...
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")

  k <- max(0,floor(rank))
  nu <- abs(rank - k)     ## fractional part of supplied edf
  if (type==1) { ## round up is more than .05 above lower
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  }

  if (nu>0) k1 <- k+1 else k1 <- k

  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}

  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]

  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
    if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
    b12 <- .5*nu*(1-nu)
    if (b12<0) b12 <- 0
    b12 <- sqrt(b12)
    B <- matrix(c(1,b12,b12,nu),2,2)
    ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
    B <- ev%*%B%*%ev
    eb <- eigen(B,symmetric=TRUE)
    rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
    vec1 <- vec
    vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
    vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
      t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  ## there is an ambiguity in the choise of test statistic, leading to slight
  ## differences in the p-value computation depending on which of 2 alternatives
  ## is arbitrarily selected. Following allows both to be computed and p-values
  ## averaged (can't average test stat as dist then unknown)
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2)
  d1 <- t(vec1)%*%(R%*%p)
  d1 <- sum(d1^2)
  ##d <- d1 ## uncomment to avoid averaging

  rank1 <- rank ## rank for lower tail pval computation below

  ## note that for <1 edf then d is not weighted by EDF, and instead is
  ## simply refered to a chi-squared 1

  if (nu>0) { ## mixture of chi^2 ref dist
    if (k1==1) rank1 <- val <- 1 else {
      val <- rep(1,k1) ##ed$val[1:k1]
      rp <- nu+1
      val[k] <- (rp + sqrt(rp*(2-rp)))/2
      val[k1] <- (rp - val[k])
    }

    if (res.df <= 0) pval <- (liu2(d,val) + liu2(d1,val))/2 else ##  pval <- davies(d,val)$Qq else
      pval <- (simf(d,val,res.df) + simf(d1,val,res.df))/2
  } else { pval <- 2 }
  ## integer case still needs computing, also liu/pearson approx only good in
  ## upper tail. In lower tail, 2 moment approximation is better (Can check this
  ## by simply plotting the whole interesting range as a contour plot!)
  if (pval > .5) {
    if (res.df <= 0) pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 else
      pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
  }
  list(stat=d,pval=min(1,pval),rank=rank)

} ## end of testStat

## Other functions within testStat()
liu2 <- function(x, lambda, h = rep(1,length(lambda)),lower.tail=FALSE) {
  # Evaluate Pr[sum_i \lambda_i \chi^2_h_i < x] approximately.
  # Code adapted from CompQuadForm package of Pierre Lafaye de Micheaux
  # and directly from....
  # H. Liu, Y. Tang, H.H. Zhang, A new chi-square approximation to the
  # distribution of non-negative definite quadratic forms in non-central
  # normal variables, Computational Statistics and Data Analysis, Volume 53,
  # (2009), 853-856. Actually, this is just Pearson (1959) given that
  # the chi^2 variables are central.
  # Note that this can be rubbish in lower tail (e.g. lambda=c(1.2,.3), x = .15)

  #  if (FALSE) { ## use Davies exact method in place of Liu et al/ Pearson approx.
  #    require(CompQuadForm)
  #    r <- x
  #    for (i in 1:length(x)) r[i] <- davies(x[i],lambda,h)$Qq
  #    return(pmin(r,1))
  #  }

  if (length(h) != length(lambda)) stop("lambda and h should have the same length!")

  lh <- lambda*h
  muQ <- sum(lh)

  lh <- lh*lambda
  c2 <- sum(lh)

  lh <- lh*lambda
  c3 <- sum(lh)

  s1 <- c3/c2^1.5
  s2 <- sum(lh*lambda)/c2^2

  sigQ <- sqrt(2*c2)

  t <- (x-muQ)/sigQ

  if (s1^2>s2) {
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  } else {
    a <- 1/s1
    delta <- 0
    l <- c2^3/c3^2
  }

  muX <- l+delta
  sigX <- sqrt(2)*a

  return(pchisq(t*sigX+muX,df=l,ncp=delta,lower.tail=lower.tail))

}

simf <- function(x,a,df,nq=50) {
  ## suppose T = sum(a_i \chi^2_1)/(chi^2_df/df). We need
  ## Pr[T>x] = Pr(sum(a_i \chi^2_1) > x *chi^2_df/df). Quadrature
  ## used here. So, e.g.
  ## 1-pf(4/3,3,40);simf(4,rep(1,3),40);1-pchisq(4,3)
  p <- (1:nq-.5)/nq
  q <- qchisq(p,df)
  x <- x*q/df
  pr <- sum(liu2(x,a)) ## Pearson/Liu approx to chi^2 mixture
  pr/nq
}


## Gamma Mixture Fit of Permutation p-value
# Gamma / Two Component Gamma Mixture Permutation p-value (parametric p-value)

### Mix gamma pdf
dgammamix <- function(x, lambda, gamma.pars, k = dim(gamma.pars)[2]) {
  if(k == 2) d <- lambda[1]*dgamma(x, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1]) + lambda[2]*dgamma(x, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2])
  else  if(k == 3) d <- lambda[1]*dgamma(x, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1]) +
      lambda[2]*dgamma(x, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2]) + lambda[3]*dgamma(x, shape = gamma.pars[1, 3], scale = gamma.pars[2, 3])
  else d <- NA
  d
}

### Mix gamma cdf
pgammamix <- function(x, lambda, gamma.pars, k = dim(gamma.pars)[2], lower.tail = FALSE) {
  if(k == 2) p <- lambda[1]*pgamma(x, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1], lower.tail = lower.tail) + lambda[2]*pgamma(x, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2], lower.tail = lower.tail)
  else if (k == 3) p <- lambda[1]*pgamma(x, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1], lower.tail = lower.tail) +
      lambda[2]*pgamma(x, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2], lower.tail = lower.tail) + lambda[3]*pgamma(x, shape = gamma.pars[1, 3], scale = gamma.pars[2, 3], lower.tail = lower.tail)
  else p <- NA
  p
}

### Mix gamma rv
rgammamix <- function(n = 1, lambda, gamma.pars, k = dim(gamma.pars)[2]) {
  if(k == 2) {
    z <- rbinom(n, size = 1, lambda[1])
    r <- z*rgamma(n, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1]) + (1 - z)*rgamma(n, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2])
  }
  else if (k == 3) {
    z <- rmultinom(n, size = 1, prob = c(lambda[1], lambda[2], lambda[3]))
    y1 <- rgamma(n, shape = gamma.pars[1, 1], scale = gamma.pars[2, 1])
    y2 <- rgamma(n, shape = gamma.pars[1, 2], scale = gamma.pars[2, 2])
    y3 <- rgamma(n, shape = gamma.pars[1, 3], scale = gamma.pars[2, 3])

    y <- rbind(y1, y2, y3)

    r <- colSums(z*y)
  }
  else r <- NA

  r
}

### Suppress print in mixtools
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

### Calculate p-value
cal_pvalue <- function(x, y, epsilon = 1e-4, p.thresh = 0.01, plot.fit = FALSE) {
  x <- as.vector(x)

  ## gamma fit
  fit1 <- suppressWarnings(fitdist(x, "gamma"))
  test.res <- goftest::ad.test(x, "pgamma", shape = fit1$estimate[1], rate = fit1$estimate[2])

  p <- pgamma(y, shape = fit1$estimate[1], rate = fit1$estimate[2], lower.tail = FALSE)

  ## gamma mix fit
  if(test.res$p.value < 1) {
    fit2 <- quiet(gammamixEM(x, maxit = 1000, k = 2, epsilon = epsilon, mom.start = TRUE, maxrestarts = 20))

    ## three random start
    for (i in 1:3) {
      set.seed(i)
      fit.temp <- quiet(gammamixEM(x, maxit = 1000, k = 2, epsilon = epsilon, maxrestarts = 20, lambda = c(i/10, 1-i/10)))
      if(fit2$loglik < fit.temp$loglik) fit2 <- fit.temp
    }

    test.res2 <- goftest::ad.test(x, "pgammamix", lambda = fit2$lambda, gamma.pars = fit2$gamma.pars, lower.tail = TRUE)

    ## LRT gamma vs gamma mix
    LRT.p <- pchisq(2*(fit2$loglik - fit1$loglik), df = 3, lower.tail = FALSE)

    if(LRT.p < p.thresh) {
      p <- pgammamix(y, lambda = fit2$lambda, gamma.pars = fit2$gamma.pars)
      if(test.res2$p.value < p.thresh) warning(paste0("Final fit does not pass AD test, with p-value ", test.res2$p.value))
    }
  }
  p
}
