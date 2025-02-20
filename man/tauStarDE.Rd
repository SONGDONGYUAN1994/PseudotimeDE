% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tauStarDE.R
\name{tauStarDE}
\alias{tauStarDE}
\title{Perform Differential Expression Test by tauStar on One Gene}
\usage{
tauStarDE(pseudotime, count.v)
}
\arguments{
\item{count.v}{Rhe expression data.}

\item{gene}{A vector of pseudotime.}
}
\value{
A list with the components:
\describe{
  \item{\code{taustar}}{The estimated Bergsma-Dassios Correlation Coefficient}
  \item{\code{taustar.pv}}{The p-value}
}
}
\description{
Test if one gene is differentially expressed along pseudotime.
}
\author{
Yuheng Lai
}
