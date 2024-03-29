% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted.R
\name{weighted_quantile}
\alias{weighted_quantile}
\title{Naive weighted sample quantile estimator}
\usage{
weighted_quantile(x, weights, prob, type = 1, na_rm = FALSE)
}
\arguments{
\item{x}{a data vector}

\item{weights}{a weights vector}

\item{prob}{requested quantiles}

\item{type}{Hyndman-Fan quantile type}

\item{na_rm}{na_rm}
}
\description{
Weighted sample quantile estimators Types 1-3 according to the Hyndman-Fan
typology.

Sample quantiles may not always be the most efficient, and a Type 1 estimator
is not even unbiased, but given the large sample sizes of PISA and the BRR
methodology (in part designed to allow for naive quantile computations),
these are not relevant concerns and it is much faster than the Harell-Davis
estimator.

Type 2 and 3 implementations are approximate algorithms that do not
completely match the behavior of \code{stats::quantile}:
\itemize{
\item type 1 is the inverse cdf
\item type 2 takes the average even when there is no tie
\item type 3 takes the nearest quantile instead of the nearest even
}
}
