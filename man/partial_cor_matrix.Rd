% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted.R
\name{partial_cor_matrix}
\alias{partial_cor_matrix}
\title{Partial correlation matrix}
\usage{
partial_cor_matrix(data, center = TRUE, method = "unbiased", na_rm = FALSE)
}
\arguments{
\item{data}{matrix}

\item{center}{center}

\item{method}{method}

\item{na_rm}{na_rm, equivalent to \code{cor(use = "complete.obs")}, as pairwise complete etc. does not make sense for partial correlations}
}
\description{
Partial correlation matrix
}
