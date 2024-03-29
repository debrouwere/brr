% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brr.R
\name{confint.brr}
\alias{confint.brr}
\title{Confidence intervals for model parameters of the balanced repeated model fit replications}
\usage{
\method{confint}{brr}(
  replications,
  level = 0.95,
  perturbation = 0.5,
  imputation = TRUE,
  extra = FALSE,
  na_rm = FALSE
)
}
\arguments{
\item{replications}{BRR replications}

\item{level}{the confidence level required}

\item{perturbation}{perturbation}

\item{imputation}{Incorporate imputation variance due to plausible values. TRUE by default but can be disabled when there's only a single outcome.}

\item{extra}{include standard error and variance components in the output}
}
\description{
Confidence intervals for model parameters of the balanced repeated model fit replications
}
