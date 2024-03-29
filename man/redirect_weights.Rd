% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional.R
\name{redirect_weights}
\alias{redirect_weights}
\title{Wrap a function that requires the name of the weights column so that instead
it accepts a vector of weights, and/or rename the weights argument.}
\usage{
redirect_weights(f, by_reference = FALSE, name = "weights")
}
\arguments{
\item{f}{f}

\item{by_reference}{by_reference}

\item{name}{name}
}
\value{
a wrapped function
}
\description{
This function makes it possible to modify fitters with unusual arguments
so they match the signature expected by \code{brr:brr}.
}
\examples{
df <- data.frame(...)
weights <- label_vector(df$weights, "weights")

lm_by_ref <- function(formula, data, weights) {
  weights <- data[[weights]]
  lm(formula, data, weights = weights)
}

redirect_weights(lm_by_ref, by_reference = TRUE)(y ~ x, data = df, weights = weights)

lm_wt <- function(formula, data, wt) {
  lm(formula, data, weights = wt)
}

redirect_weights(lm_wt, name = "wt")(y ~ x, data = df, weights = weights)
}
