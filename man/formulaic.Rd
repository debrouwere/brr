% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functional.R
\name{formulaic}
\alias{formulaic}
\title{Wrap a statistical function in a formula interface}
\usage{
formulaic(fn)
}
\arguments{
\item{fn}{}
}
\description{
Wraps a function that has a \code{data} argument in first position with a
\code{wrapper(formula, data, ...)} interface. The formula argument will be used to
pull a column from the data, which will then be passed on to the wrapped
function.

\code{formulaic} is useful to pass a "plain" statistical function to \code{brr} that
does not have a formula interface, with much less overhead compared to
\code{weighted_aggregate}
}
