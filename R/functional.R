library("tidyverse")
library("rlang")

#' Wrap a statistical function in a formula interface
#'
#' @description
#' Wraps a function that has a `data` argument in first position with a
#' `wrapper(formula, data, ...)` interface. The formula argument will be used to
#' pull a column from the data, which will then be passed on to the wrapped
#' function.
#'
#' `formulaic` is useful to pass a "plain" statistical function to `brr` that
#' does not have a formula interface, with much less overhead compared to
#' `weighted_aggregate`
#'
#' @param fn
#'
#' @export
formulaic <- function(fn) {
  function(formula, data, ...) {
    outcome <- all.vars(rlang::f_lhs(formula))[1]
    fn(data[[outcome]], ...)
  }
}

#' Label a vector
#'
#' @description
#' In conjunction with `brr::redirect_weights`, allows the use of `brr:brr`
#' fitters or other functions that have a different name for their weights
#' argument, and/or fitters that take a weights column name rather than the
#' actual weights.
#'
#' @param v vector
#' @param label name
#'
#' @return a vector with a label attribute
#' @export
label_vector <- function(v, label) {
  attr(v, "label") <- label
  v
}

#' Wrap a function that requires the name of the weights column so that instead
#' it accepts a vector of weights, and/or rename the weights argument.
#'
#' @description
#' This function makes it possible to modify fitters with unusual arguments
#' so they match the signature expected by `brr:brr`.
#'
#' @param f f
#' @param by_reference by_reference
#' @param name name
#'
#' @return a wrapped function
#' @export
#'
#' @examples
#' df <- data.frame(...)
#' weights <- label_vector(df$weights, "weights")
#'
#' lm_by_ref <- function(formula, data, weights) {
#'   weights <- data[[weights]]
#'   lm(formula, data, weights = weights)
#' }
#'
#' redirect_weights(lm_by_ref, by_reference = TRUE)(y ~ x, data = df, weights = weights)
#'
#' lm_wt <- function(formula, data, wt) {
#'   lm(formula, data, weights = wt)
#' }
#'
#' redirect_weights(lm_wt, name = "wt")(y ~ x, data = df, weights = weights)
redirect_weights <- function(f, by_reference = FALSE, name = "weights") {
  function(...) {
    # arguments as actual objects
    weights <- list(...)$weights
    # arguments as unevaluated expressions (important for formula environment)
    args <- as.list(substitute(list(...)))[-1L]

    if (by_reference) {
      wt <- attr(weights, "label", exact = TRUE)
    } else {
      wt <- weights
    }

    args[[as.symbol("weights")]] <- NULL
    args[[as.symbol(name)]] <- wt
    do.call(f, args)
  }
}
