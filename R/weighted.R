library("tidyverse")

READ_WEIGHTS <- c("w_read_fstuwt", str_glue("w_read_fstr{d}", d = 1:80))
MATH_WEIGHTS <- c("w_math_fstuwt", str_glue("w_math_fstr{d}", d = 1:80))
SCIE_WEIGHTS <- c("w_scie_fstuwt", str_glue("w_scie_fstr{d}", d = 1:80))
FINAL_WEIGHTS <- c("w_fstuwt")
REPLICATE_WEIGHTS <- str_glue("w_fstr{d}", d = 1:80)
WEIGHTS <- c(FINAL_WEIGHTS, REPLICATE_WEIGHTS)

#' Compute weighted summary statistics of data subsets (or the entire dataset) for one or more outcome variables
#'
#' @description
#' This function is powerful but can get a bit slow. To compute an aggregate for a single outcome
#' without any groups to condition on, consider using `brr::formulaic` instead, e.g.
#' `formulaic(weighted_mean)(outcome ~ 1, data, weights)` to create a formula interface for
#' a simple statistical function.
#'
#' @param formula formula
#' @param statistic statistical function or fitter with signature `fitter(data, weights, ...)`
#' @param data data
#' @param weights a vector of weights
#' @param ... additional arguments, passed onto `statistic``
#'
#' @export
#'
#' @examples
#' df <- tibble(x = 1:10, y = 11:20, w = 1 + 1:10 / 10, u = rep(1, 10), g = rep(c(1, 2), each = 5))
#' weighted_aggregate(x + y ~ g, weighted_mean, df, df$w)
#' # or more succinct...
#' weighted_mean_by <- partial(weighted_aggregate, statistic = weighted_mean)
#' weighted_mean_by(x + y ~ g, df, df$w)
weighted_aggregate <- function(formula, statistic, data, weights, ...) {
  outcomes <- all.vars(rlang::f_lhs(formula))
  groupers <- all.vars(rlang::f_rhs(formula))
  data <- bind_cols(data, .weights = weights)
  outcomes <- data |>
    group_by(across({{ groupers }})) |>
    summarize(across({{ outcomes }}, ~ statistic(.x, weights = .data$.weights, ...)))
  outcomes <- new_tibble(outcomes, class = "aggregate")
  attr(outcomes, "groups") <- groupers
  outcomes
}

#' Weighted mean conditional on zero or more groups
#'
#' @param formula formula
#' @param data data
#' @param weights a vector of weights
#' @param ... additional arguments, passed onto `weighted_mean`, such as `na_rm`
#'
#' @export
weighted_mean_by <- partial(weighted_aggregate, statistic = weighted_mean)

#' Weighted median conditional on zero or more groups
#'
#' @param formula formula
#' @param data data
#' @param weights a vector of weights
#' @param ... additional arguments, passed onto `weighted_median`, such as `na_rm`
#'
#' @export
weighted_median_by <- partial(weighted_aggregate, statistic = weighted_quantile, probs = 0.50)

#' Weighted mean with a formula interface
#'
#' @description
#' See `brr::formulaic` for more details. Useful because `brr::brr` expects a statistic
#' with a formula interface.
#'
#' @param formula formula of the form `outcome ~ 1`, with a single outcome and no left-hand side variables
#' @param data data
#' @param weights a vector of weights
#' @param ... additional arguments, passed onto `weighted_mean`, such as `na_rm`
#'
#' @export
pull_weighted_mean <- formulaic(weighted_mean)

#' Weighted median with a formula interface
#'
#' @description
#' See `brr::formulaic` for more details. Useful because `brr::brr` expects a statistic
#' with a formula interface.
#'
#' @param formula formula of the form `outcome ~ 1`, with a single outcome and no left-hand side variables
#' @param data data
#' @param weights a vector of weights
#' @param ... additional arguments, passed onto `weighted_median`, such as `na_rm`
#'
#' @export
pull_weighted_median <- formulaic(weighted_median)

#' Tidy an aggregate object
#'
#' @param results a tibble with aggregated outcomes as produced by `weighted_aggregate`
#' @param index the estimate to extract, if there is more than one outcome

#' @export
#' @importFrom generics tidy
tidy.aggregate <- function(results, index = 1) {
  groupers <- attr(results, "groups", exact = TRUE)

  if (!length(groupers)) {
    groupers <- c("term")
    results <- bind_cols(list(term = "(Intercept)"), results)
  }

  groups <- lapply(results[, groupers], str_replace_na, replacement = "<NA>")
  names(groups) <- NULL
  groups$sep <- ":"
  results$terms <- do.call(str_c, groups)

  tidied <- results[, c(ncol(results), length(groupers) + index)]
  colnames(tidied) <- c("term", "estimate")
  tidied
}

#' Weighted variance
#'
#' @description Authored by Gavin Simpson.
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param na_rm remove NA values
weighted_var <- function(x, weights, na_rm = FALSE) {
  if (na_rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na_rm)
}

#' Weighted standard deviation
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param na_rm remove NA values
#'
#' @export
weighted_sd <- function(x, weights, na_rm = FALSE) {
  sqrt(weighted_var(x = x, weights = weights, na_rm = na_rm))
}

#' Weighted scaling and centering of a vector
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param na_rm in line with `base::scale`, removes NAs by default
#'
#' @export
weighted_scale <- function(x, weights, na_rm = TRUE) {
  mu <- weighted_mean(x, weights, na_rm = na_rm)
  sigma <- weighted_sd(x, weights, na_rm = na_rm)
  (x - mu) / sigma
}

#' Scaling and centering of a vector by the weighted mean and sd of a target vector
#'
#' @description
#' Sometimes this is useful to make sure that a bunch of different vectors all receive the same rescaling
#'
#' @param x a data vector
#' @param x_target a target data vector
#' @param weights a target weights vector
#' @param na_rm in line with `base::scale`, removes NAs by default
#'
#' @export
weighted_scale_by <- function(x, x_target, weights, na_rm = TRUE) {
  mu <- weighted_mean(x_target, weights, na_rm = na_rm)
  sigma <- weighted_sd(x_target, weights, na_rm = na_rm)
  (x - mu) / sigma
}

#' Weighted mean
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param na_rm na_rm
#' @param ...
#'
#' @export
weighted_mean <- function(x, weights, na_rm = FALSE, ...) {
  # weighted.mean chokes on NA weights even with na.rm=TRUE
  if (na_rm) {
    weights <- replace_na(weights, 0)
  }
  weighted.mean(x, weights, na.rm = na_rm, ...)
}

#' A weighted sum
#'
#' @description
#' Really just a sum of products. Weights are not normalized but used as-is.
#'
#' @param data a data vector
#' @param weights a weights vector
#' @param na_rm na_rm
#'
#' @export
weighted_sum <- function(x, weights, na_rm = FALSE) {
  sum(x * weights, na.rm = na_rm)
}

weighted_quantile_generic <- function(x, probs, cdf.gen, na.rm = FALSE, weights = NA) {
  if (na.rm) {
    weights <- weights[i <- !is.na(x)]
    x <- x[i]
  }

  n <- length(x)
  if (any(is.na(weights))) {
    weights <- rep(1 / n, n)
  }
  # Kish's effective sample size
  nw <- sum(weights)^2 / sum(weights^2)

  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]

  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))

  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}

#' Distribution-free weighted quantile estimator
#'
#' @description
#' Authored by Andrey Akinshin according to the proposal in Harrell & Davis 1982
#
# "Sample quantiles have many desirable properties. However, they also have drawbacks.
# They are not particularly efficient estimators of location for distributions such as
# the normal, good estimators of the variance of sample quantiles do not exist for general
# distributions, sample quantiles may not be jackknifed, and the sample median differs
# in form and in efficiency depending on the sample size being even or odd."
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param probs requested quantiles
#' @param na_rm na_rm
#'
#' @export
weighted_quantile_hd82 <- function(x, weights, probs, na_rm = FALSE) {
  cdf.gen <- function(n, p) {
    return(function(cdf.probs) {
      pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
    })
  }
  weighted_quantile_generic(x, probs, cdf.gen, na.rm = na_rm, weights = weights)
}

#' Naive weighted sample quantile estimator
#'
#' @description
#' Weighted sample quantile estimators Types 1-3 according to the Hyndman-Fan
#' typology.
#'
#' Sample quantiles may not always be the most efficient, and a Type 1 estimator
#' is not even unbiased, but given the large sample sizes of PISA and the BRR
#' methodology (in part designed to allow for naive quantile computations),
#' these are not relevant concerns and it is much faster than the Harell-Davis
#' estimator.
#'
#' Type 2 and 3 implementations are approximate algorithms that do not
#' completely match the behavior of `stats::quantile`:
#'
#'   * type 1 is the inverse cdf
#'   * type 2 takes the average even when there is no tie
#'   * type 3 takes the nearest quantile instead of the nearest even
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param prob requested quantiles
#' @param type Hyndman-Fan quantile type
#' @param na_rm na_rm
#'
#' @export
weighted_quantile <- function(x, weights, prob, type = 1, na_rm = FALSE) {
  if (na_rm) {
    i <- !is.na(x)
    weights <- weights[i]
    x <- x[i]
  }

  ixs <- order(x)
  cdf <- cumsum(weights[ixs]) / sum(weights)

  i <- which.max(cdf > prob)
  ix <- ixs[i]
  ix0 <- ixs[i - 1]

  switch(type,
    {
      x[ix0]
    },
    {
      (x[ix0] + x[ix]) / 2
    },
    {
      ifelse(cdf[ix] - cdf[ix0] >= 2 * prob, x[ix0], x[ix])
    }
  )
}

weighted_median <- partial(weighted_quantile, prob = 0.5)


#' Weighted variance explained
#'
#' @param formula formula
#' @param data data
#' @param weights weights
#' @param ...
#'
#' @export
weighted_r2 <- function(formula, data, weights, ...) {
  fit <- lm(formula, data = data, weights = weights)
  summary(fit)$r.squared
}


#' Convert from reliability (complement of the error proportion of the variance) to the error variance
#'
#' @description
#' In psychometrics, reliability is commonly expressed as 1 minus the variance of the error to the total variance of the variable,
#' but software packages typically expect an error variance instead.
#'
#' @param x x
#' @param r r
#' @param ...
#'
#' @export
reliability_to_variance <- function(x, r, ...) {
  var(x, ...) * (1 - r)
}

#' Convert from weighted reliability (complement of the error proportion of the variance) to the weighted error variance
#'
#' @description
#' In psychometrics, reliability is commonly expressed as 1 minus the variance of the error to the total variance of the variable,
#' but software packages typically expect an error variance instead.
#'
#' @param x a data vector
#' @param weights a weights vector
#' @param r reliability between 0.0 and 1.0
#' @param ...
#'
#' @export
weighted_reliability_to_variance <- function(x, weights, r, ...) {
  weighted_var(x, w, ...) * (1 - r)
}


#' Weighted covariance matrix
#'
#' @param data data
#' @param weights a non-negative and non-zero vector of weights for each observation. Its length must equal the number of rows of x
#' @param center either a logical or a numeric vector specifying the centers to be used when computing covariances
#' @param method string specifying how the result is scaled, `unbiased` or `ML`
#'
#' @export
weighted_cov_matrix <- function(data, weights, center = TRUE, method = "unbiased", na_rm = FALSE) {
  if (na_rm) {
    mask <- vctrs::vec_detect_complete(data)
    data <- data[mask, ]
    weights <- weights[mask]
  }
  cov.wt(data, wt = weights, center = center, method = method)$cov
}

#' Weighted correlation matrix
#'
#' @param data data
#' @param weights a non-negative and non-zero vector of weights for each observation. Its length must equal the number of rows of x
#' @param center either a logical or a numeric vector specifying the centers to be used when computing correlations
#' @param method string specifying how the result is scaled, `unbiased` or `ML`
#'
#' @export
weighted_cor_matrix <- function(data, weights, center = TRUE, method = "unbiased", na_rm = FALSE) {
  if (na_rm) {
    mask <- vctrs::vec_detect_complete(data)
    data <- data[mask, ]
    weights <- weights[mask]
  }
  cov.wt(data, wt = weights, cor = TRUE, center = center, method = method)$cor
}

# `ginv` and `marginal_to_partial` are inspired by code in William Revelle's `psych package`
ginv <- function(x, tol = sqrt(.Machine$double.eps)) {
  decomposition <- svd(x)
  d <- decomposition$d
  u <- decomposition$u
  v <- decomposition$v
  nonzero <- d > max(tol * decomposition$d[1], 0)
  v[, nonzero, drop = FALSE] %*% (1 / d[nonzero] * t(u[, nonzero, drop = FALSE]))
}

marginal_to_partial <- function(x) {
  inverse <- ginv(x)
  smc <- 1 - 1 / diag(inverse)
  residual <- -inverse
  diag(residual) <- 1 / (1 - smc)
  residual <- cov2cor(residual)
  rownames(residual) <- colnames(residual) <- colnames(x)
  residual
}

# `na_rm` is equivalent to `cor(use = "complete.obs")`
# (pairwise complete etc. doesn't make sense for a partial correlation matrix)
partial_cor_matrix <- function(data, center = TRUE, method = "unbiased", na_rm = FALSE) {
  weights <- rep(1, nrow(data))
  m <- weighted_cor_matrix(data, weights = weights, center = center, method = method, na_rm = na_rm)
  marginal_to_partial(m)
}

weighted_partial_cor_matrix <- function(data, weights, center = TRUE, method = "unbiased", na_rm = FALSE) {
  m <- weighted_cor_matrix(data, weights = weights, center = center, method = method, na_rm = na_rm)
  marginal_to_partial(m)
}

#' Weighted covariance
#'
#' @param x a numeric vector, matrix or data frame
#' @param y a numeric vector, matrix or data frame
#' @param weights a vector of weights for each observation
#' @param na_rm remove observations for which either the weight or one or both variables are NA
#' @param scale calculate a correlation instead
#'
#' @export
weighted_cov <- function(x, y, weights, na_rm = FALSE, scale = FALSE) {
  data <- tibble(x = x, y = y, .weights = weights)
  if (na_rm) data <- drop_na(data)
  xy <- select(data, -.weights)
  w <- pull(data, .weights)
  weighted_cov_matrix(xy, weights = w)[1,2]
}

#' Weighted correlation
#'
#' @param x a numeric vector, matrix or data frame
#' @param y a numeric vector, matrix or data frame
#' @param weights a vector of weights for each observation
#' @param na_rm remove observations for which either the weight or one or both variables are NA
#'
#' @export
weighted_cor <- function(x, y, weights, na_rm = FALSE) {
  data <- tibble(x = x, y = y, .weights = weights)
  if (na_rm) data <- drop_na(data)
  xy <- select(data, -.weights)
  w <- pull(data, .weights)
  weighted_cor_matrix(xy, weights = w)[1,2]
}
