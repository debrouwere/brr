library("dplyr")
library("tidyr")
library("broom")
library("broom.mixed")
library("cli")
library("rlang")

is_tidy <- function(data) {
  is_tibble(data) & all(c("term", "estimate") %in% colnames(data))
}

as_tidy <- function(object) {
  if (is_tidy(object)) {
    object
  } else if (is.numeric(object)) {
    tibble(
      term = "(Intercept)",
      estimate = object
    )
  } else {
    broom::tidy(object)
  }
}

#' Indicate that a data frame contains balanced repeated replications
#'
#' @param x a data frame or tibble with fits, as produced by `brr::brr`
#'
#' @export
new_brr <- function(fits) {
  fits |>
    select(imputation, weights, term, estimate) |>
    structure(class = c("brr", "tbl_df", "tbl", "data.frame"))
}

#' Indicate that a data frame contains balanced repeated replications
#'
#' @description
#' Alias to `brr::new_brr`
#'
#' @param x a data frame or tibble with fits, as produced by `brr::brr`
#'
#' @export
as_brr <- function(replications) {
  new_brr(replications)
}

#' Test whether an object is a collection of replications
#'
#' @param x any object
#'
#' @export
is_brr <- function(x) {
  "brr" %in% class(x)
}

#' Extract model coefficients from the balanced repeated replications of a model fit
#'
#' @param replications BRR replications
#' @param final_weights the weights index that represents the final weights (usually this is 1)
#' @param simplify return a named vector instead of `tibble(term, estimate)`
#'
#' @export
coef.brr <- function(replications, final_weights = 1, simplify = TRUE) {
  final <- replications |> filter(weights == {{ final_weights }})

  coefs <- final |>
    summarize(estimate = mean(estimate), .by = term)

  # simplification is the default, to match the format of `brr.lm`
  if (simplify) {
    coefs |>
      pivot_wider(names_from = "term", values_from = "estimate") |>
      unlist()
  } else {
    coefs
  }
}


#' Tidy a BRR object (experimental)
#'
#' @param replications A brr object, as produced by `brr::brr`
#'
#' @export
#' @importFrom generics tidy
tidy.brr <- function(replications) {
  coef.brr(replications, simplify = FALSE)
}

#' Imputation, estimation and total variance of balanced repeated replications
#'
#' @param replications BRR replications
#' @param perturbation perturbation
#' @param final_weights the weights index that represents the final weights (usually this is 1)
#'
#' @return a list of mean and variance component vectors
#' @export
VarCorr.brr <- function(replications, perturbation = 0.50, final_weights = 1) {
  t0 <- replications |> filter(weights == {{ final_weights }})
  t <- replications |> filter(weights != {{ final_weights }})

  # number of imputations, replications, effective replications and the replication design effect
  # for variance calculations using Fay's method
  m <- length(unique(t$imputation))
  w <- length(unique(t$weights))
  w_eff <- w * (1 - perturbation)^2
  w_deff <- w / w_eff

  # means by plausible value
  m0 <- t0 |>
    summarize(mean = mean(estimate), .by = c(term, imputation))

  # final means
  mm0 <- t0 |>
    summarize(mean = mean(estimate), .by = term)

  t0 <- full_join(t0, mm0, by = c("term")) |>
    mutate(squared_deviation = (estimate - mean)^2)
  t <- full_join(t, m0, by = c("term", "imputation")) |>
    mutate(squared_deviation = (estimate - mean)^2)

  # see the technical manual with regards to the factor `1 + (1 / m)`
  imputation_variance <- t0 |>
    summarize(imputation = sum(squared_deviation) * (1 + 1 / m) / (m - 1), .by = term)

  estimation_variance <- t |>
    summarize(estimation = mean(squared_deviation) * w_deff, .by = term)

  full_join(imputation_variance, estimation_variance, by = "term")
}

#' Subtract sets of replications from each other
#'
#' @param e1 lefthand set of replications
#' @param e2 righthand set of replications
#'
#' @description Useful to compare scores between countries or between cycles when passed on to
#'   `confint`
#'
#' @export
`-.brr` <- function(e1, e2) {
  e <- inner_join(e1, e2, by = c("imputation", "weights", "term"))
  mutate(e, estimate = estimate.x - estimate.y, .keep = "unused")
}

vc <- function(varcorr_tbl) {
  select(varcorr_tbl, -term)
}

interval <- function(estimate, variance, level) {
  se <- sqrt(variance)
  critical <- qnorm(1 - (1 - level) / 2)
  tibble(
    se = se,
    lower = estimate - critical * se,
    upper = estimate + critical * se,
  )
}

#' Confidence intervals for model parameters of the balanced repeated model fit replications
#'
#' @param replications one or more sets of BRR replications
#' @param level the confidence level required
#' @param links a vector of link errors, usually zero for the reference set and positive for all
#'   others
#' @param reference a function that selects the reference set from among the sets of replications
#' @param simplify if only a single set of replications is passed, do not wrap its intervals in a list
#' @param extra include variance components in the output
#'
#' @description If `extra = TRUE`, output includes a confidence interval for each variance
#'   component. These intervals are cumulative: they include the uncertainty related to all previous
#'   components as well.
#'
#'   If multiple sets of replications are passed to `confint`, confidence intervals for all sets
#'   except for the last one (or the one picked by the function provided to the `reference`
#'   argument) will include the error of both current set and reference set. This allows for easy
#'   calculation of the confidence intervals for comparisons between countries or between cycles.
#'
#'   To use this functionality for comparisons between countries within a cycle, the `links`
#'   argument must be explicitly set to 0. (An alternative would be to subtract the replication
#'   estimates from one set from the other, which can be done using the arithmetic operator `-`, and
#'   to ask for a confidence interval for the resulting set using `confint(reps2 - reps1)`. The
#'   results that this produces are potentially more intuitive.)
#'
#'   For comparisons between cycles, the link error for each comparison between the reference set
#'   and subsequent sets must be provided by the user. This information is available as part of the
#'   official PISA reports.
#'
#'   Passing multiple sets of replications without setting the `links` argument will result in `NA`
#'   confidence intervals.
#'
#' @export
confint.brr <- function(..., level = 0.95, links = NA, reference = last, simplify = TRUE, extra = FALSE) {
  replications <- list2(...)

  # there can be no link error within a single set of replications (unless explicitly demanded)
  if(is.na(links) & length(replications) == 1) links <- 0.0

  if (length(replications) != length(links) & length(links) != 1) {
    cli_abort("{length(replications)} replications but only {length(links)} links")
  }

  variances <- map(replications, VarCorr.brr)
  comparisons <- map2(variances, links, function(variance, link) {
    bind_cols(
      term = variance$term,
      vc(variance) + as.integer(link != 0.0) * vc(reference(variances)),
      link = link^2,
    )
  })
  means <- map(replications, \(rr) unname(coef.brr(rr)))
  margins <- map2(comparisons, means, function(variance, means) {
    variance |>
      as_tibble() |>
      mutate(
        estimate = means,
        link = interval(estimate, imputation + estimation + link, level),
        estimation = interval(estimate, imputation + estimation, level),
        imputation = interval(estimate, imputation, level)
      ) |>
      select(term, estimate, imputation, estimation, link)
  })

  if (!extra) {
    margins <- map(margins, function(m) {
      m |>
        select(term, estimate, link) |>
        unnest(link)
    })
  }

  if (simplify & length(replications) == 1) {
    margins[[1]]
  } else {
    margins
  }
}

#' Summarize balanced repeated replications
#'
#' @description
#' This function is a placeholder and has not yet been implemented.
#' In the meanwhile, use `coef` and `confint`.
#'
#' @param replications
#'
#' @export
summary.brr <- function(replications) {
  # * confint(verbose=TRUE)
  # * brr_n
  # * some diagnostics (single outcome? NA values?)
  # * we can't do that much more because people could be using any underlying model
  #   so we can't really say anything about R2 or effective sample size or ...
}
