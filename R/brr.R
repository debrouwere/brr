library("tidyverse")
library("tidyselect")
library("rlang")
library("broom")
library("cli")

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

#' Balanced repeated replications
#'
#' @param formula formula
#' @param statistic statistical function or fitter with signature `fitter(data, weights, ...)`
#' @param data data
#' @param final_weights depending on `.select` either a tidy selection (the default);
#'   a vector of column names; or a vector, one-column data frame or one-column tibble
#' @param replicate_weights depending on `.select` either a tidy selection (the default);
#'   a vector of column names; or a data frame or tibble
#' @param replications number of replications to perform (for convenience, you can also just pass fewer columns to `replicate_weights`)
#' @param .progress show a progress bar
#' @param .select select weights from `data` using selection helpers such as `starts_with`, `matches` etc. from the tidyselect package (tidy); using column names (names); or by passing the data itself (noop)
#' @param ...
#'
#' @export
brr <- function(formula, statistic, data, final_weights, replicate_weights, r = 80, .progress = TRUE, .select = "tidy", ...) {
  outcomes <- all.vars(rlang::f_lhs(formula))
  predictors <- rlang::f_text(formula)

  # select weights from `data` using selection helpers such as `starts_with`,
  # `matches` etc. from the tidyselect package, using nonstandard evaluation
  switch(.select,
    tidy = {
      final_weights_selectors <- rlang::enquo(final_weights)
      final_ixs <- eval_select(final_weights_selectors, data = data)
      final_weights <- data[, final_ixs]

      replicate_weights_selectors <- rlang::enquo(replicate_weights)
      replicate_ixs <- eval_select(replicate_weights_selectors, data = data)
      replicate_weights <- data[, replicate_ixs]
    },
    names = {
      final_weights <- data[, final_weights]
      replicate_weights <- data[, replicate_weights]
    },
    {
      # by default, we take `final_weights` and `replicate_weights` as-is,
      # which should be in the form of data frames or tibbles
    }
  )

  weights <- bind_cols(final_weights, replicate_weights[, 1:r])
  final_weights_colname <- colnames(final_weights)[1]

  conditions <- expand_grid(
    outcome = outcomes,
    weights = colnames(weights)
  )

  formulae <- str_c(outcomes, " ~ ", predictors)
  names(formulae) <- outcomes
  conditions$formula <- formulae[conditions$outcome]

  progressor <- cli_progress_bar("Balanced repeated replication", total = nrow(conditions))

  replicate <- function(condition) {
    if (.progress) cli::cli_progress_update(id = progressor)
    as_tidy(statistic(
      formula = as.formula(condition$formula),
      data = data,
      weights = label_vector(weights[[condition$weights]], condition$weights),
      ...
    ))
  }

  replications <- conditions |>
    rowwise() |>
    reframe(
      outcome = outcome,
      weights = weights,
      formula = formula,
      results = replicate(.data)
    ) |>
    unnest_wider(results)

  cli::cli_progress_done()


  # t0: W_FSTUWT are the final weights, used to compute point estimates and imputation variance
  #     (imputation variance arises due to IRT scaling, where we assume proficiency is latent)
  # t:  W_FSTR* are the replicate weights, used to compute the estimation variance
  #     (the estimates themselves are then discarded)

  structure(list(
    t0 = replications |> filter(weights == {{ final_weights_colname }}),
    t  = replications |> filter(weights != {{ final_weights_colname }})
  ), class = "brr")
}

#' Diagnose common problems with BRR replications
#'
#' @param replications
#'
#' @export
brr_diagnose <- function(replications) {
  na_outcomes <- c()
  na_weights <- c()
  map(replications$t, function(tx) {
    data <- tx |> drop_na()
    na_outcomes <- setdiff(unique(data$outcome), unique(data$outcome))
    na_weights <- setdiff(unique(data$weights), unique(data$weights))
    if (length(na_outcomes)) message("Missing replication outcomes: ", str_flatten_comma(na_outcomes))
    if (length(na_weights)) message("Missing replication weights: ", str_flatten_comma(na_weights))
    data
  })
}

#' Extract model coefficients from the balanced repeated replications of a model fit
#'
#' @description
#' The `na_rm` argument should be used sparingly because it may hide model or fit errors,
#' but it can be particularly useful for repeated cross-sectional analyses of PISA data
#' from 2000-now, allowing you to specify a model with 10 plausible values even though
#' 2000, 2003, 2006, 2009 and 2012 assessments only include 5.
#'
#' @param replications BRR replications
#' @param na_rm remove replications with NA estimates
#' @param simplify return a named vector instead of `tibble(term, estimate)`
#'
#' @export
coef.brr <- function(replications, na_rm = FALSE, simplify = TRUE) {
  if (na_rm) {
    results <- replications$t0 |> drop_na(estimate)
  } else {
    results <- replications$t0
  }

  coefs <- results |>
    group_by(term) |>
    summarize(
      estimate = mean(estimate, na.rm = na_rm)
    )

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
tidy.brr <- function(replications, na_rm = FALSE) {
  coef.brr(replications, na_rm = na_rm, simplify = FALSE)
}

brr_n <- function(t, perturbation = 0.50) {
  n <- list()
  n$outcomes <- length(unique(t$outcome))
  n$replications <- length(unique(t$weights))
  # denominator for variance calculations using Fay's method
  n$effective_replications <- n$replications * (1 - perturbation)^2
  n$brr_design_effect <- n$replications / n$effective_replications
  n
}

#' Imputation, estimation and total variance of balanced repeated replications
#'
#' @description
#' The `na_rm` argument should be used sparingly because it may hide model or fit errors,
#' but it can be particularly useful for repeated cross-sectional analyses of PISA data
#' from 2000-now, allowing you to specify a model with 10 plausible values even though
#' 2000, 2003, 2006, 2009 and 2012 assessments only include 5.
#'
#' @param replications BRR replications
#' @param perturbation perturbation
#' @param imputation Incorporate imputation variance due to plausible values. TRUE by default but can be disabled when there's only a single outcome.
#' @param na_rm na_rm
#'
#' @return A list of mean and variance component vectors.
#' @export
brr_var <- function(replications, perturbation = 0.50, imputation = TRUE, na_rm = FALSE) {
  if (na_rm) {
    t0 <- replications$t0 |> drop_na(estimate)
    t <- replications$t |> drop_na(estimate)
  } else {
    t0 <- replications$t0
    t <- replications$t
  }

  n <- brr_n(t, perturbation)

  # means by outcome
  m0 <- t0 |>
    group_by(term, outcome) |>
    summarize(mean = mean(estimate)) |>
    ungroup()

  # final means
  mm0 <- t0 |>
    group_by(term) |>
    summarize(mean = mean(estimate))

  t0 <- full_join(t0, mm0, by = c("term")) |>
    mutate(squared_deviation = (estimate - mean)^2)
  t <- full_join(t, m0, by = c("term", "outcome")) |>
    mutate(squared_deviation = (estimate - mean)^2)

  # FIXME: n$outcomes may well be inaccurate if na_rm is selected!
  imputation_variance <- t0 |>
    group_by(term) |>
    summarize(imputation = sum(squared_deviation) / (n$outcomes - 1))

  estimation_variance <- t |>
    group_by(term) |>
    summarize(estimation = mean(squared_deviation) * n$brr_design_effect)

  # note that first we obtain the imputation variance by dividing by n-1,
  # but then here we multiply by n+(1/n), which gets you almost but not
  # quite back to 1; why PISA recommends this approach eludes me but it is
  # what it says in the technical manual
  variance <- full_join(imputation_variance, estimation_variance, by = "term")
  variance$total <- variance$estimation + (1 + 1 / n$outcomes) * variance$imputation
  variance
}

#' Imputation, estimation and total standard deviation of balanced repeated replications
#'
#' @param replications BRR replications
#' @param perturbation perturbation
#' @param imputation Incorporate imputation variance due to plausible values. TRUE by default but can be disabled when there's only a single outcome.
#'
#' @return A list of mean and variance component vectors.
#' @export
brr_sd <- function(replications, perturbation = 0.50, imputation = TRUE, na_rm = FALSE) {
  s2 <- brr_var(replications, perturbation, imputation, na_rm)
  s <- s2 |> mutate(across(!term, ~ sqrt(.x)))
  s
}

#' Confidence intervals for model parameters of the balanced repeated model fit replications
#'
#' @param replications BRR replications
#' @param level the confidence level required
#' @param perturbation perturbation
#' @param imputation Incorporate imputation variance due to plausible values. TRUE by default but can be disabled when there's only a single outcome.
#' @param extra include standard error and variance components in the output

#' @export
confint.brr <- function(replications, level = 0.95, perturbation = 0.50, imputation = TRUE, extra = FALSE, na_rm = FALSE) {
  # means <- replications$t0 |> select(where(is.numeric)) |> summarize(across(everything(), mean))
  variances <- brr_var(replications, perturbation = perturbation, imputation = imputation, na_rm = na_rm)
  if (imputation) {
    variance <- variances$total
  } else {
    variance <- variance$estimation
  }
  means <- unname(coef.brr(replications, na_rm = na_rm))
  critical <- qnorm(1 - (1 - level) / 2)
  margins <- sqrt(variance) * critical
  names <- variances$term

  base <- tibble(
    term = names,
    estimate = means,
    lower = means - margins,
    upper = means + margins
  )

  additional <- tibble(
    imputation_var = variances$imputation,
    estimation_var = variances$estimation,
    total_var = variances$total,
    se = sqrt(variances$total)
  )

  if (extra) {
    bind_cols(base, additional)
  } else {
    base
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
