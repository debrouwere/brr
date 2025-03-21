library("tidyverse")
library("tidyselect")
library("rlang")
library("cli")

#' Indicate that a data frame contains balanced repeated replications
#'
#' @param x a data frame or tibble with fits, as produced by `brr::brr`
new_brr <- function(fits) {
  structure(fits, class = c("brr", "tbl_df", "tbl", "data.frame"))
}

#' Balanced repeated replications
#'
#' @param formula a formula or a vector of outcomes
#' @param statistic statistical function or fitter with signature `fitter(data, weights, ...)`
#' @param data data
#' @param final_weights a tidy selection or a vector of column names
#' @param replicate_weights a tidy selection or a vector of column names
#' @param r number of replications to perform (for convenience, you can also just pass fewer columns to `replicate_weights`)
#' @param .progress show a progress bar
#' @param ...
#'
#' @export
brr <- brrw <- function(formula, statistic, data, final_weights, replicate_weights, r = 80, .progress = TRUE, ...) {
  is_marginal <- !rlang::is_formula(formula)
  
  # fast path for simple aggregates of the entire dataset
  if (is_marginal) {
    outcomes <- formula
    predictors <- NA
  } else {
    outcomes <- all.vars(rlang::f_lhs(formula))
    predictors <- rlang::f_text(formula)
  }
  
  # TODO: check whether outcomes are present in the dataset
  
  # select weights from `data` using selection helpers such as `starts_with`,
  # `matches` etc. from the tidyselect package, using nonstandard evaluation
  final_weights_quosure <- rlang::enquo(final_weights)
  final_weights <- data |> select({{ final_weights_quosure }})
  replicate_weights_quosure <- rlang::enquo(replicate_weights)
  replicate_weights <- data |> select({{ replicate_weights_quosure }})
  weights <- bind_cols(final_weights, replicate_weights[, 1:r])
  
  # harmonization of wide and long format for plausible values
  # FIXME: \d+ is too brittle, would also match pv1read1, st01q01 etc.
  imputations <- as.integer(str_extract(unique(outcomes), "\\d+"))
  
  conditions <- tibble(outcome = outcomes, imputation = imputations) |>
    expand_grid(weights = colnames(weights))
  conditions$is_final <- conditions$weights == colnames(final_weights)[1]
  conditions$formula <- if (!is_marginal) {
    formulae <- str_c(outcomes, " ~ ", predictors)
    names(formulae) <- outcomes
    formulae[conditions$outcome]
  } else {
    NA
  }
  
  if (.progress) {
    progressor <- cli::cli_progress_bar("Balanced repeated replication", total = nrow(conditions))
    tick <- invisibly(\() cli::cli_progress_update(id = progressor))
  } else {
    tick <- identity
  }
  
  # TODO: can we support long and wide format without having these parallel code
  # paths all over the place? Convert internally or something like that?
  #
  # I'm also wondering whether `x` and `w`, although convenient, don't lead
  # to an unnecessary data copy?
  replicate <- if (is_marginal) {
    function(condition) {
      x <- data |> pull(condition$outcome)
      w <- weights |> pull(condition$weights)
      statistic(x = x, weights = w, ...)
    }
  } else {
    function(condition) {
      f <- as.formula(condition$formula)
      x <- data
      w <- label_vector(weights[[condition$weights]], condition$weights)
      statistic(formula = f, data = x, weights = w, ...)
    }
  }
  
  replicate_tidily <- compose(tick, as_tidy, replicate)
  
  replications <- conditions |>
    rowwise() |>
    reframe(
      outcome = outcome,
      formula = formula,
      weights = weights,
      imputation = imputation,
      is_final = is_final,
      results = replicate_tidily(.data)
    ) |>
    unnest_wider(results)
  
  if (.progress) cli::cli_progress_done()
  
  # t0: W_FSTUWT are the final weights, used to compute point estimates and imputation variance
  #     (imputation variance arises due to matrix sampling)
  # t:  W_FSTR* are the replicate weights, used to compute the estimation variance
  #     (the estimates themselves are then discarded)
  new_brr(replications)
}