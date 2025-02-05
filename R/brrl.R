library("tidyverse")

list_glue <- function(l) {
  map_if(l, is_character, \(v) as.character(list_c(map(v, str_glue))))
}

expand_brr <- function(conditions) {
  glued <- list_glue(conditions)
  expand_grid(!!!glued)
}

pivot_brr <- function(data, outcomes) {
  data |> pivot_longer(
    all_of(outcomes),
    names_to = 'i',
    values_to = common_suffix(outcomes)
  )
}

#' Indicate that a data frame contains balanced repeated replications
#'
#' @param x a data frame or tibble with fits, as produced by `brr::brr`
new_brr <- function(fits) {
  structure(fits, class = c("brr", "tbl_df", "tbl", "data.frame"))
}

#' Balanced repeated replications
#'
#' @description
#' `brrl` allows multiple imputation of covariates (not just plausible values)
#' without performance overhead relative to the analysis of a single imputation
#' as with the original `brr` function.
#'
#' Because `weights` are shared by all imputations, they are typically provided
#' by a data frame with 1/5th or 1/10th the amount of rows as the `data` set.
#'
#' See `pisa.rx.parquet` for a multiply imputed dataset that is compatible with
#' `brrl`, or alternatively use `pivot_longer(..., names_pattern = 'pv(\\d+)(math)', names_to = c('i', '.value'))`
#' on a wide format PISA dataset and merge these plausible values with imputations
#' such as those generated with `mice`.
#'
#' @param statistic statistical function or fitter with signature `fitter(data, weights)`
#' @param data plausible values and covariates in long format, one imputation per row
#' @param weights weights shared by all imputations, of which the first weight will be used as final weight
#' @param i number of imputations to to analyze, if only a subset of `data` should be used
#' @param r number of replications to perform, if only a subset of `weights` should be used
#' @param .progress show a progress bar
#' @param ...
#'
#' @export
brrl <- function(statistic, data, weights, conditions, i = NULL, r = NULL, .progress = TRUE, .verbose = TRUE) {
  conditions <- list_glue(conditions)

  r <- if (!is.null(r)) { r } else { length(conditions$weights) - 1 }
  weight_cols <- head(conditions$weights, n = r + 1)
  final_weight_cols <- weight_cols[1]
  replicate_weight_cols <- weight_cols[-1]
  w <- length(weight_cols)

  m <- if (!is.null(i)) {
    i
  } else if (!is.null(conditions$outcome)) {
    length(conditions$outcome)
  } else {
    length(conditions$i)
  }

  if (!is.null(conditions$outcome)) {
    data <- pivot_brr(data, conditions$outcome)
  }

  s <- m * w

  if (.progress) {
    progressor <- cli::cli_progress_bar("Balanced repeated replication", total = s)
    tick <- invisibly(\() cli::cli_progress_update(id = progressor))
  } else {
    tick <- identity
  }

  tidy_statistic <- compose(tick, as_tidy, statistic)

  replicate <- function(ixs, ...) {
    list_rbind(map(weight_cols, function(weight_col) {
      tidy_statistic(data = slice(data, ixs), weights = weights[[weight_col]])
    }), names_to = 'weights')
  }

  # NOTE: although this approach lends itself to parallelization, e.g. using
  # `furrr::future_map` on each imputation, so much data is transferred
  # (needlessly but unavoidably) that the performance gains are barely worth it
  ixs_by_i <- data |> group_by(i) |> group_rows() |> head(n = m)
  replications <- map(ixs_by_i, replicate) |>
    list_rbind(names_to = 'imputation')

  if (.progress) cli::cli_progress_done()

  new_brr(replications)
}
