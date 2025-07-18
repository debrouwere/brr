library("tidyverse")
library("vctrs")

list_glue <- function(l) {
  map_if(l, is_character, \(v) as.character(list_c(map(v, str_glue))))
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
#' See `pisa.rx.parquet` for a multiply imputed dataset that is compatible with
#' `brrl`, or alternatively use `pivot_longer(..., names_pattern = 'pv(\\d+)(math)', names_to = c('i', '.value'))`
#' on a wide format PISA dataset and merge these plausible values with imputations
#' such as those generated with `mice`.
#'
#' @param statistic statistical function or fitter with signature `fitter(data, weights)`
#' @param data plausible values, covariates and weights in long format, one imputation per row
#' @param conditions the replicates to be analyzed, e.g. `list(i = 1:5, weights = c("w_student_final", "w_student_r1")`
#' @param i number of imputations to to analyze, if only a subset of `data` should be used
#' @param r number of replications to perform, if only a subset of `weights` should be used
#' @param .progress show a progress bar
#'
#' @export
brrl <- function(statistic, data, conditions, i = NULL, r = NULL, .progress = TRUE) {
  conditions <- list_glue(conditions)

  r <- if (!is.null(r)) { r } else { length(conditions$weights) - 1 }
  weight_cols <- head(conditions$weights, n = r + 1)
  final_weight_cols <- weight_cols[1]
  replicate_weight_cols <- weight_cols[-1]

  w <- length(weight_cols)
  m <- if (!is.null(i)) { i } else { length(conditions$i) }
  s <- m * w

  if (.progress) {
    progressor <- cli::cli_progress_bar("Balanced repeated replication", total = s)
    tick <- invisibly(\() cli::cli_progress_update(id = progressor))
  } else {
    tick <- identity
  }

  tidy_statistic <- compose(tick, as_tidy, statistic)

  replicate <- function(imputation) {
    vec_rbind(!!!lapply(weight_cols, function(weight_col) {
      tidy_statistic(data = imputation, weights = imputation[[weight_col]])
    }), .names_to = "weights")
  }

  ixs_by_i <- vec_group_loc(data[["i"]])$loc |> head(n = m)
  imputations <- vec_chop(data, ixs_by_i)
  replications <- vec_rbind(!!!lapply(imputations, replicate), .names_to = "imputation")

  if (.progress) cli::cli_progress_done()

  new_brr(replications)
}
