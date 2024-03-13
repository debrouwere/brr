library("tidyverse")

as_factors <- function(df, colnames) {
  for (colname in colnames) {
    if (!is.factor(df[[colname]])) {
      old_type <- typeof(df[[colname]])
      message(str_glue("Converting {colname} from {old_type} to factor."))
      df[, colname] <- factor(df[[colname]])
    }
  }
  df
}

NA_LEVEL <- "<NA>"

na_value_to_level <- function(df, colnames) {
  for (colname in colnames) {
    df[, colname] <- fct_na_value_to_level(df[[colname]], level = NA_LEVEL)
  }
  df
}


#' Calculate weight adjustment factors so that the proportion of each strata in `df` matches the proportion in `df_target`
#'
#' @param df data frame for which to calculate the adjustment factors
#' @param df_target data frame that serves as the target or baseline
#' @param factors factor variables on which to poststratify
#' @param weights column name of a column that contains observation weights
#'
#' @export
poststrata <- function(df, df_target, factors, weights = NULL) {
  # Note that even though we are only computing adjustment factors that are supposed to be applied
  # to pre-existing weights, we must still know about those weights and cannot simply rely on
  # the unweighted `n` of observations. For example, if a country has, by proportion, twice as many
  # first-generation immigrants as the target country, but also happens to undersample them by a factor
  # of two relative to the target, the unweighted adjustment factor would be 1, when it ought to be 0.5.

  if (is.null(weights)) {
    count_wt <- count
  } else {
    count_wt <- partial(count, wt = .data[[weights]])
  }

  df <- df |>
    as_factors(factors) |>
    na_value_to_level(factors)
  df_target <- df_target |>
    as_factors(factors) |>
    na_value_to_level(factors)

  total <- count_wt(df)
  counts <- df |> count_wt(across({{ factors }}))
  counts$p <- counts$n / total$n

  target_total <- count_wt(df_target)
  target_counts <- df_target |> count_wt(across({{ factors }}))
  target_counts$p_target <- target_counts$n / target_total$n

  counts <- left_join(counts, target_counts[, c(factors, "p_target")], by = factors)

  # let's drop NA cells for now; in the future perhaps raise a warning
  # (this should not happen due to actual missing values, which should get their own cell,
  # but it can happen when there are fewer factor levels for the reference block than for
  # other blocks, e.g. a <NA> or other level in df but not df_target)
  counts <- drop_na(counts)
  counts$w <- counts$p_target / counts$p

  counts
}

#' Clip strata
#'
#' @description
#' also known (inaccurately) as trimming and truncating
#'
#' @param strata strata with weight adjustment factors, as produced by `poststrata`
#' @param upper maximum value to which to clip weight adjustment factors
#'
#' @export
clip_strata <- function(strata, upper = 5) {
  # (note that it would also be possible to clip the weights of actual
  # observations rather than the weight adjustment factors in the strata,
  # but that is not covered by this function)

  n <- strata$n
  w <- strata$w
  # adjust upper to account for renormalization
  total <- sum(n * w)
  excess <- sum(n * pmax(0, w - upper))
  upper <- (1 - excess / total) * upper
  # clip
  w_clipped <- pmin(w, upper)
  # renormalize
  strata$w <- w_clipped * sum(w * n) / sum(w_clipped * n)
  strata
}

#' Drop strata that include NA factor levels
#'
#' @param strata strata with weight adjustment factors, as produced by `poststrata`
#'
#' @export
drop_na_strata <- function(strata) {
  n <- strata$n
  w <- strata$w
  # clean
  factors <- strata |> select(where(is.factor))
  has_missing <- apply(sweep(factors, 2, "<NA>", "=="), 1, any)
  w_clean <- strata$w
  w[has_missing] <- 0.0
  strata[has_missing, "w"] <- 0.0
  # renormalize
  strata$w <- w * sum(w * n) / sum(w_clean * n)
  strata
}

#' Poststratify a dataset with strata weight adjustment factors
#'
#' @param df a data frame
#' @param strata strata with weight adjustment factors, as produced by `poststrata(df, ...)`
#' @param name rename the weight adjustment factors to `name`
#' @export
poststratify <- function(df, strata, name = "weight") {
  weight_col <- c("w")
  names(weight_col) <- name

  factors <- strata |>
    select(where(is.factor)) |>
    colnames()
  df <- df |> as_factors(factors)

  # <NA> is a stratum level, but when joining with the original data
  # it needs to be converted back to an actual NA value
  strata <- strata |>
    rename(all_of(weight_col)) |>
    mutate(across(where(is.factor), ~ na_if(.x, NA_LEVEL)))

  # two NA values are considered identical for the purposes of joining poststratification weights
  # (this is currently the default for `left_join`, added for clarity)
  left_join(df, strata[, c(factors, name)], by = factors, na_matches = "na")
}


# helper function for strata_collapse
total_distance <- function(distances, left_ixs, right_ixs) {
  ixs <- expand.grid(left = left_ixs, right = right_ixs)
  # as.vector performs columnwise flattening, which matches the ordering of expand.grid
  ixs$distance <- as.vector(distances[left_ixs, right_ixs])
  sum(ixs$distance)
}

distance_between_levels <- function(factor) {
  k <- length(levels(factor))
  coordinates <- attr(factor, "coordinates", exact = TRUE)

  if (!is.null(coordinates)) {
    d <- abs(matrix(coordinates, ncol = k, nrow = k, byrow = TRUE) - coordinates)
  } else if (is.ordered(factor)) {
    coordinates <- 1:k
    d <- abs(matrix(coordinates, ncol = k, nrow = k, byrow = TRUE) - coordinates)
  } else {
    d <- 1 - diag(1, nrow = k)
  }
  rownames(d) <- colnames(d) <- levels(factor)
  d
}

normalize <- function(x) {
  x / max(x)
}



#' Collapse small strata cells
#'
#' @description
#' For every stratification factor that is added, the amount of strata
#' multiplies by the k levels of that factor, so even with large samples it is
#' quite common to end up with a handful of combinations that match only a
#' handful of observations. This can lead to extreme weight adjustment factors,
#' which technically are not biased but nevertheless lead to very unstable
#' estimates. We can deal with extreme weights directly (by clipping to a
#' maximum of e.g. 5 or 10), indirectly by collapsing small cells until each
#' cell has at least n_min observations, or both. (You should also consider
#' dropping observations with NA categorizations for stratification factors.)
#'
#' The collapse algorithm below is fairly advanced in a couple of ways:
#'
#'   * the default strategy (distance) will look for opportunities to collapse small,
#'     similar cells across factor, rather than collapsing cells within one factor
#'     at a time as most algorithms do
#'   * it allows for any two levels within a factor to have a prespecified distance,
#'     a distance inferred from the ordering of an ordered factor, or an equal
#'     distance between all levels
#'
#' For a prespecified distance, the algorithm looks for a `coordinates`
#' attribute:
#'
#'   x <- factor(
#'     x=c(1, 1, 3, 2),
#'     levels=c('disagree', 'disagree slightly', 'agree nor disagree', 'agree')
#'     )
#'   attr(x, 'coordinates') <- c(1, 2, 3, 5)
#'
#' The two strategies that are implemented by `collapse_strata` are:
#'
#'   * the distance algorithm prefers to merge closely related cells regardless
#'     of the factor in which levels are to be collapsed
#'   * the adjacency algorithm prefers to merge cells by collapsing levels
#'     in the last factor (by column position) before moving to the next-to-last
#'     factor and so on; more or less as described in Little 1993
#'
#' (To test the adjacency strategy, which is still somewhat experimental,
#' specify strata in reverse order and see how that changes what cells get merged.)
#'
#' @param strata strata with weight adjustment factors, as produced by `poststrata(df, ...)`
#' @param n_min cells with fewer than `n_min` (weighted) observations are collapsed together with a nearby cell
#' @param strategy `distance` (recommended) or `adjacency` (classic)
#'
#' @export
collapse_strata <- function(strata, n_min = 10, strategy = "distance") {
  if (!(strategy %in% c("distance", "adjacency"))) {
    stop(str_glue("{strategy} is not a valid collapse_strata strategy, choose from: distance, adjacency"))
  }

  # `m` is the n of potentially merged strata (this way, the original stratum n
  # remains available after merging strata)
  strata$index <- 1:nrow(strata)
  strata$label <- as.character(strata$index)
  strata$m <- strata$n
  strata$pm <- 1.0

  # the distance between two cells is the amount of factors
  # for which they have a different level
  factors <- strata |> select(where(is.factor))
  # distances <- apply(factors, 1, function(factor) {
  #  rowSums(sweep(factors, 2, factor, '!='))
  # })
  within_distances <- map(colnames(factors), function(name) {
    normalize(distance_between_levels(factors[[name]]))
  })

  distances <- apply(factors, 1, function(left) {
    apply(factors, 1, function(right) {
      sum(map_dbl(1:ncol(factors), function(j) {
        within_distances[[j]][left[j], right[j]]
      }))
    })
  })

  # the adjacency between two cells is their distance, adjusted downwards
  # when any mismatch occurs in a factor near the front of the list and
  # adjusted upwards when the mismatch occurs in a factor near the end of the list
  nearness <- matrix(ncol(factors):1, ncol = ncol(factors), nrow = nrow(factors), byrow = TRUE)
  adjacencies <- apply(factors, 1, function(factor) {
    rowSums(sweep(factors, 2, factor, "!=") * nearness)
  })

  repeat {
    candidates_left <- strata |>
      filter(m < n_min) |>
      arrange(m)

    if (nrow(candidates_left) < 1) {
      break
    }

    left <- candidates_left |> first()
    left_factors <- left |>
      select(where(is.factor)) |>
      unlist()
    candidates_right <- strata
    candidates_right_factors <- candidates_right |> select(where(is.factor))

    # for merged strata, similarity is a fractional value: the average
    # of the distance of each left-right pair of cells weighted by their joint probability;
    # for example, with unit weights: (A,Z), (B,Z) <-> (A,X) is 1.5 dissimilar
    left_ixs <- strata |>
      filter(label == left$label) |>
      pull("index")
    candidate_ixs <- map(strata$label, function(label) {
      strata |>
        filter(label == {{ label }}) |>
        pull("index")
    })
    # broadcast weights to obtain joint probabilities, then do
    # an elementwise multiplication with the distance matrix
    weights <- matrix(strata$pm, nrow = nrow(strata), ncol = nrow(strata), byrow = TRUE) * strata$pm
    weighted_distances <- distances * weights
    candidates_right$dissimilarity <- map_dbl(candidate_ixs, function(right_ixs) {
      total_distance(weighted_distances, left_ixs, right_ixs)
    })

    if (strategy == "adjacency") {
      weighted_adjacencies <- adjacencies * weights
      candidates_right$adjacency <- map_dbl(candidate_ixs, function(right_ixs) {
        total_distance(weighted_adjacencies, left_ixs, right_ixs)
      })
    } else {
      candidates_right$adjacency <- 0
    }

    # we prefer to merge with other small strata to minimize the amount of merge operations,
    # though we have no strong theoretical argument for why a small number of merges
    # with a big impact on the resulting weights would be preferable over a larger number of
    # merges with a smaller impact; anyhow this only applies to ties and we will always
    # try to merge with the most similar or most adjacent cell first
    right <- candidates_right |>
      filter(label != left$label, dissimilarity >= 1) |>
      arrange(pick(adjacency, dissimilarity, m)) |>
      first()

    m_merged <- left$m + right$m
    w_merged <- (left$w * left$m + right$w * right$m) / m_merged

    labels <- c(left$label, right$label)
    ixs <- strata$label %in% labels
    strata[ixs, "m"] <- m_merged
    strata[ixs, "w"] <- w_merged
    strata[ixs, "pm"] <- strata[ixs, "n"] / strata[ixs, "m"]

    # by updating the label to account for the merge, we make sure that any
    # subsequent merges involve all rows involved in earlier merges
    strata[ixs, "label"] <- str_flatten(labels, ",")
  }

  strata
}
