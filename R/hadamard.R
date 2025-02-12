library("tidyverse")

# just some rudimentary helper functions -- there are various R packages that
# are more efficient and full-featured, but would be overkill for the handful of
# scenarios in which we need to construct our own Hadamard matrix in the context
# of a BRR analysis

is_whole <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  !is.character(all.equal.numeric(x %% 1, 0.0))
}

is_even <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  !is.character(all.equal.numeric(x %% 2, 0.0))
}

is_power <- function(x, base, tolerance = sqrt(.Machine$double.eps)) {
  remainder <- log(x, base = base) %% 1
  !is.character(all.equal.numeric(remainder, 0.0))
}

#' Hadamard matrix of order 2^k with Sylvester's method
#'
#' @param order order of the matrix
#'
#' @export
sylvester_matrix <- function(order) {
  if (order < 1 | order != 1 & !is_power(order, 2)) {
    stop("The order of a Sylvester Hadamard matrix must be a positive integer power of two.")
  }

  k <- log(order, base = 2)

  if (k == 0) {
    matrix(1)
  } else if (k == 1) {
    matrix(c(1, 1, 1, -1), nrow = 2, ncol = 2)
  } else if (k > 1) {
    kronecker(sylvester(2), sylvester(2^(k - 1)))
  }
}

# quadratic residues from the finite field of odd primes
quadratic_residues <- function(p) {
  map_int(seq(1, p - 1, by = 2), function(i) {
    i^2 %% p
  })
}

paley_prime_matrix <- function(p) {
  q <- quadratic_residues_gfp(p)
  upper_ixs <- seq(1, (p - 1) / 2)
  ixs <- list(rep(1:p, times = p), rep(1:p, each = p))

  matrix(pmap_int(ixs, function(i, j) {
    loc <- (i - j) %% p
    if (loc == 0) {
      0
    } else if (loc %in% q[upper_ixs]) {
      1
    } else {
      -1
    }
  }), nrow = p, ncol = p)
}

# primes where mod(prime) %% 4 == 3, listed here up to 200
paley_primes <- c(
  3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83,
  103, 107, 127, 131, 139, 151, 163, 167, 179, 191, 199
  )

#' Hadamard matrix with Paley's method
#'
#' @param order order of the matrix
#'
#' @description
#' Paley's method works for orders that are one larger than a prime number,
#' where that prime divided by 4 has a remainder of 3, so e.g. 80 because
#' 79 is a prime number and 79 divided by 4 has a remainder of 3.
#'
#' Adapted from code originally written Appavoo Dhandapani and Revan Siddesha.
#'
#' @export
paley_matrix <- function(order) {
  q <- order - 1

  if (!(q %in% paley_primes)) {
    options <- str_flatten_comma(paley_primes + 1)
    stop(str_glue("The order of a Paley Hadamard matrix must be one of: {options}."))
  }

  body <- paley_prime_matrix(q)
  top <- matrix(c(0, rep(1, q)), nrow = 1, ncol = order)
  left <- matrix(c(rep(-1,q)),nrow = q, ncol = 1)
  complete <- rbind(top , cbind(left, body))
  complete + diag(order)
}

#' Hadamard matrix
#'
#' @param order order of the matrix
#'
#' @description
#' Compute a Hadamard matrix with an order of any power of two, plus
#' a handful of other useful sizes including 20, 60, 80, 140, 180, 200.
#'
#' @export
hadamard_matrix <- function(order) {
  if (order < 1 | order != 1 & !is_even(order)) {
    stop("The order of a Hadamard matrix must be a positive integer.")
  }

  if (order == 1 | is_even(order) & is_whole(log(order, base = 2))) {
    sylvester_matrix(order)
  } else if ((order - 1) %in% paley_primes) {
    paley_matrix(order)
  } else {
    stop(str_glue("No algorithm available for a Hadamard matrix of order {order}"))
  }
}

#' Balanced repeated replication weight matrix with Fay's correction
#'
#' @param ncol order of the matrix
#' @param nrow length of the matrix
#' @param ncol width of the matrix
#' @param scale scale by which to multiply the original entries (-1 and 1)
#' @param location location by which to shift the scaled entries
#' @param weights sampling weights or a single weight to further adjust the weights of each row
#'
#' @description
#' This replication weight matrix generalizes the Hadamard matrix with arbitrary
#' weights (or weight adjustment factors) instead of -1 and 1. The default
#' scale and location parameters produce weights of 0.5 and 1.5, as with Fay's
#' method.
#'
#' Using the `weights` argument,these replications weights can then further be
#' multiplied with the sampling weights so they can be used in a balanced
#' repeated replication.
#'
#' If `nrow` is different from `ncol`, or if sampling weights are provided,
#' the matrix will be repeated or abbreviated to attain the desired amount of
#' rows.
#'
#' @export
fay_matrix <- function(ncol = 80, nrow = ncol, scale = 0.5, location = 1, weights = 1) {
  if (length(weights) > 1) nrow <- length(weights)

  block <- location + hadamard_matrix(ncol) * scale

  if (nrow < ncol) {
    weights * block[1:nrow, ]
  } else if (nrow > ncol) {
    n <- ceiling(nrow / ncol)
    ixs <- head(rep(1:ncol, times = n), n = nrow)
    weights * block[ixs, ]
  } else {
    weights * block
  }
}
