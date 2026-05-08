# Purpose: Backward-compatible spline helpers.

ns.basis <- function(X, nk) {
  X <- as.matrix(X)
  basis <- NULL

  for (i in 1:ncol(X)) {
    X_i <- X[, i]
    knots <- quantile(X_i, seq(0, 1, length = nk))

    j <- 0
    while (length(unique(knots)) != nk) {
      j <- j + 1
      knots <- unique(quantile(X_i, seq(0, 1, length = nk + j)))
    }

    d_k <- (trunc.cub(X_i, knots[nk - 1])) / (knots[nk] - knots[nk - 1])
    evals <- sapply(1:(nk - 2), function(ii) {
      d_i <- (trunc.cub(X_i, knots[ii])) / (knots[nk] - knots[ii])
      d_i - d_k
    })

    basis <- cbind(basis, cbind(X_i, evals))
  }

  return(basis)
}

#' Truncated cubic basis helper.
#'
#' @param X Variable of interest.
#' @param x Knot location.
#'
#' @export trunc.cub
trunc.cub <- function(X, x) {
  ((X > x) * (X - x))^3
}
