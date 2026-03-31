# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param num_knots Number of knots.
#' @param knots Optional list of knot vectors, one per column of `X`.
#' @param return_knots Logical; if `TRUE`, attach the knot locations used to the
#' returned matrix as an attribute.
#' @export
#' @return Matrix containing natural spline basis.
#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
#'
NaturalSplineBasis <- function(X, num_knots, knots = NULL, return_knots = FALSE) {
  X <- as.matrix(X)
  basis_parts <- vector("list", ncol(X))
  knots_used <- vector("list", ncol(X))

  for (i in 1:ncol(X)) {
    X_i <- X[, i]

    if (is.null(knots)) {
      current_knots <- quantile(X_i, seq(0, 1, length = num_knots))
      j <- 0

      while (length(unique(current_knots)) != num_knots) {
        j <- j + 1
        current_knots <- unique(quantile(X_i, seq(0, 1, length = num_knots + j)))
      }
    } else {
      current_knots <- knots[[i]]
    }
    knots_used[[i]] <- current_knots

    # Compute the natural spline basis.
    d_k <- (TruncatedCubic(X_i, current_knots[num_knots - 1])) /
      (current_knots[num_knots] - current_knots[num_knots - 1])
    evals <- sapply(1:(num_knots - 2), function(ii) {
      d_i <- (TruncatedCubic(X_i, current_knots[ii])) / (current_knots[num_knots] - current_knots[ii])
      basis.new <- d_i - d_k
    })
    evals <- matrix(evals, nrow = length(X_i))

    # Bind original variable and basis.
    basis_parts[[i]] <- cbind(X_i, evals)
  }

  # Return basis including everything.
  basis <- do.call(cbind, basis_parts)
  if (return_knots) {
    attr(basis, "knots") <- knots_used
  }

  return(basis)
}


TruncatedCubic <- function(x, knot_location) {
  return(((x > knot_location) * (x - knot_location))^3)
}

# Use a slightly more conservative ridge path for spline-expanded designs,
# which typically have wider bases than the polynomial branches.
SplineRidgeRegression <- function(X, y, weights = NULL) {
  RidgeRegression(
    X = X,
    y = y,
    weights = weights,
    exponents = seq(0.8, 2.8, length.out = 120)
  )
}
