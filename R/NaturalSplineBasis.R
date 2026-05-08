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

# Use a slightly more conservative ridge path for spline-expanded designs.
SplineRidgeRegression <- function(X, y, weights = NULL, penalty_factor = NULL) {
  RidgeRegression(
    X = X,
    y = y,
    weights = weights,
    exponents = seq(0.8, 2.8, length.out = 120),
    penalty_factor = penalty_factor
  )
}

spline_indicator_term <- function(C, ...) {
  as.matrix(C)
}

build_spline_predictors <- function(basis_matrix,
                                    D = NULL,
                                    X = NULL,
                                    X_interaction = NULL) {
  parts <- list(as.matrix(basis_matrix))

  if (!is.null(D)) {
    parts[[length(parts) + 1]] <- as.matrix(D)
  }
  if (!is.null(X)) {
    parts[[length(parts) + 1]] <- as.matrix(X)
  }
  if (!is.null(X_interaction)) {
    parts[[length(parts) + 1]] <- as.matrix(X_interaction)
  }

  as.matrix(do.call(cbind, parts))
}

build_spline_penalty_factor <- function(basis_matrix,
                                        D = NULL,
                                        X = NULL,
                                        X_interaction = NULL,
                                        unpenalize_binary_X = FALSE) {
  ridge_design_penalty_factor(
    n_basis = ncol(as.matrix(basis_matrix)),
    D = D,
    X = X,
    X_interaction = X_interaction,
    unpenalize_binary_X = unpenalize_binary_X
  )
}

fit_spline_model <- function(X,
                             y,
                             weights = NULL,
                             penalty_factor = NULL,
                             ...) {
  X <- as.matrix(X)
  fit <- SplineRidgeRegression(
    X = X,
    y = y,
    weights = weights,
    penalty_factor = penalty_factor
  )
  list(type = "linear", coefficients = fit$coefficients)
}

predict_spline_model <- function(model,
                                 X,
                                 ...) {
  boot::inv.logit(as.matrix(cbind(1, X)) %*% model$coefficients)
}
