# Updated: 2021-08-30

#' Ridge regression.
#'
#' @param X Covariate matrix.
#' @param y Numeric outcome vector.
#' @param weights Numeric vector of weights.
#' @param lambda Penalization parameter.
#' @param family Exponential family of interest.
#' @importFrom glmnet glmnet
#' @export
#' @return Vector containing regression coefficients.
#'

RidgeRegression <- function(X, y, weights = NULL, lambda = 1e-04,
                            family = "binomial") {
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  # Multiple lambdas to make sure convergence happens.
  gamma <- tryCatch(
    (coef(glmnet::glmnet(X, y,
      weights = weights, alpha = 0,
      lambda = seq(lambda, 100 * lambda, length.out = 100),
      family = binomial()
    ))),
    error = function(e) rep(NA, ncol(X) + 1)
  )

  # Take the value at specified value of lambda.
  gamma <- gamma[, ncol(gamma)]

  return(gamma)
}
