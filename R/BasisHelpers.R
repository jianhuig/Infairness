# Purpose: Internal basis-construction helpers.

polynomial <- function(data, order) {
  polynomial_Z <- rep(1, nrow(data))
  for (i in 1:order) {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z, polynomial_i)
  }
  return(polynomial_Z)
}

gbic_ppo <- function(Y_new, covariates_matrix, additional_matrix, order_poly) {
  d <- ncol(Y_new)
  n <- nrow(Y_new)
  p <- ncol(covariates_matrix)
  q <- ncol(additional_matrix)

  polynomial_matrix <- polynomial(covariates_matrix, order_poly)
  polynomial_matrix <- cbind(polynomial_matrix, additional_matrix)

  gammahat <- lm(Y_new ~ polynomial_matrix - 1)$coefficients
  residual_square <- (Y_new - polynomial_matrix %*% gammahat)^2
  sigmahat_square <- colSums(residual_square) / (n - p * order_poly - q - 1)
  design_matrix <- t(polynomial_matrix) %*% polynomial_matrix / n

  if (class(try(solve(design_matrix), silent = TRUE))[1] == "try-error") {
    return(9999999)
  }

  trace_det_AB <- numeric()
  for (j in 1:d) {
    Bhat_j <- t(polynomial_matrix * residual_square[, j]) %*%
      polynomial_matrix / n
    cova_constrast <- 1 / sigmahat_square[j] * solve(design_matrix) %*% Bhat_j
    trace_AB <- sum(diag(cova_constrast))
    det_AB <- det(cova_constrast)
    trace_det_AB <- c(trace_det_AB, trace_AB - log(abs(det_AB)))
  }

  1 / n * (
    d * (n - p * order_poly - q - 1) +
      n * sum(log(sigmahat_square)) +
      d * log(n) * (p * order_poly + q + 1) +
      sum(trace_det_AB)
  )
}
