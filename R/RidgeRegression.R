#' Ridge regression.
#'
#' @param X Covariate matrix.
#' @param y Numeric outcome vector.
#' @param coef Coefficient for scaling lambda values.
#' @param weights Numeric vector of weights.
#' @param family Exponential family of interest.
#' @param exponents Sequence of exponents for generating lambda values.
#' @importFrom glmnet glmnet cv.glmnet
#' @export
#' @return Vector containing regression coefficients.
#'

# RidgeRegression <- function(X,
#                             y,
#                             coef,
#                             weights = NULL,
#                             family = "binomial",
#                             exponents = seq(1.05, 2.5, length.out = 100)) {
#   if (is.null(weights)) {
#     weights <- rep(1, length(y))
#   }
# 
#   n <- length(y) # Number of samples
# 
#   # Generate original lambda values as o(n^{-1}), scaled by larger constants
#   original_lambdas <- sapply(exponents, function(exp) coef / n^exp)
# 
#   # Generate warm start values
#   warm_start_lambdas <- seq(100 * original_lambdas[1], original_lambdas[1],
#     length.out = 100
#   )[-100]
# 
#   # Combine warm start lambdas and original lambdas
#   lambdas <- c(warm_start_lambdas, original_lambdas)
# 
#   cv_fit <- tryCatch(
#     glmnet::cv.glmnet(X, y,
#       weights = weights, alpha = 0, lambda = lambdas,
#       family = family
#     ),
#     error = function(e) {
#       message("Error in fitting Ridge Regression: ", e$message)
#       return(NULL)
#     }
#   )
# 
# 
#   if (is.null(cv_fit)) {
#     return(rep(NA, ncol(X) + 1)) # Return NA vector on error
#   }
# 
#   # Extract cross-validated MSE or deviance for all lambdas
#   cv_mse <- cv_fit$cvm
# 
#   # Extract the lambda values used in cross-validation
#   lambda_values <- cv_fit$lambda
# 
#   # Find the index of the minimum cross-validation error (MSE or deviance) for
#   # lambdas in original lambdas only
#   valid_indices <- which(lambda_values %in% original_lambdas)
#   best_index <- valid_indices[which.min(cv_mse[valid_indices])]
# 
#   # Extract the best lambda from the original lambdas
#   best_lambda <- lambda_values[best_index]
# 
#   # Find the corresponding exponent for the best lambda from original lambdas
#   best_lambda_index <- which(original_lambdas == best_lambda)
#   best_exp <- exponents[best_lambda_index]
# 
#   # Print the exponent for the best lambda
#   message("Exponent corresponding to best lambda: ", best_exp)
# 
#   # Extract coefficients at the best lambda
#   gamma <- coef(cv_fit, s = best_lambda)
# 
#   return(as.numeric(gamma))
# }


RidgeRegression <- function(X, y, coef = NULL, weights = NULL, 
                          exponents = seq(1 + 0.05, 2, length.out = 100)) {
  if (is.null(weights)) weights <- rep(1, length(y))
  n <- length(y)
  p <- ncol(X)
  
  # Data-adaptive scaling for lambda
  c <- max(abs(t(X) %*% (y - mean(y)))) / n
  
  # Lambda grid: o(n^{-1}) (exponents > 1)
  original_lambdas <- c / (n^exponents)
  
  # Warm start grid (larger lambdas for convergence)
  warm_start_lambdas <- exp(seq(log(100 * original_lambdas[1]), 
                                log(original_lambdas[1]), 
                                length.out = 50))
  
  # Fit model (use warm start + original lambdas)
  cv_fit <- glmnet::cv.glmnet(X, y, weights = weights, alpha = 0, 
                              lambda = c(warm_start_lambdas, original_lambdas), 
                              family = "binomial")
  
  # --- Key Fix: Restrict best lambda to original_lambdas ---
  # Get all CV results
  lambda_values <- cv_fit$lambda
  cv_mse <- cv_fit$cvm
  
  # Find index of best lambda *only in original_lambdas*
  valid_indices <- which(lambda_values %in% original_lambdas)
  best_index <- valid_indices[which.min(cv_mse[valid_indices])]
  best_lambda <- lambda_values[best_index]
  
  # Get corresponding exponent
  best_exponent <- exponents[which.min(abs(original_lambdas - best_lambda))]
  #print(sprintf("Selected exponent: %.3f (lambda = %.2e)", best_exponent, best_lambda))
  #print(sprintf("n * lambda = %.3f (should -> 0 as n -> ∞)", n * best_lambda))
  
  return(as.numeric(coef(cv_fit, s = best_lambda)))
}