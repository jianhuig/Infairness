#' Ridge regression.
#'
#' @param X Covariate matrix.
#' @param y Numeric outcome vector.
#' @param weights Numeric vector of weights.
#' @param exponents Sequence of exponents for generating lambda values.
#' @param penalty_factor Optional per-column penalty factors passed to
#' `glmnet::cv.glmnet()`. Use 0 for unpenalized columns and 1 for the default
#' ridge penalty.
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


ridge_is_binary_column <- function(x,
                                   tol = sqrt(.Machine$double.eps)) {
  values <- unique(stats::na.omit(as.numeric(x)))

  length(values) > 0L &&
    length(values) <= 2L &&
    all(abs(values - round(values)) < tol)
}

ridge_binary_main_effect_penalty <- function(X,
                                             unpenalize_binary = FALSE) {
  if (is.null(X)) {
    return(numeric(0))
  }

  X <- as.matrix(X)
  if (!unpenalize_binary || ncol(X) == 0L) {
    return(rep(1, ncol(X)))
  }

  as.numeric(!vapply(
    seq_len(ncol(X)),
    function(j) ridge_is_binary_column(X[, j]),
    logical(1)
  ))
}

ridge_design_penalty_factor <- function(n_basis = 0L,
                                        D = NULL,
                                        X = NULL,
                                        X_interaction = NULL,
                                        unpenalize_binary_X = FALSE) {
  parts <- list(rep(1, n_basis))

  if (!is.null(D)) {
    parts[[length(parts) + 1L]] <- rep(1, ncol(as.matrix(D)))
  }
  if (!is.null(X)) {
    parts[[length(parts) + 1L]] <- ridge_binary_main_effect_penalty(
      X,
      unpenalize_binary = unpenalize_binary_X
    )
  }
  if (!is.null(X_interaction)) {
    parts[[length(parts) + 1L]] <- rep(1, ncol(as.matrix(X_interaction)))
  }

  unname(unlist(parts, use.names = FALSE))
}

RidgeRegression <- function(X, y, weights = NULL,
                            exponents = seq(1.05, 2.5, length.out = 100),
                            penalty_factor = NULL) {
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1, length(y))
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, ncol(X))
  }
  if (length(penalty_factor) != ncol(X)) {
    stop("`penalty_factor` must have one value per column of `X`.")
  }

  n <- length(y)
  
  # Data-adaptive scaling
  c0 <- max(abs(crossprod(X, y - mean(y)))) / n
  
  # Vanishing penalty grid
  original_lambdas <- c0 / (n^exponents)
  
  k <- 10
  foldid <- integer(n)
  foldid[y == 1] <- sample(rep(1:k, length.out = sum(y == 1)))
  foldid[y == 0] <- sample(rep(1:k, length.out = sum(y == 0)))
  
  cv_fit <- glmnet::cv.glmnet(
    X, y,
    weights = weights,
    alpha = 0,
    foldid = foldid,
    standardize = TRUE,
    lambda = sort(original_lambdas, decreasing = TRUE),
    penalty.factor = penalty_factor,
    family = "binomial"
  )
  
  lambda_values <- cv_fit$lambda
  cv_mse <- cv_fit$cvm
  
  # safer matching
  valid_indices <- sapply(original_lambdas, function(lam) {
    which.min(abs(lambda_values - lam))
  })
  valid_indices <- unique(valid_indices)
  
  best_index <- valid_indices[which.min(cv_mse[valid_indices])]
  best_lambda <- lambda_values[best_index]
  best_exponent <- exponents[which.min(abs(original_lambdas - best_lambda))]
  
  list(
    coefficients = as.numeric(stats::coef(cv_fit, s = best_lambda)),
    best_lambda = best_lambda,
    best_exponent = best_exponent
  )
}
