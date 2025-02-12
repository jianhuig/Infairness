#' Polynomial basis with order minimizing generalized BIC
#'
#' @param Y binary response variable
#' @param covariates_matrix covariates matrix
#' @param additional_matrix additional matrix
#' @param gamma highest polynomial order, default is 10
#' @param weights perturbation weights
#' @return optimal polynomial order
#'
#' @export

find_alpha_glm <- function(Y, covariates_matrix, additional_matrix, gamma = 10,
                           weights = NULL) {
  n <- length(Y)

  # cross fitting
  # supervised theta
  theta_part1 <- glm(
    Y[1:round(n / 2)] ~ cbind(
      covariates_matrix[1:round(n / 2), ],
      additional_matrix[1:round(n / 2), ]
    ),
    family = binomial(link = "logit"),
    weights = weights[1:round(n / 2)]
  )
  theta_part2 <- glm(
    Y[(round(n / 2) + 1):n] ~ cbind(
      covariates_matrix[(round(n / 2) + 1):n, ],
      additional_matrix[(round(n / 2) + 1):n, ]
    ),
    family = binomial(link = "logit"),
    weights = weights[(round(n / 2) + 1):n]
  )

  exp_numerator_part1 <- as.vector(exp(cbind(
    1, covariates_matrix[1:round(n / 2), ],
    additional_matrix[1:round(n / 2), ]
  ) %*% theta_part1$coefficients))

  exp_numerator_part2 <- as.vector(exp(cbind(
    1, covariates_matrix[(round(n / 2) + 1):n, ],
    additional_matrix[(round(n / 2) + 1):n, ]
  ) %*% theta_part2$coefficients))

  L_first_derivative_part1 <- (exp_numerator_part1 / (1 + exp_numerator_part1)
    - Y[1:round(n / 2)]) * cbind(
    1, covariates_matrix[1:round(n / 2), ],
    additional_matrix[1:round(n / 2), ]
  )

  L_first_derivative_part2 <- (exp_numerator_part2 / (1 + exp_numerator_part2)
    - Y[(round(n / 2) + 1):n]) * cbind(
    1, covariates_matrix[(round(n / 2) + 1):n, ],
    additional_matrix[(round(n / 2) + 1):n, ]
  )

  L_first_derivative <- rbind(L_first_derivative_part1, L_first_derivative_part2)

  alpha <- which.min(sapply(
    1:gamma,
    function(t) {
      GBIC_ppo(
        L_first_derivative,
        covariates_matrix,
        additional_matrix,
        t
      )
    }
  ))

  return(alpha)
}


# Polynomial basis
polynomial <- function(data, order) {
  #---------------------------Arguments----------------------------------------#
  # Purpose: This function is to produce the polynomial basis of X including
  #           the intercept vector.
  #
  # Input:
  #       data: A matrix, whose each row is an observation of predictor vector.
  #       order: The polynomial order.
  #----------------------------------------------------------------------------#
  polynomial_Z <- rep(1, nrow(data))
  for (i in 1:order)
  {
    polynomial_i <- data^i
    polynomial_Z <- cbind(polynomial_Z, polynomial_i)
  }
  return(polynomial_Z)
}


############### (1.2) the data driven selector GBIC_ppo
GBIC_ppo <- function(Y_new, covariates_matrix, additional_matrix, order_poly) {
  #-----------------------Arguments--------------------------------------------------#
  # Purpose: GBIC_ppo criterion function
  #
  # Input:
  #       Y_new: The first derivatives of the loss function L.
  #       covariates_matrix: A matrix, each row is an observation of predictor vector.
  #       order_poly: The polynomial order.
  #
  #---------------------------------------------------------------------------------#
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


  if (class(try(solve(design_matrix), silent = T))[1] == "try-error") {
    GBIC <- 9999999
  } else {
    trace_det_AB <- numeric()
    for (j in 1:d) {
      Bhat_j <- t(polynomial_matrix * residual_square[, j]) %*%
        polynomial_matrix / n
      cova_constrast <- 1 / sigmahat_square[j] * solve(design_matrix) %*% Bhat_j
      trace_AB <- sum(diag(cova_constrast))
      det_AB <- det(cova_constrast)
      trace_det_AB <- c(trace_det_AB, trace_AB - log(abs(det_AB)))
    }
    GBIC <- 1 / n * (d * (n - p * order_poly - q - 1) +
      n * sum(log(sigmahat_square)) +
      d * log(n) * (p * order_poly + q) +
      sum(trace_det_AB))
  }

  return(GBIC)
}
