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

find_alpha_glm <- function(Y, covariates_matrix, additional_matrix = NULL,
                           gamma = 1:10, weights = NULL) {
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
  ) %*% theta_part2$coefficients))

  exp_numerator_part2 <- as.vector(exp(cbind(
    1, covariates_matrix[(round(n / 2) + 1):n, ],
    additional_matrix[(round(n / 2) + 1):n, ]
  ) %*% theta_part1$coefficients))

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

  gbics <- sapply(
    gamma,
    function(t) {
      GBIC_ppo(
        L_first_derivative,
        covariates_matrix,
        additional_matrix,
        t
      )
    }
  )
  
  gbics[is.na(gbics)] <- Inf
  
  if (all(gbics >= 9e6)){
    alpha <- 1
  } else {
    alpha <- which.min(gbics)
  }
  
  return(alpha)
}


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

compute_gbic <- function(Y, covariates_matrix, weights = NULL,
                         P_total, gamma_ebic = 0.5,
                         penalize_intercept = FALSE) {
  n0 <- length(Y)
  if (is.null(weights)) weights <- rep(1, n0)
  
  # (Optional but recommended) random split
  idx <- sample.int(n0)
  i1 <- idx[1:floor(n0/2)]
  i2 <- idx[(floor(n0/2)+1):n0]
  
  # Fit logistic on each half (assumes covariates_matrix has intercept in col 1)
  fit1 <- tryCatch(
    glm(Y[i1] ~ covariates_matrix[i1, -1, drop=FALSE],
        family = binomial(), weights = weights[i1]),
    error = function(e) NULL
  )
  fit2 <- tryCatch(
    glm(Y[i2] ~ covariates_matrix[i2, -1, drop=FALSE],
        family = binomial(), weights = weights[i2]),
    error = function(e) NULL
  )
  if (is.null(fit1) || is.null(fit2)) return(Inf)
  
  # Cross-fitted means (use plogis for stability)
  mu1 <- as.numeric(plogis(covariates_matrix[i1, ] %*% fit2$coefficients))
  mu2 <- as.numeric(plogis(covariates_matrix[i2, ] %*% fit1$coefficients))
  
  # Score matrix Y_new = (mu - y) * X
  Y_new <- rbind((mu1 - Y[i1]) * covariates_matrix[i1, ],
                 (mu2 - Y[i2]) * covariates_matrix[i2, ])
  
  d <- ncol(Y_new)
  n <- nrow(Y_new)
  X <- covariates_matrix
  k <- ncol(X)
  df <- n - k
  if (df <= 5) return(Inf)
  
  gammahat <- tryCatch(lm(Y_new ~ X - 1)$coefficients, error = function(e) NULL)
  if (is.null(gammahat)) return(Inf)
  
  resid2 <- (Y_new - X %*% gammahat)^2
  sigma2 <- colSums(resid2) / df
  if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) return(Inf)
  
  A <- crossprod(X) / n
  Ainv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(Ainv)) return(Inf)
  
  trace_det_sum <- 0
  for (j in 1:d) {
    B <- crossprod(X * resid2[, j], X) / n
    M <- (Ainv %*% B) / sigma2[j]
    
    tr <- sum(diag(M))
    ld <- determinant(M, logarithm = TRUE)
    if (!is.finite(tr) || !is.finite(ld$modulus)) return(Inf)
    
    trace_det_sum <- trace_det_sum + (tr - as.numeric(ld$modulus))
  }
  
  # ---- EBIC penalty ----
  if (penalize_intercept) {
    k_eff <- k
    P_eff <- P_total
  } else {
    # common convention: don't penalize intercept
    k_eff <- k - 1
    P_eff <- P_total - 1
  }
  
  if (k_eff < 0 || P_eff <= 0 || k_eff > P_eff) return(Inf)
  
  penalty <- d * ( log(n) * k_eff + 2 * gamma_ebic * lchoose(P_eff, k_eff) )
  
  GBIC <- ( d*df + n*sum(log(sigma2)) + penalty + trace_det_sum ) / n
  if (!is.finite(GBIC)) GBIC <- Inf
  GBIC
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
      d * log(n) * (p * order_poly + q + 1) +
      sum(trace_det_AB))
  }

  return(GBIC)
}
