# Purpose: Public basis-selection utilities.

#' Polynomial order selection via generalized BIC.
#'
#' @param Y Binary response variable.
#' @param covariates_matrix Covariates matrix.
#' @param additional_matrix Additional matrix.
#' @param gamma Candidate polynomial orders.
#' @param weights Optional perturbation weights.
#'
#' @return Optimal polynomial order.
#' @export
find_alpha_glm <- function(Y, covariates_matrix, additional_matrix = NULL,
                           gamma = 1:10, weights = NULL) {
  n <- length(Y)

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
      additional_matrix[(round(n / 2) + 1):n]
    ),
    family = binomial(link = "logit"),
    weights = weights[(round(n / 2) + 1):n]
  )

  exp_numerator_part1 <- as.vector(exp(cbind(
    1,
    covariates_matrix[1:round(n / 2), ],
    additional_matrix[1:round(n / 2), ]
  ) %*% theta_part2$coefficients))

  exp_numerator_part2 <- as.vector(exp(cbind(
    1,
    covariates_matrix[(round(n / 2) + 1):n, ],
    additional_matrix[(round(n / 2) + 1):n, ]
  ) %*% theta_part1$coefficients))

  L_first_derivative_part1 <- (exp_numerator_part1 / (1 + exp_numerator_part1) -
    Y[1:round(n / 2)]) * cbind(
    1,
    covariates_matrix[1:round(n / 2), ],
    additional_matrix[1:round(n / 2), ]
  )

  L_first_derivative_part2 <- (exp_numerator_part2 / (1 + exp_numerator_part2) -
    Y[(round(n / 2) + 1):n]) * cbind(
    1,
    covariates_matrix[(round(n / 2) + 1):n, ],
    additional_matrix[(round(n / 2) + 1):n, ]
  )

  L_first_derivative <- rbind(L_first_derivative_part1, L_first_derivative_part2)

  gbics <- sapply(
    gamma,
    function(t) {
      gbic_ppo(
        L_first_derivative,
        covariates_matrix,
        additional_matrix,
        t
      )
    }
  )

  gbics[is.na(gbics)] <- Inf

  if (all(gbics >= 9e6)) {
    alpha <- 1
  } else {
    alpha <- which.min(gbics)
  }

  return(alpha)
}

#' Compute generalized BIC for a candidate design matrix.
#'
#' @param Y Binary response variable.
#' @param covariates_matrix Candidate design matrix.
#' @param weights Optional perturbation weights.
#' @param P_total Total number of available predictors.
#' @param gamma_ebic EBIC penalty parameter.
#' @param penalize_intercept Whether to include the intercept in the penalty.
#'
#' @return Generalized BIC score.
#' @export
compute_gbic <- function(Y, covariates_matrix, weights = NULL,
                         P_total, gamma_ebic = 0.5,
                         penalize_intercept = FALSE) {
  n0 <- length(Y)
  if (is.null(weights)) {
    weights <- rep(1, n0)
  }

  idx <- sample.int(n0)
  i1 <- idx[1:floor(n0 / 2)]
  i2 <- idx[(floor(n0 / 2) + 1):n0]

  fit1 <- tryCatch(
    glm(
      Y[i1] ~ covariates_matrix[i1, -1, drop = FALSE],
      family = binomial(),
      weights = weights[i1]
    ),
    error = function(e) NULL
  )
  fit2 <- tryCatch(
    glm(
      Y[i2] ~ covariates_matrix[i2, -1, drop = FALSE],
      family = binomial(),
      weights = weights[i2]
    ),
    error = function(e) NULL
  )
  if (is.null(fit1) || is.null(fit2)) {
    return(Inf)
  }

  mu1 <- as.numeric(plogis(covariates_matrix[i1, ] %*% fit2$coefficients))
  mu2 <- as.numeric(plogis(covariates_matrix[i2, ] %*% fit1$coefficients))

  Y_new <- rbind(
    (mu1 - Y[i1]) * covariates_matrix[i1, ],
    (mu2 - Y[i2]) * covariates_matrix[i2, ]
  )

  d <- ncol(Y_new)
  n <- nrow(Y_new)
  X <- covariates_matrix
  k <- ncol(X)
  df <- n - k
  if (df <= 5) {
    return(Inf)
  }

  gammahat <- tryCatch(lm(Y_new ~ X - 1)$coefficients, error = function(e) NULL)
  if (is.null(gammahat)) {
    return(Inf)
  }

  resid2 <- (Y_new - X %*% gammahat)^2
  sigma2 <- colSums(resid2) / df
  if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    return(Inf)
  }

  A <- crossprod(X) / n
  Ainv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(Ainv)) {
    return(Inf)
  }

  trace_det_sum <- 0
  for (j in 1:d) {
    B <- crossprod(X * resid2[, j], X) / n
    M <- (Ainv %*% B) / sigma2[j]

    tr <- sum(diag(M))
    ld <- determinant(M, logarithm = TRUE)
    if (!is.finite(tr) || !is.finite(ld$modulus)) {
      return(Inf)
    }

    trace_det_sum <- trace_det_sum + (tr - as.numeric(ld$modulus))
  }

  if (penalize_intercept) {
    k_eff <- k
    P_eff <- P_total
  } else {
    k_eff <- k - 1
    P_eff <- P_total - 1
  }

  if (k_eff < 0 || P_eff <= 0 || k_eff > P_eff) {
    return(Inf)
  }

  penalty <- d * (log(n) * k_eff + 2 * gamma_ebic * lchoose(P_eff, k_eff))

  GBIC <- (d * df + n * sum(log(sigma2)) + penalty + trace_det_sum) / n
  if (!is.finite(GBIC)) {
    GBIC <- Inf
  }
  GBIC
}
