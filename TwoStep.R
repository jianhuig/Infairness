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
GBIC_ppo <- function(Y_new, covariates_matrix, order_poly) {
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

  polynomial_matrix <- polynomial(covariates_matrix, order_poly)

  gammahat <- lm(Y_new ~ polynomial_matrix - 1)$coefficients
  residual_square <- (Y_new - polynomial_matrix %*% gammahat)^2
  sigmahat_square <- colSums(residual_square) / (n - p * order_poly - 1)
  design_matrix <- t(polynomial_matrix) %*% polynomial_matrix / n


  if (class(try(solve(design_matrix), silent = T))[1] == "try-error") {
    GBIC <- 9999999
  } else {
    trace_det_AB <- numeric()
    for (j in 1:d) {
      Bhat_j <- t(polynomial_matrix * residual_square[, j]) %*% polynomial_matrix / n
      cova_constrast <- 1 / sigmahat_square[j] * solve(design_matrix) %*% Bhat_j
      trace_AB <- sum(diag(cova_constrast))
      det_AB <- det(cova_constrast)
      trace_det_AB <- c(trace_det_AB, trace_AB - log(det_AB))
    }
    GBIC <- 1 / n * (d * (n - p * order_poly - 1) + n * sum(log(sigmahat_square)) + d * log(n) * p * order_poly + sum(trace_det_AB))
  }

  return(GBIC)
}

find_alpha <- function(Y, covariates_matrix, gamma = 10) {
  n <- length(Y)

  # cross fitting
  # supervised theta
  theta_part1 <- glm(
    Y[1:round(n / 2)] ~ covariates_matrix[1:round(n / 2), ],
    family = binomial(link = "logit")
  )
  theta_part2 <- glm(
    Y[(round(n / 2) + 1):n] ~ covariates_matrix[(round(n / 2) + 1):n, ],
    family = binomial(link = "logit")
  )

  exp_numerator_part1 <- as.vector(exp(cbind(
    1, covariates_matrix[1:round(n / 2), ]
  ) %*% theta_part1$coefficients))

  exp_numerator_part2 <- as.vector(exp(cbind(
    1, covariates_matrix[(round(n / 2) + 1):n, ]
  ) %*% theta_part2$coefficients))

  L_first_derivative_part1 <- (exp_numerator_part1 / (1 + exp_numerator_part1)
    - Y[1:round(n / 2)]) * cbind(
    1, covariates_matrix[1:round(n / 2), ]
  )

  L_first_derivative_part2 <- (exp_numerator_part2 / (1 + exp_numerator_part2)
    - Y[(round(n / 2) + 1):n]) * cbind(
    1, covariates_matrix[(round(n / 2) + 1):n, ]
  )

  L_first_derivative <- rbind(L_first_derivative_part1, L_first_derivative_part2)

  alpha <- which.min(sapply(
    1:gamma,
    function(t) {
      print(GBIC_ppo(
        L_first_derivative,
        covariates_matrix,
        t
      ))
    }
  ))

  return(alpha)
}

step_I <- function(Y_lab,
                   A_lab,
                   A_unlab,
                   basis,
                   S_lab = NULL,
                   S_unlab = NULL,
                   X_lab = NULL,
                   X_unlab = NULL) {
  nclass <- sort(unique(A_lab))

  if (basis == "S") {
    m_unlabeled <- S_unlab # imputed values
    m_labeled <- S_lab # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha(
        Y = Y_lab[A_lab == a],
        covariates_matrix = S_lab[A_lab == a] %>% as.matrix()
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_lab %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlab %>% as.data.frame(), alpha)

      gamma <- glm(
        Y_lab[A_lab == a] ~ basis_labeled[A_lab == a, -1] %>%
          as.matrix(),
        family = binomial()
      )$coef

      imputed_unlabeled <- boot::inv.logit(basis_unlabeled[A_unlab == a, ] %>%
        as.matrix() %*% gamma)

      m_unlabeled[A_unlab == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(basis_labeled[A_lab == a, ] %>%
        as.matrix() %*% gamma)

      m_labeled[A_lab == a] <- imputed_labeled
    }
  } else if (basis == "S,X") {
    m_unlabeled <- S_unlab # imputed values
    m_labeled <- S_lab # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha(
        Y = Y_lab[A_lab == a],
        covariates_matrix = cbind(S_lab[A_lab == a], X_lab[A_lab == a, ])
        %>% as.matrix()
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(cbind(S_lab, X_lab) %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(cbind(S_unlab, X_unlab) %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = basis_labeled[A_lab == a, -1] %>% as.matrix(),
            y = Y_lab[A_lab == a],
            coef = log(ncol(basis_labeled) + 1)
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(basis_unlabeled[A_unlab == a, ] %>%
        as.matrix() %*% gamma)

      m_unlabeled[A_unlab == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(basis_labeled[A_lab == a, ] %>%
        as.matrix() %*% gamma)

      m_labeled[A_lab == a] <- imputed_labeled
    }
  } else if (basis == "X") {
    m_unlabeled <- S_unlab # imputed values
    m_labeled <- S_lab # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha(
        Y = Y_lab[A_lab == a],
        covariates_matrix = X_lab[A_lab == a, ] %>% as.matrix()
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(X_lab %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(X_unlab %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = basis_labeled[A_lab == a, -1] %>% as.matrix(),
            y = Y_lab[A_lab == a],
            coef = log(ncol(basis_labeled) + 1)
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(basis_unlabeled[A_unlab == a, ] %>%
        as.matrix() %*% gamma)

      m_unlabeled[A_unlab == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(basis_labeled[A_lab == a, ] %>%
        as.matrix() %*% gamma)

      m_labeled[A_lab == a] <- imputed_labeled
    }
  }

  return(list(m_labeled = m_labeled, m_unlabeled = m_unlabeled, alphas = alphas))
}


step_II <- function(Y_lab, A_lab, A_unlab, basis, m_lab, m_unlab, D_lab,
                    D_unlab, S_lab = NULL, S_unlab = NULL) {
  nclass <- sort(unique(A_lab))

  if (basis == "S" | basis == "S,X") {
    imp_unlab <- m_unlab # imputed values
    imp_lab <- m_lab # imputed values

    for (a in nclass) {
      gamma <- glm(Y_lab[A_lab == a] ~ D_lab[A_lab == a] - 1,
        offset = m_lab[A_lab == a],
        family = binomial()
      )$coef
      
      print(gamma)

      imp_unlab[A_unlab == a] <-
        boot::inv.logit(
          cbind(D_unlab[A_unlab == a]) %*% gamma +
            m_unlab[A_unlab == a]
        )

      imp_lab[A_lab == a] <-
        boot::inv.logit(
          cbind(D_lab[A_lab == a]) %*% gamma +
            m_lab[A_lab == a]
        )
    }
  } else if (basis == "X") {
    imp_unlab <- m_unlab # imputed values
    imp_lab <- m_lab # imputed values

    for (a in nclass) {
      gamma <- glm(Y_lab[A_lab == a] ~ D_lab[A_lab == a] + S_lab[A_lab == a],
        offset = m_lab[A_lab == a],
        family = binomial()
      )$coef
      
      print(summary(glm(Y_lab[A_lab == a] ~ D_lab[A_lab == a] + S_lab[A_lab == a] ,
                        offset = m_lab[A_lab == a],
                        family = binomial()
      )))

      imp_unlab[A_unlab == a] <- boot::inv.logit(cbind(
        1,
        D_unlab[A_unlab == a],
        S_unlab[A_unlab == a]
      ) %*% gamma +
        m_unlab[A_unlab == a])

      imp_lab[A_lab == a] <- boot::inv.logit(cbind(
        1,
        D_lab[A_lab == a],
        S_lab[A_lab == a]
      ) %*% gamma + m_lab[A_lab == a])
    }
  }
  return(list(m_labeled = imp_lab, m_unlabeled = imp_unlab))
}

get_all_metric <- function(A, S, D, imp1, imp2, basis) {
  class <- sort(unique(A))
  out <- c()
  for (a in class) {
    mu_Y <- mean(imp1[A == a])
    mu_D <- mean(D[A == a])
    mu_S2 <- mean(S[A == a]^2)
    if (basis %in% c("S", "S,X")) {
      mu_SY <- mean(S[A == a] * imp1[A == a])
    } else if (basis == "X") {
      mu_SY <- mean(S[A == a] * imp2[A == a])
    }
    mu_DY <- mean(D[A == a] * imp2[A == a])

    tpr <- mu_DY / mu_Y
    fnr <- 1 - tpr
    fpr <- (mu_D - mu_DY) / (1 - mu_Y)
    tnr <- 1 - fpr
    ppv <- mu_DY / mu_D
    npv <- (1 - mu_D - mu_Y + mu_DY) / (1 - mu_D)
    f1 <- 2 * mu_DY / (mu_D + mu_Y)
    acc <- 1 - mu_Y - mu_D + 2 * mu_DY
    bs <- mu_S2 - 2 * mu_SY + mu_Y
    out <- c(out, tpr, tnr, fpr, fnr, npv, ppv, acc, f1, bs)
  }
  out <- cbind(matrix(out, ncol = 2, byrow = FALSE), NA)
  out[, 3] <- out[, 1] - out[, 2]
  colnames(out) <- c(paste0("Group", class), "Delta")
  rownames(out) <- c("TPR", "TNR", "FPR", "FNR", "NPV", "PPV", "ACC", "F1", "BS")
  tibble::rownames_to_column(as.data.frame(out), "Metric")
}

Infairness_twostep <- function(Y,
                               S,
                               A,
                               threshold = 0.5,
                               X = NULL,
                               basis = "S",
                               W = NULL) {
  if (!is.null(W)) {
    W_label <- 4 * rbeta(sum(!is.na(Y)), 1 / 2, 3 / 2)
  } else {
    W_label <- NULL
  }

  labeled_ind <- which(!is.na(Y))
  Y_lab <- Y[labeled_ind]
  A_lab <- A[labeled_ind]
  A_unlab <- A[-labeled_ind]

  S_lab <- S[labeled_ind]
  S_unlab <- S[-labeled_ind]

  if (!is.null(X)) {
    X_lab <- X[labeled_ind, , drop = FALSE]
    X_unlab <- X[-labeled_ind, , drop = FALSE]
  } else {
    X_lab <- NULL
    X_unlab <- NULL
  }

  D <- ifelse(S > threshold, 1, 0)
  D_lab <- D[labeled_ind]
  D_unlab <- D[-labeled_ind]

  m_I <- step_I(
    Y_lab,
    A_lab,
    A_unlab,
    basis,
    S_lab,
    S_unlab,
    X_lab,
    X_unlab
  )

  m_II <- step_II(Y_lab, A_lab, A_unlab, basis,
    m_lab = m_I$m_labeled,
    m_unlab = m_I$m_unlabeled, D_lab, D_unlab, S_lab, S_unlab
  )

  return(get_all_metric(
    A = A_unlab,
    S = S_unlab,
    D = D_unlab,
    imp1 = m_I$m_unlabeled,
    imp2 = m_II$m_unlabeled,
    basis = basis
  ))
}
