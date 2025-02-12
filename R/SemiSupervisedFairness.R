# Purpose: Estimation of fairness metrics in semi-supervised setting.

#' Semi-supervised fairness estimation.
#'
#' @param Y Outcome, can contain missing values.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.
#' @param X Optional covariates matrix to be adjusted in the semi-supervised setting.
#' Default is NULL.
#' @param basis Different ways to construct W. Default is "ridge(S + X)". Other
#' options are "Poly(S)", "Poly(S) + X", "Spline(S)", "Spline(S) + X", "glm(S)",
#' "glm(S + X)", "ridge(S)", "ridge(S + X)"
#' @param nknots Number of knots (only used when basis = "Spline(S)" or
#' "Spline(S) + X"). Default is 3. Ignored otherwise.
#' @return List of estimated fairness metrics and their variances.
#' @export
#'

SSFairness <- function(Y,
                       S,
                       A,
                       threshold = 0.5,
                       X = NULL,
                       basis = "Poly(S)",
                       nknots = 3,
                       W = NULL) {
  if (!is.null(W)) {
    W_label <- 4 * rbeta(sum(!is.na(Y)), 1 / 2, 3 / 2)
  } else {
    W_label <- NULL
  }

  labeled_ind <- which(!is.na(Y))

  # Create Indicator
  C <- ifelse(S > threshold, 1, 0)

  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  C_labeled <- C[labeled_ind]


  # Unlabeled data.
  S_unlabeled <- S[-labeled_ind]
  A_unlabeled <- A[-labeled_ind]
  C_unlabeled <- C[-labeled_ind]

  if (!is.null(X)) {
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
  }

  nclass <- sort(unique(A))

  # Polynomial basis with selection
  if (basis == "Poly(S)") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha_glm(
        Y = Y_labeled[A_labeled == a],
        covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
        additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            coef = log(ncol(basis_labeled) + 1),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, -1],
          C_unlabeled[A_unlabeled == a]
        )
      ) %*% gamma)

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a, -1],
          C_labeled[A_labeled == a]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "Poly(S) + X") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha_glm(
        Y = Y_labeled[A_labeled == a],
        covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
        additional_matrix = cbind(
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        ) %>% as.matrix()
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            coef = log(ncol(basis_labeled) + ncol(X_labeled) + 1),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, -1],
          C_unlabeled[A_unlabeled == a],
          X_unlabeled[A_unlabeled == a, ]
        )
      ) %*% gamma)

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a, -1],
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "Spline(S)") {
    alphas <- c(1, 1)

    # Natural cubic spline basis with S included
    basis_exp <- NaturalSplineBasis(c(S_labeled, S_unlabeled),
      num_knots = nknots
    )

    # Labeled Basis
    basis_labeled <- basis_exp[1:length(S_labeled), ]

    # Unlabeled Basis
    basis_unlabeled <- basis_exp[(length(S_labeled) + 1):nrow(basis_exp), ]

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values

    nclass <- sort(unique(A))
    for (a in nclass) {
      ## Step I, Ridge Regression
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(basis_labeled[A_labeled == a, ], C_labeled[A_labeled == a]),
            y = Y_labeled[A_labeled == a],
            coef = log(ncol(basis_exp) + 1),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          return(NA)
        }
      )

      index_vector <- A_unlabeled == a

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, ],
          C_unlabeled[index_vector]
        )
      ) %*% gamma)
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a, ],
          C_labeled[A_labeled == a]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "Spline(S) + X") {
    alphas <- c(1, 1)

    # Natural cubic spline basis with S included
    basis_exp <- NaturalSplineBasis(c(S_labeled, S_unlabeled),
      num_knots = nknots
    )

    # Labeled Basis
    basis_labeled <- basis_exp[1:length(S_labeled), ]

    # Unlabeled Basis
    basis_unlabeled <- basis_exp[(length(S_labeled) + 1):nrow(basis_exp), ]

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values

    nclass <- sort(unique(A))
    for (a in nclass) {
      ## Step I, Ridge Regression
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, ], C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ),
            y = Y_labeled[A_labeled == a],
            coef = log(ncol(basis_exp) + 1),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          return(NA)
        }
      )

      index_vector <- A_unlabeled == a

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, ],
          C_unlabeled[index_vector], X_unlabeled[A_unlabeled == a, ]
        )
      ) %*% gamma)
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a, ],
          C_labeled[A_labeled == a], X_labeled[A_labeled == a, ]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "glm(S)") {
    # Labeled Basis
    basis_labeled <- S_labeled

    # Unlabeled Basis
    basis_unlabeled <- S_unlabeled

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c(1, 1)

    for (a in nclass) {
      # Just glm
      gamma <- glm.fit(
        x = cbind(1, S_labeled[A_labeled == a], C_labeled[A_labeled == a]),
        y = Y_labeled[A_labeled == a],
        weights = W_label[A_labeled == a],
        family = binomial(link = "logit")
      )$coefficients

      index_vector <- A_unlabeled == a

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, S_unlabeled[A_unlabeled == a],
          C_unlabeled[index_vector]
        )
      ) %*% gamma)
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, S_labeled[A_labeled == a],
          C_labeled[A_labeled == a]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "glm(S + X)") {
    alphas <- c(1, 1)

    # Labeled Basis
    basis_labeled <- S_labeled

    # Unlabeled Basis
    basis_unlabeled <- S_unlabeled

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values

    for (a in nclass) {
      # Just glm
      gamma <- glm.fit(
        x = cbind(
          1, S_labeled[A_labeled == a],
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        ),
        y = Y_labeled[A_labeled == a],
        weights = W_label[A_labeled == a],
        family = binomial(link = "logit")
      )$coefficients

      index_vector <- A_unlabeled == a

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, S_unlabeled[A_unlabeled == a],
          C_unlabeled[index_vector],
          X_unlabeled[A_unlabeled == a, ]
        )
      ) %*% gamma)
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, S_labeled[A_labeled == a],
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "ridge(S)") {
    # Labeled Basis
    basis_labeled <- S_labeled

    # Unlabeled Basis
    basis_unlabeled <- S_unlabeled

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c(1, 1)

    for (a in nclass) {
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a],
              C_labeled[A_labeled == a]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            coef = log(2),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a],
          C_unlabeled[A_unlabeled == a]
        )
      ) %*% gamma)

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a],
          C_labeled[A_labeled == a]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  } else if (basis == "ridge(S + X)") {
    # Labeled Basis
    basis_labeled <- S_labeled

    # Unlabeled Basis
    basis_unlabeled <- S_unlabeled

    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c(1, 1)

    for (a in nclass) {
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            coef = log(1 + ncol(X_labeled) + 1),
            weights = W_label[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a],
          C_unlabeled[A_unlabeled == a],
          X_unlabeled[A_unlabeled == a, ]
        )
      ) %*% gamma)

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_labeled[A_labeled == a],
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        )
      ) %*% gamma)
      m_labeled[A_labeled == a] <- imputed_labeled
    }
  }

  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = W
  )

  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold,
    method = "semi-supervised"
  )

  return(list(est = est, var = var))
}
