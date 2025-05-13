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
#' @param basis Different ways to construct W. Options are "Poly(S)", "Poly(S, X)"
#' @return List of estimated fairness metrics and their variances.
#' @export
#'

SSFairness <- function(Y,
                       S,
                       A,
                       threshold = 0.5,
                       X = NULL,
                       basis = "Poly(S)") {

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
  if (basis == "Beta") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      S_calibrated <- FitParametricCalibration(Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_unlabeled,
        method = "Beta"
      )
      m_unlabeled[A_unlabeled == a] <- S_calibrated[A_unlabeled == a]
    }
  }

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
            coef = log(ncol(basis_labeled) + 1)
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

  } else if (basis == "Poly(S, X)") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values

    for (a in nclass) {
      alpha <- find_alpha_glm(
        Y = Y_labeled[A_labeled == a],
        covariates_matrix = cbind(
          S_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        ) %>% as.matrix(),
        additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
      )
      # print(alpha)
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(cbind(S_labeled, X_labeled) %>% 
                                    as.data.frame(), alpha)
      basis_unlabeled <- polynomial(cbind(S_unlabeled, X_unlabeled) %>% 
                                      as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            coef = log(ncol(basis_labeled) + 1)
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

  } 
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold
  )

  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold,
    method = "semi-supervised"
  )

  return(list(est = est, var = var))
}
