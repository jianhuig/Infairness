# Purpose: Imputation-quality helpers for semi-supervised estimation.

compute_imputation_quality <- function(Y,
                                       S,
                                       A,
                                       threshold = 0.5,
                                       X = NULL,
                                       basis = c("Poly(S)", "Poly(S)"),
                                       k = 10,
                                       alphas = NULL,
                                       ...) {
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
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }

  m_unlabeled <- S_unlabeled
  m_labeled <- S_labeled
  alphas <- c()

  for (a in nclass) {
    basis_a <- basis[which(nclass == a)]

    if (basis_a == "Poly(S)") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
          )
        },
        error = function(e) 1
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
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a]
          )
        ) %*% gamma$coefficients
      )

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      C_a <- C_labeled[A_labeled == a]
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)

      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)
        alpha_fold <- alphas[a + 1]
        basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
        basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(basis_train[, -1], C_a[train_id]) %>% as.matrix(),
              y = Y_a[train_id]
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

        m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
          as.matrix(cbind(1, basis_test[, -1], C_a[test_id])) %*% gamma
        )
      }
    } else if (basis_a == "Poly(S) + X") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix()
          )
        },
        error = function(e) 1
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
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      X_a <- X_labeled[A_labeled == a, , drop = FALSE]
      C_a <- C_labeled[A_labeled == a]
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)

      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)
        alpha_fold <- alphas[a + 1]
        basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
        basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id],
                X_a[train_id, ]
              ) %>% as.matrix(),
              y = Y_a[train_id]
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

        m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
          as.matrix(
            cbind(1, basis_test[, -1], C_a[test_id], X_a[test_id, ])
          ) %*% gamma
        )
      }
    } else if (basis_a == "Interaction") {
      X_int <- S * as.matrix(X)
      X_int <- cbind(X, X_int)
      X_labeled_int <- X_int[labeled_ind, , drop = FALSE]
      X_unlabeled_int <- X_int[-labeled_ind, , drop = FALSE]

      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled_int[A_labeled == a, ]
            ) %>% as.matrix()
          )
        },
        error = function(e) 1
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
              X_labeled_int[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled_int[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      X_a <- X_labeled_int[A_labeled == a, , drop = FALSE]
      C_a <- C_labeled[A_labeled == a]
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)

      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)
        alpha_fold <- tryCatch(
          {
            find_alpha_glm(
              Y = Y_a[train_id],
              covariates_matrix = as.matrix(S_a[train_id]),
              additional_matrix = as.matrix(cbind(X_a[train_id, ], C_a[train_id]))
            )
          },
          error = function(e) 1
        )

        basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
        basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id],
                X_a[train_id, ]
              ) %>% as.matrix(),
              y = Y_a[train_id],
              weights = NULL
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

        m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
          as.matrix(
            cbind(1, basis_test[, -1], C_a[test_id], X_a[test_id, ])
          ) %*% gamma
        )
      }
    } else if (basis_a == "Beta") {
      alphas <- c(alphas, NA_real_)

      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_unlabeled,
        method = "Beta"
      )
      m_unlabeled[A_unlabeled == a] <- S_calibrated[A_unlabeled == a]

      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_labeled,
        method = "Beta"
      )
      m_labeled[A_labeled == a] <- S_calibrated[A_labeled == a]
    } else if (basis_a == "kernel") {
      alphas <- c(alphas, NA_real_)

      bandwidth <- sqrt(var(S_labeled[A_labeled == a])) /
        (length(Y_labeled[A_labeled == a])^0.45)
      mhat <- npreg(
        S_labeled[A_labeled == a],
        Y_labeled[A_labeled == a],
        S_unlabeled,
        bandwidth
      )
      m_unlabeled[A_unlabeled == a] <- mhat[A_unlabeled == a]

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)
        bandwidth_fold <- sqrt(var(S_a[train_id])) / (length(Y_a[train_id])^0.45)
        mhat_fold <- npreg(S_a[train_id], Y_a[train_id], S_a[test_id], bandwidth_fold)
        m_labeled[which(A_labeled == a)[test_id]] <- mhat_fold
      }
    } else {
      stop("Basis not recognized.")
    }
  }
  
  # Goodness of fit
  bs <- c()
  log_loss <- c()
  auc <- c()
  ece <- c()
  for (a in nclass) {
    # Brier score
    bs_a <- mean((Y_labeled[A_labeled == a] - m_labeled[A_labeled == a])^2)
    bs <- c(bs, bs_a)
    # Log-loss
    log_loss_a <- -mean(
      Y_labeled[A_labeled == a] * log(m_labeled[A_labeled == a]) +
        (1 - Y_labeled[A_labeled == a]) * log(1 - m_labeled[A_labeled == a])
    )
    log_loss <- c(log_loss, log_loss_a)
    # AUC
    auc_a <- pROC::auc(
      Y_labeled[A_labeled == a],
      m_labeled[A_labeled == a]
    ) %>% as.numeric()
    auc <- c(auc, auc_a)
    # ECE (use quantile bins; last bin right-closed)
    p_a <- m_labeled[A_labeled == a]
    y_a <- Y_labeled[A_labeled == a]
    bins <- unique(quantile(p_a, probs = seq(0, 1, by = 0.1), na.rm = TRUE))

    ece_a <- 0
    for (b in 1:(length(bins) - 1)) {
      if (b == length(bins) - 1) {
        ind <- which(p_a >= bins[b] & p_a <= bins[b + 1])
      } else {
        ind <- which(p_a >= bins[b] & p_a < bins[b + 1])
      }
      if (length(ind) > 0) {
        bin_acc <- mean(y_a[ind])
        bin_conf <- mean(p_a[ind])
        ece_a <- ece_a + (length(ind) / length(p_a)) * abs(bin_acc - bin_conf)
      }
    }
    ece <- c(ece, ece_a)
  }
  metrics_by_group <- cbind(
    Brier_Score = bs,
    Log_Loss = log_loss,
    AUC = auc,
    ECE = ece
  ) %>% as.data.frame()

  rownames(metrics_by_group) <- paste0("Group_", nclass)
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = NULL
  )
  
  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold,
    method = "semi-supervised"
  )
  
  return(list(est = est, var = var, alpha = alphas, imp_quality = metrics_by_group,
    m_labeled = m_labeled, m_unlabeled = m_unlabeled))
}

#' Imputation Quality
#'
#' Convenience wrapper around the shared semi-supervised imputation routine.
#' For new code, prefer `SSFairness(..., return_imputation_quality = TRUE)`.
#' @export
ImputeQuality <- function(Y,
                          S,
                          A,
                          threshold = 0.5,
                          X = NULL,
                          basis = c("Poly(S)", "Poly(S)"),
                          k = 10,
                          alphas = NULL,
                          ...) {
  compute_imputation_quality(
    Y = Y,
    S = S,
    A = A,
    threshold = threshold,
    X = X,
    basis = basis,
    k = k,
    alphas = alphas,
    ...
  )
}
