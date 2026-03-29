# Purpose: Estimation of fairness metrics in semi-supervised setting.

#' Imputation Quality
#' @export


ImputeQuality <- function(Y,
                          S,
                          A,
                          threshold = 0.5,
                          X = NULL,
                          basis = "Poly(S)",
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

  if (basis == "Poly(S)") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values
    
    for (a in nclass) {
      
      alpha <- tryCatch({find_alpha_glm(
        Y = Y_labeled[A_labeled == a],
        covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
        additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
      )}, error = function(e) { 1 })
      
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
      
      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, -1],
          C_unlabeled[A_unlabeled == a]
        )
      ) %*% gamma$coefficients)
      
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
    }

    for (a in nclass) {
      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      C_a <- C_labeled[A_labeled == a]

      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)

      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)

        alpha_fold <- find_alpha_glm(
          Y = Y_a[train_id],
          covariates_matrix = S_a[train_id] %>% as.matrix(),
          additional_matrix = C_a[train_id] %>% as.matrix()
        )
        alpha_fold <- alphas[a + 1] # same alpha as above

        basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
        basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id]
              ) %>% as.matrix(),
              y = Y_a[train_id]
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients
        # gamma <- coef(glmnet(
        #   y = Y_a[train_id], x = cbind(basis_train[, -1],
        #                         C_a[train_id]),
        #   family = "binomial",
        #   alpha = 0,
        #   lambda = gamma$best_lambda
        # ))

        imputed_test <- boot::inv.logit(as.matrix(
          cbind(
            1, basis_test[, -1],
            C_a[test_id]
          )
        ) %*% gamma)

        m_labeled[ which(A_labeled == a)[test_id] ] <- imputed_test
      }
    }
    
  } else if (basis == "Poly(S) + X") {
    # Augmentation in each class
    m_unlabeled <- S_unlabeled # imputed values
    m_labeled <- S_labeled # imputed values
    alphas <- c() # store alpha values
  
    for (a in nclass) {
      alpha <- tryCatch({find_alpha_glm(
        Y = Y_labeled[A_labeled == a],
        covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
        additional_matrix = cbind(
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        ) %>% as.matrix()
      )}, error = function(e) { 1 })
  
      # print(alpha)
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
      
      imputed_unlabeled <- boot::inv.logit(as.matrix(
        cbind(
          1, basis_unlabeled[A_unlabeled == a, -1],
          C_unlabeled[A_unlabeled == a],
          X_unlabeled[A_unlabeled == a, ]
        )
      ) %*% gamma$coefficients)
      
      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
    }

    for (a in nclass) {
      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      X_a <- X_labeled[A_labeled == a, , drop = FALSE]
      C_a <- C_labeled[A_labeled == a]
      
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)

        alpha_fold <- tryCatch({
          find_alpha_glm(
            Y = Y_a[train_id],
            covariates_matrix = as.matrix(S_a[train_id]),
            additional_matrix = as.matrix(cbind(X_a[train_id, ], C_a[train_id]))
          )
        }, error = function(e) { 1 })

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

        imputed_test <- boot::inv.logit(as.matrix(
          cbind(
            1, basis_test[, -1],
            C_a[test_id],
            X_a[test_id, ]
          )
        ) %*% gamma)
        m_labeled[ which(A_labeled == a)[test_id] ] <- imputed_test
      }
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
