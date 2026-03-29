# Purpose: Internal helpers for model-selection workflows.

step_two_fit <- function(Y,
                         S,
                         A,
                         threshold = 0.5,
                         basis1 = NULL,
                         basis2 = NULL,
                         k = 10) {
  labeled_ind <- which(!is.na(Y))

  Y_labeled <- Y[labeled_ind]
  A_labeled <- A[labeled_ind]
  S_labeled <- S[labeled_ind]

  A_unlabeled <- A[-labeled_ind]
  S_unlabeled <- S[-labeled_ind]

  nclass <- sort(unique(A))

  m_unlabeled <- rep(NULL, length(A_unlabeled))
  m_labeled <- rep(NULL, length(A_labeled))

  g1 <- nclass[1]
  g2 <- nclass[2]

  for (a in nclass) {
    if (a == g1) {
      X <- as.matrix(basis1)
    } else {
      X <- as.matrix(basis2)
    }

    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]

    gamma <- tryCatch(
      {
        RidgeRegression(
          X = X_labeled[A_labeled == a, -1],
          y = Y_labeled[A_labeled == a],
          weights = NULL
        )
      },
      error = function(e) {
        print("Ridge Regression produced an error")
        print(e)
      }
    )

    imputed_unlabeled <- boot::inv.logit(
      as.matrix(X_unlabeled[A_unlabeled == a, ]) %*% gamma$coefficients
    )

    m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
  }

  for (a in nclass) {
    if (a == g1) {
      X <- as.matrix(basis1)
    } else {
      X <- as.matrix(basis2)
    }

    X_labeled <- X[labeled_ind, , drop = FALSE]

    Y_a <- Y_labeled[A_labeled == a]
    X_labeled_a <- X_labeled[A_labeled == a, , drop = FALSE]

    fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)

    for (i in 1:k) {
      train_id <- which(fold != i)
      test_id <- which(fold == i)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = X_labeled_a[train_id, -1],
            y = Y_a[train_id],
            weights = NULL
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )$coefficients

      imputed_test <- boot::inv.logit(
        as.matrix(X_labeled_a[test_id, , drop = FALSE]) %*% gamma
      )

      m_labeled[which(A_labeled == a)[test_id]] <- imputed_test
    }
  }

  est <- get_metric(
    Y = m_unlabeled,
    S = S_unlabeled,
    A = A_unlabeled,
    threshold = threshold,
    W = NULL
  )

  var <- Influence_curve(
    est,
    Y_labeled,
    S_labeled,
    A_labeled,
    m_labeled,
    threshold,
    method = "semi-supervised"
  )

  return(list(est = est, var = var))
}
