StepII <- function(Y,
                   S,
         A,
         threshold = 0.5,
         basis1 = NULL,
         basis2 = NULL,
         k = 10) {
  labeled_ind <- which(!is.na(Y))
  
  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  A_labeled <- A[labeled_ind]
  S_labeled <- S[labeled_ind]
  
  # Unlabeled data.
  A_unlabeled <- A[-labeled_ind]
  S_unlabeled <- S[-labeled_ind]

  nclass <- sort(unique(A))
  
  # Augmentation in each class
  m_unlabeled <- rep(NULL, length(A_unlabeled))
  m_labeled <- rep(NULL, length(A_labeled))
  
  for (a in nclass) {
    if(a == 0){
      X = basis1
    } else {
      X = basis2
    }
    
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
    
    
    gamma <- tryCatch(
      {
        RidgeRegression(
          X = X_labeled[A_labeled == a, -1],
          y = Y_labeled[A_labeled == a],
          coef = log(ncol(X_labeled)),
          weights = NULL
        )
      },
      error = function(e) {
        print("Ridge Regression produced an error")
        print(e)
      }
    )
    
    imputed_unlabeled <- boot::inv.logit(as.matrix(X_unlabeled[A_unlabeled == a, ]
    ) %*% gamma$coefficients)
    
    m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
  }
  
  for (a in nclass) {
    if(a == 0){
      X = basis1
    } else {
      X = basis2
    }
    
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
    
    Y_a <- Y_labeled[A_labeled == a]
    X_labeled_a <- X_labeled[A_labeled == a, , drop = FALSE]
    X_unlabeled_a <- X_unlabeled[A_labeled == a, , drop = FALSE]
    
    fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
    
    for (i in 1:k) {
      train_id <- which(fold != i)
      test_id <- which(fold == i)
      
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = X_labeled_a[train_id, -1],
            y = Y_a[train_id],
            coef = log(ncol(X_labeled_a)),
            weights = NULL
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )$coefficients
    

      imputed_test <- boot::inv.logit(
        as.matrix(X_labeled_a[test_id, , drop = FALSE]
        )
       %*% gamma)
      
      m_labeled[which(A_labeled == a)[test_id] ] <- imputed_test
    }
  }
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = NULL
  )
  
  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold,
    method = "semi-supervised"
  )
  
  return(list(est = est, var = var, m_labeled = m_labeled, m_unlabeled = m_unlabeled))
}

