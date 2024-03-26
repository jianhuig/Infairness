# Purpose: Estimation of fairness metrics in semi-supervised setting.
# Updated: 2022-08-30

#' Semi-supervised fairness estimation.
#'
#' @param Y Outcome, can contain missing values.
#' @param S Model score.
#' @param A Group indicator.
#' @param labeled_ind Index for labeled data.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.
#' @param augment_score Indicator to apply basis expansion to the score.
#' Default value is TRUE.
#' @param nknots Number of knots for basis expansion.
#' Default value is 3.
#' @param ridge Add Ridge penalty for basis expansion.
#' Default is TRUE.
#' @param lambda Ridge penalty parameter vector.
#' @param W Weight for labeled data.
#' Default value is NULL.
#' @param ridge Whether to perform ridge regression
#' @export
#'
Infairness <- function(Y,
                       S,
                       A,
                       threshold = 0.5,
                       nknots = 3,
                       penalty = 1.5,
                       W = NULL,
                       ridge = TRUE) {
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
        # One step ridge
        # Tuning parameter; log(k)/n^penalty
        lambda <- log(ncol(basis_exp)) / length(Y_labeled[A == a])^penalty

        ## Step I, Ridge Regression
        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(basis_labeled[A_labeled == a, ], C_labeled[A_labeled == a]),
              y = Y_labeled[A_labeled == a],
              lambda = lambda,
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
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = W
  )
  var <- Influence_curve(est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold, method = "semi-supervised")

  return(list(est = est, var = var))
}
