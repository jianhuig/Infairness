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
#' @param method Smoothing method. Options are  natural cubic spline "spline" , kernel smoothing "ks" , and quadratic polynomial "quad"
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
                       method = "spline",
                       ridge = TRUE,
                       step = 1) {
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

  if (method == "spline") {
    # Natural cubic spline basis with S included
    basis_exp <- NaturalSplineBasis(c(S_labeled, S_unlabeled),
      num_knots = nknots
    )

    # Labeled Basis
    basis_labeled <- basis_exp[1:length(S_labeled), ]

    # Unlabeled Basis
    basis_unlabeled <- basis_exp[(length(S_labeled) + 1):nrow(basis_exp), ]
  } else if (method == "ns") {
    # Natural cubic spline basis
    basis_exp <- splines::ns(c(S_labeled, S_unlabeled), df = 2)

    # Labeled Basis
    basis_labeled <- basis_exp[1:length(S_labeled), ]

    # Unlabeled Basis
    basis_unlabeled <- basis_exp[(length(S_labeled) + 1):nrow(basis_exp), ]
  } else if (method == "quad") {
    # Natural cubic spline basis
    basis_exp <- stats::poly(c(S_labeled, S_unlabeled), df = 2)

    # Labeled Basis
    basis_labeled <- basis_exp[1:length(S_labeled), ]

    # Unlabeled Basis
    basis_unlabeled <- basis_exp[(length(S_labeled) + 1):nrow(basis_exp), ]
  } else if (method == "ks") {
    # Transformed Score
    Strans <- ecdf(S)(S)
    Strans_labeled <- Strans[labeled_ind]
    Strans_unlabeled <- Strans[-labeled_ind]
  } else if (method == "beta") {
    S_labeled[S_labeled == 1] <- S_labeled[S_labeled == 1] - 1e6
    S_labeled[S_labeled == 0] <- S_labeled[S_labeled == 0] + 1e6
    S_unlabeled[S_unlabeled == 0] <- S_unlabeled[S_unlabeled == 0] + 1e6
    S_unlabeled[S_unlabeled == 0] <- S_unlabeled[S_unlabeled == 0] + 1e6
  }

  # Augmentation in each class
  m_unlabeled <- S_unlabeled # imputed values
  m_labeled <- S_labeled # imputed values

  nclass <- sort(unique(A))
  for (a in nclass) {
    if (method == "spline") {
      if (step == 2) {
        # 2 Step procedure

        if (ridge) {
          # Tuning parameter; log(k)/n^penalty
          lambda <- log(ncol(basis_exp)) / length(Y_labeled[A == a])^penalty

          ## Step I, Ridge Regression
          gamma <- tryCatch(
            {
              RidgeRegression(
                X = basis_labeled[A_labeled == a, ],
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
        } else {
          gamma <- coef(glm(
            Y_labeled[A_labeled == a] ~
              basis_labeled[A_labeled == a, ],
            weights = W_label[A_labeled == a],
            family = "binomial"
          ))
        }

        if (length(gamma) == 1) {
          return(NA)
        }
        recal_labeled <- as.matrix(
          cbind(1, basis_labeled[A_labeled == a, ])
        ) %*% gamma

        recal_unlabeled <- as.matrix(
          cbind(1, basis_unlabeled[A_unlabeled == a, ])
        ) %*% gamma

        ## Step II: Robust augmentation
        impute_model <- suppressWarnings(
          glm(Y_labeled[A_labeled == a] ~ C_labeled[A_labeled == a],
            offset = recal_labeled,
            family = binomial,
            weights = W_label[A_labeled == a]
          )$coefficients
        )

        imputed_labeled <- boot::inv.logit(
          impute_model[1] +
            impute_model[2] * C_labeled[A_labeled == a] +
            recal_labeled
        )

        imputed_unlabeled <- boot::inv.logit(
          impute_model[1] +
            impute_model[2] * C_unlabeled[A_unlabeled == a] +
            recal_unlabeled
        )

        m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
      } else {
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
    } else if (method == "quad") {
      gamma <- coef(glm(
        Y_labeled[A_labeled == a] ~
          basis_labeled[A_labeled == a, ],
        weights = W_label[A_labeled == a],
        family = "binomial"
      ))

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(as.matrix(
        cbind(1, basis_unlabeled[A_unlabeled == a, ])
      ) %*% gamma)
    } else if (method == "platt") {
      gamma <- coef(glm(
        Y_labeled[A_labeled == a] ~
          S_labeled[A_labeled == a],
        weights = W_label[A_labeled == a],
        family = "binomial"
      ))
      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(as.matrix(
        cbind(1, S_unlabeled[A_unlabeled == a])
      ) %*% gamma)
    } else if (method == "beta") {
      gamma <- coef(glm(
        Y_labeled[A_labeled == a] ~
          log(S_labeled[A_labeled == a]) +
          log(1 - S_labeled[A_labeled == a]),
        weights = W_label[A_labeled == a],
        family = "binomial"
      ))
      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(as.matrix(
        cbind(
          1, log(S_unlabeled[A_unlabeled == a]),
          log(1 - S_unlabeled[A_unlabeled == a])
        )
      ) %*% gamma)
    } else {
      # compute bandwith
      bw <- sd(Strans_labeled[A_labeled == a]) / (length(Strans_labeled[A_labeled == a])^0.45)

      # NW estimator
      mhat <- npreg(Strans_labeled[A_labeled == a], Y_labeled[A_labeled == a],
        Strans_unlabeled[A_unlabeled == a],
        bw = bw, Wt = W_label[A_labeled == a]
      )
      m_unlabeled[A_unlabeled == a] <- mhat
    }
  }

  return(get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = W
  ))
}
