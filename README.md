# SSFairnessAudit

R package for semi-supervised group fairness auditing when labels are only
partially observed.

## Installation

```r
devtools::install_github("jianhuig/Infairness")
```

## Example

The example below compares a supervised analysis using only labeled
observations with a semi-supervised analysis that uses all scores and a
binary covariate `W` through the `"Spline Interaction"` basis.

```r
library(SSFairnessAudit)

## Set up sample sizes.
n_labeled <- 200
N_unlabeled <- 10000
N_total <- n_labeled + N_unlabeled

generate_example_data <- function(N) {
  A <- rbinom(N, 1, 0.6)
  W <- rbinom(N, 1, 0.5)
  Y <- rbinom(N, 1, 0.3)
  S <- numeric(N)

  beta_par <- data.frame(
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    W = c(0, 1, 0, 1, 0, 1, 0, 1),
    Y = c(1, 1, 0, 0, 1, 1, 0, 0),
    alpha = c(14, 4.5, 2.5, 7, 7, 7, 5, 5),
    beta = c(4, 10, 10, 5, 5, 5, 7, 7)
  )

  for (j in seq_len(nrow(beta_par))) {
    row <- beta_par[j, ]
    idx <- A == row$A & W == row$W & Y == row$Y
    S[idx] <- rbeta(sum(idx), row$alpha, row$beta)
  }

  data.frame(Y = Y, S = S, A = A, W = W)
}

dat <- generate_example_data(N_total)

threshold <- 0.5

## Keep labels only for a random labeled subset.
labeled_ind <- sample(seq_len(N_total), n_labeled)
Y <- rep(NA_real_, N_total)
Y[labeled_ind] <- dat$Y[labeled_ind]
X <- data.frame(W = dat$W)

## Supervised fairness audit using labeled observations only.
fair_supervised <- SupervisedFairness(
  Y = dat$Y[labeled_ind],
  S = dat$S[labeled_ind],
  A = dat$A[labeled_ind],
  threshold = threshold
)

## Semi-supervised fairness audit using spline-by-W interactions.
fair_ss <- SSFairness(
  Y = Y,
  S = dat$S,
  A = dat$A,
  X = X,
  threshold = threshold,
  basis = "Spline Interaction",
  nknots = 3
)

fair_supervised$est
fair_ss$est
```
