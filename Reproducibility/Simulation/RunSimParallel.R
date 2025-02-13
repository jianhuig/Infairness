library(doParallel)
library(dplyr)
library(SSFairnessAudit)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(dplyr)
  library(SSFairnessAudit)
})

args <- commandArgs(TRUE)
if (args[1] == "balanced") {
  prev <- 0.5 # Prevalence of the protected attribute.
} else {
  prev <- 0.4 # Prevalence of the protected attribute.
}
print(prev)

nclass <- 2

n <- 500 * nclass # Labeled data size.
N <- 1e4 * nclass # Unlabeled data size.
p <- 10
rho <- 0.4
threshold <- 0.5
model <- args[2]
set.seed(1234)
indep <- DataGeneration(
  n_labeled = 1e3,
  N_unlabeled = 0,
  prot_att_prevalence = prev,
  model = model,
  p = 10,
  rho = rho
)

model_0 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 0) %>% select(-Y_miss, -A))
model_1 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 1) %>% select(-Y_miss, -A))

clusterExport(cl = cl, varlist = list(
  "indep", "nclass", "n", "N",
  "p", "rho", "threshold", "model", "prev", "model_0", "model_1"
))

nsim <- 1e4
result <- pbapply::pblapply(cl = cl, X = 1:nsim, FUN = function(j) {
  # generate the main dataset
  dat <- DataGeneration(
    n_labeled = n,
    N_unlabeled = N,
    prot_att_prevalence = prev,
    model = model,
    p = 10,
    rho = rho
  )

  # using indepdent data to train S
  dat$S <- rep(NA, nrow(dat))
  # for (a in unique(indep$A)) {
  #   model_S <- glm(Y ~ .,
  #                  family = binomial(),
  #                  data = indep %>% filter(A == a) %>% select(-Y_miss, -A))
  #
  #   # Predict S for the corresponding group in dat
  #   dat$S[dat$A == a] <- predict(model_S, newdata = dat %>% filter(A == a), type = "response")
  # }
  dat$S[dat$A == 0] <- predict(model_0, newdata = dat %>% filter(A == 0), type = "response")
  dat$S[dat$A == 1] <- predict(model_1, newdata = dat %>% filter(A == 1), type = "response")

  # prepare main data
  dat$C <- ifelse(dat$S > threshold, 1, 0)
  labeled <- dat %>% filter(!is.na(Y_miss))
  unlabeled <- dat %>% filter(is.na(Y_miss))

  # check AUC
  # pROC::auc(dat$Y, dat$S)
  # pROC::auc(dat$Y[dat$A==0], dat$S[dat$A==0])
  # pROC::auc(dat$Y[dat$A==1], dat$S[dat$A==1])

  # oracle
  oracle <- Audit_Fairness(
    Y = dat$Y,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "supervised"
  )

  # sup
  sup <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "supervised"
  )

  # naive
  naive <- Audit_Fairness(
    Y = dat$S,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "supervised"
  )


  # GLM(S)
  ss_s <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "glm(S)",
    X = NULL
  )

  ss_s_all <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "glm(S + X)",
    X = dat %>% select(starts_with("X"))
  )

  ss_s_ridge <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "ridge(S)",
    X = NULL
  )

  ss_s_all_ridge <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "ridge(S + X)",
    X = dat %>% select(starts_with("X"))
  )

  # GLM(Spline(S))
  ss_spline <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "Spline(S)",
    X = NULL
  )

  ss_spline_all <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "Spline(S) + X",
    X = dat %>% select(starts_with("X")) %>% as.matrix()
  )


  # Poly(S)
  ss_poly <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "Poly(S)",
    X = NULL
  )


  # Poly(S) + all
  ss_poly_all <- Audit_Fairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    threshold = threshold,
    method = "semi-supervised",
    basis = "Poly(S) + X",
    X = dat %>% select(starts_with("X")) %>% as.matrix()
  )


  list(
    oracle = oracle,
    sup = sup,
    naive = naive,
    ss_s = ss_s,
    ss_s_all = ss_s_all,
    ss_s_ridge = ss_s_ridge,
    ss_s_all_ridge = ss_s_all_ridge,
    ss_spline = ss_spline,
    ss_spline_all = ss_spline_all,
    ss_poly = ss_poly,
    ss_poly_all = ss_poly_all
  )
})

stopCluster(cl)
saveRDS(result, file = paste0("0101_", args[1], "_", args[2], ".rds"))
