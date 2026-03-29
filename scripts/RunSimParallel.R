library(doParallel)
library(dplyr)
library(Infairness)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(dplyr)
  library(Infairness)
  library(glmnet)
})

args <- commandArgs(TRUE)
# args = c("balanced", "correct")
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
print(model)

set.seed(1234)
indep <- DataGeneration(
  n_labeled = 3000,
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
clusterEvalQ(cl, source(file.path("scripts", "TwoStep.R")))

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
  oracle <- SupervisedFairness(
    Y = dat$Y,
    S = dat$S,
    A = dat$A,
    threshold = threshold
  )$est

  # sup
  sup <- SupervisedFairness(
    Y = labeled$Y,
    S = labeled$S,
    A = labeled$A,
    threshold = threshold
  )$est
  # naive
  # naive <- SupervisedFairness(
  #   Y = dat$S,
  #   S = dat$S,
  #   A = dat$A,
  #   threshold = threshold
  # )


  # GLM(S)
  ss_s <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = NULL,
    basis = "glm(S)"
  )$est

  ss_s_all <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = dat %>% select(starts_with("X")) %>% as.matrix(),
    basis = "glm(S + X)"
  )$est

  ss_s_ridge <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = NULL,
    basis = "ridge(S)"
  )$est

  ss_s_all_ridge <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = dat %>% select(starts_with("X")) %>% as.matrix(),
    basis = "ridge(S + X)"
  )$est

  # GLM(Spline(S))
  # ss_spline <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "Spline(S)"
  # )
  # 
  # ss_spline_all <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(starts_with("X")) %>% as.matrix(),
  #   basis = "Spline(S) + X"
  # )

  # scale all x to 0 and 1
  # dat <- dat %>%
  #  mutate(across(starts_with("X"), ~ (.- min(.)) / (max(.) - min(.))))

  # Poly(S)
  ss_poly <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = NULL,
    basis = "Poly(S)"
  )$est

  # Poly(S) + X
  # ss_poly_x <- Infairness(Y = dat$Y_miss,
  #                         S = dat$S,
  #                         A = dat$A,
  #                         X = dat %>% select(X_1) %>% as.matrix(),
  #                         basis = "Poly(S) + X")


  # Poly(S) + all
  ss_poly_all <- Infairness(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = dat %>% select(starts_with("X")) %>% as.matrix(),
    basis = "Poly(S) + X"
  )$est


  # ss_platt <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "Platt"
  # )
  # 
  # # Beta Calibration
  # ss_Beta <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "Beta"
  # )

  ss_twostep_s <- Infairness_twostep(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = NULL,
    basis = "S"
  )

  ss_twostep_all <- Infairness_twostep(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = dat %>% select(contains("X")),
    basis = "S,X"
  )

  ss_twostep_x <- Infairness_twostep(
    Y = dat$Y_miss,
    S = dat$S,
    A = dat$A,
    X = dat %>% select(contains("X")),
    basis = "X"
  )


  list(
    oracle = oracle,
    sup = sup,
    #naive = naive,
    ss_s = ss_s,
    ss_s_all = ss_s_all,
    ss_s_ridge = ss_s_ridge,
    ss_s_all_ridge = ss_s_all_ridge,
    #ss_spline = ss_spline,
    #ss_spline_all = ss_spline_all,
    ss_poly = ss_poly,
    # ss_poly_x = ss_poly_x,
    ss_poly_all = ss_poly_all,
    #ss_platt = ss_platt,
    #ss_Beta = ss_Beta,
    ss_twostep_s = ss_twostep_s,
    ss_twostep_all = ss_twostep_all,
    ss_twostep_x = ss_twostep_x
  )
  
})
stopCluster(cl)
saveRDS(result, file = paste0("0302_", args[1], "_", args[2], ".rds"))
