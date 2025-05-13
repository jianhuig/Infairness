# Reliable fairness auditing with semi-supervised inference


# Installation

```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jianhuig/SS-Fairness-Audit")
```

# Example
```{R}
library(SSFairnessAudit)
library(dplyr)

set.seed(123)

dat <- DataGeneration(n_labeled = 1e3,
                      N_unlabeled = 1e4,
                      prot_att_prevalence = 0.5,
                      model = "misspecified 1",
                      rho = 0.4)
                      

# Indepedent Training Dataset
indep <- DataGeneration(n_labeled = 3e3,
                      N_unlabeled = 0,
                      prot_att_prevalence = 0.5,
                      model = "misspecified 1",
                      rho = 0.4)
                      
# Train Logistic Regression Model
model_0 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 0) %>% select(-Y_miss, -A, -W_1:-W_5))
model_1 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 1) %>% select(-Y_miss, -A, -W_1:-W_5))
dat$S <- rep(NA, nrow(dat))
dat$S[dat$A == 0] <- predict(model_0, newdata = dat %>% filter(A == 0), type = "response")
dat$S[dat$A == 1] <- predict(model_1, newdata = dat %>% filter(A == 1), type = "response")
                      
# Supervised Estimation                      
sup <- Audit_Fairness(Y = dat$Y_miss,
                      S = dat$S,
                      A = dat$A,
                      threshold = 0.5,
                      method = "supervised")

# Semi-supervised Estimation
ss <- Audit_Fairness(Y = dat$Y_miss,
                     S = dat$S,
                     A = dat$A,
                     threshold = 0.5,
                     method = "semi-supervised",
                     X = dat %>% select(contains("W")))
                      

# Point Estimate
                      
sup$est
  Metric     Group0     Group1       Delta
1    TPR 0.36290323 0.52941176 -0.16650854
2    FPR 0.04846939 0.06027397 -0.01180458
3    PPV 0.70312500 0.74117647 -0.03805147
4    NPV 0.82522124 0.85964912 -0.03442788
5     F1 0.47872340 0.61764706 -0.13892365
6    ACC 0.81007752 0.83884298 -0.02876546
7     BS 0.13250175 0.11684682  0.01565493

ss$est
  Metric     Group0    Group1       Delta
1    TPR 0.40404109 0.51458165 -0.11054056
2    FPR 0.05940048 0.04898599  0.01041449
3    PPV 0.69860671 0.77169637 -0.07308966
4    NPV 0.82242849 0.85892900 -0.03650051
5     F1 0.51197820 0.61744162 -0.10546342
6    ACC 0.80422735 0.84476843 -0.04054108
7     BS 0.13552039 0.11593314  0.01958725

# Variance Estimate
sup$var
  Metric       Group0       Group1        Delta
1    TPR 1.864552e-03 2.093571e-03 0.0039581232
2    FPR 1.176533e-04 1.551809e-04 0.0002728342
3    PPV 3.261566e-03 2.256870e-03 0.0055184357
4    NPV 3.190955e-04 3.023872e-04 0.0006214827
5     F1 2.019311e-03 1.600271e-03 0.0036195822
6    ACC 2.981627e-04 2.793088e-04 0.0005774714
7     BS 9.944322e-05 8.947343e-05 0.0001889167

ss$var
  Metric       Group0       Group1        Delta
1    TPR 8.662682e-04 1.039611e-03 0.0019058788
2    FPR 7.273322e-05 1.038688e-04 0.0001766020
3    PPV 3.031006e-03 2.102220e-03 0.0051332264
4    NPV 2.671569e-04 2.643316e-04 0.0005314885
5     F1 1.182558e-03 1.083128e-03 0.0022656857
6    ACC 2.516233e-04 2.444778e-04 0.0004961012
7     BS 8.444773e-05 7.622762e-05 0.0001606753
```
