library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(purrr)
library(ggsci)

files <- c("0101_unbalanced_correct.rds", "0101_unbalanced_misspecified 1.rds", "0101_unbalanced_misspecified 2.rds", "0101_unbalanced_misspecified 3.rds")

for(f in files){

result <- readRDS(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Data/", f))

# Create the oracle estimates and variances once (outside the loop)
oracle_est <- do.call(rbind, lapply(result, function(ll) ll$oracle$est)) %>%
  group_by(Metric) %>%
  summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -c(Metric),
    names_to = "Group",
    values_to = "Est"
  )

oracle_avar <- do.call(rbind, lapply(result, function(ll) ll$oracle$var)) %>%
  group_by(Metric) %>%
  summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -c(Metric),
    names_to = "Group",
    values_to = "aVar"
  )

oracle_evar <- do.call(rbind, lapply(result, function(ll) ll$oracle$est)) %>%
  group_by(Metric) %>%
  summarise_all(function(x) var(x, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -c(Metric),
    names_to = "Group",
    values_to = "eVar"
  )

oracle <- oracle_est %>%
  left_join(oracle_avar, by = c("Metric", "Group")) %>%
  left_join(oracle_evar, by = c("Metric", "Group"))


results_list <- list()

# Loop over each method to compute bias and MSE
for (method in names(result[[1]])) {
  
  # Create the estimation part (mean)
  est <- do.call(rbind, lapply(result, function(ll) ll[[method]]$est)) %>%
    group_by(Metric) %>%
    summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
    pivot_longer(
      cols = -c(Metric),
      names_to = "Group",
      values_to = "Est"
    )
  
  # Create the variance part (variance)
  avar <- do.call(rbind, lapply(result, function(ll) ll[[method]]$var)) %>%
    group_by(Metric) %>%
    summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
    pivot_longer(
      cols = -c(Metric),
      names_to = "Group",
      values_to = "aVar"
    )
  
  evar <- do.call(rbind, lapply(result, function(ll) ll[[method]]$est)) %>%
    group_by(Metric) %>%
    summarise_all(function(x) var(x, na.rm = TRUE)) %>%
    pivot_longer(
      cols = -c(Metric),
      names_to = "Group",
      values_to = "eVar"
    )
  
  # Combine estimation and variance
  combined <- est %>%
    left_join(avar, by = c("Metric", "Group")) %>%
    left_join(evar, by = c("Metric", "Group")) %>%
    mutate(Method = method)
  
  # Calculate the bias and MSE by comparing with oracle
  combined <- combined %>%
    left_join(oracle %>% select(Metric, Group, Est) %>% rename(Oracle_Est = Est), 
              by = c("Metric", "Group")) %>%
    mutate(bias = Est - Oracle_Est, MSE = bias^2 + eVar)
  
  # Append to the results list
  results_list[[method]] <- combined
}

# Combine all the results into one data frame
final_results <- bind_rows(results_list)

labels <- c("Oracle", "Supervised", "Naive", "glm(S)", "glm(S + X)", "ridge(S)", "ridge(S + X)", "Spline(S)", "Spline(S) + X", "Poly(S)", "Poly(S) + X", "Platt", "Beta")

final_results$Method <- factor(final_results$Method, levels = names(result[[1]]), labels = labels)

# add CI
final_results <- final_results %>%
  group_by(Metric, Group, Method) %>%
  mutate(Lower = Est - qnorm(0.975) * sqrt(aVar),
         Upper = Est + qnorm(0.975) * sqrt(aVar))

final_results %>% 
  filter(!Metric %in% c("FNR", "TNR")) %>%
  filter(Method %in% c("Oracle", "Supervised", "ridge(S + X)", "Naive", "Platt", "Beta")) %>%
  mutate(
    Method = factor(Method,
                    levels = c("Oracle", "Supervised", "ridge(S + X)", "Naive", "Platt", "Beta"),
                    labels = c("Oracle", "Supervised", "SS", "Naive", "Platt Scaling", "Beta Calibration")),
    Group = factor(Group, levels = c("Group0", "Group1", "Delta"))
  ) %>%
  ggplot(aes(x = Metric, y = Est, fill = Method, group = Method)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_fill_nejm() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank()    # Proper way to remove x-axis label
  ) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) +
  labs(y = "")
  
#ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/pest_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), ".png"))

final_results %>% 
  mutate(pbias = (Est - Oracle_Est) / Oracle_Est * 100) %>%
  filter(!Metric %in% c("FNR", "TNR")) %>%
  filter(Method %in% c("Supervised", "ridge(S + X)", "Naive", "Platt", "Beta")) %>%
  mutate(
    Method = factor(Method,
                    levels = c("Supervised", "ridge(S + X)", "Naive", "Platt", "Beta"),
                    labels = c("Supervised", "SS", "Naive", "Platt Scaling", "Beta Calibration")),
    Group = factor(Group, levels = c("Group0", "Group1", "Delta"))
  ) %>%
  ggplot(aes(x = Metric, y = pbias, fill = Method, group = Method)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  facet_wrap(~ Group) +
  theme_bw() +
  scale_fill_manual(values = pal_nejm()(6)[-1]) + # Use first 5 NEJM colors+
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank()    # Proper way to remove x-axis label
  ) +
  labs(y = "Percentage Bias") +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

#ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/pbias_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), ".png"))

sup_mse <- final_results %>%
  filter(Method == "Supervised") %>%
  select(Metric, Group, MSE) %>%
  rename(sup_MSE = MSE)

# Initialize an empty data frame to store RE results
re_final <- data.frame()

# Loop over each method to compute RE relative to sup_est
for (i in seq_along(methods)) {
  
  # Get the current method and label
  method <- methods[i]
  label <- labels[i]
  
  # Get MSE of the current method
  method_mse <- final_results %>%
    filter(Method == label) %>%
    select(Metric, Group, MSE)
  
  # Compute RE (sup_MSE / method_MSE)
  re_temp <- sup_mse %>%
    left_join(method_mse, by = c("Metric", "Group")) %>%
    mutate(RE = sup_MSE / MSE, 
           Method = label)  # Apply the label
  
  # Append to the re_final data frame
  re_final <- rbind(re_final, re_temp)
}


# Convert 'Method' to a factor with the correct order for plotting, if needed
re_final$Method <- factor(re_final$Method, levels = labels)
re_final <- re_final %>% mutate(Setting = sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f))

re_df <- rbind(re_df, re_final)

re_df %>%
  filter(Method == "ridge(S + X)") %>%
  filter(!Metric %in% c("FNR", "TNR")) %>%
  mutate(
    Group = factor(Group, levels = c("Group0", "Group1", "Delta"))
  ) %>%
  ggplot(aes(x = Metric, y = RE, fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge2()) +
  scale_fill_cosmic() +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "red") +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank()    # Proper way to remove x-axis label
  ) + 
  labs(y = "Relative Efficiency\n Supervised : Semi-Supervised")
#ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/re_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), ".png"))

re_df %>%
  filter(!Metric %in% c("TNR", "FNR")) %>%
  filter(Method %in% c("ridge(S + X)", "glm(S)")) %>%
  mutate(Group = factor(Group, levels = c("Group0", "Group1", "Delta"))) %>%
  filter(Group == "Delta") %>%
  mutate(Method = factor(Method, levels = c("ridge(S + X)", "glm(S)"), labels = c("Proposed Method", "Alternative Imputation"))) %>%
  mutate(Metric = case_when(
    Metric == "ACC" ~ "Delta[ACC]",
    Metric == "BS"  ~ "Delta[BS]",
    Metric == "F1"  ~ "Delta[F1]",
    Metric == "FNR" ~ "Delta[FNR]",
    Metric == "FPR" ~ "Delta[FPR]",
    Metric == "NPV" ~ "Delta[NPV]",
    Metric == "PPV" ~ "Delta[PPV]",
    Metric == "TNR" ~ "Delta[TNR]",
    Metric == "TPR" ~ "Delta[TPR]")) %>%
  ggplot(aes(x = Metric, y = RE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  labs(y = "Relative Efficiency\n Supervised : Semi-Supervised") +
  xlab("") +
  scale_x_discrete(labels = scales::parse_format())+
  theme(legend.position = "bottom") +
  scale_fill_npg()

ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/re_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), "noX.png"))


re_df %>% 
  filter(!Metric %in% c("TNR", "FNR")) %>%
  filter(Method %in% c("ridge(S + X)", "Poly(S) + X")) %>%
  mutate(Group = factor(Group, levels = c("Group0", "Group1", "Delta"))) %>%
  filter(Group == "Delta") %>%
  mutate(Method = factor(Method, levels = c("ridge(S + X)", "Poly(S) + X"), labels = c("Proposed Method", "Polynomial Basis"))) %>%
  mutate(Metric = case_when(
    Metric == "ACC" ~ "Delta[ACC]",
    Metric == "BS"  ~ "Delta[BS]",
    Metric == "F1"  ~ "Delta[F1]",
    Metric == "FNR" ~ "Delta[FNR]",
    Metric == "FPR" ~ "Delta[FPR]",
    Metric == "NPV" ~ "Delta[NPV]",
    Metric == "PPV" ~ "Delta[PPV]",
    Metric == "TNR" ~ "Delta[TNR]",
    Metric == "TPR" ~ "Delta[TPR]")) %>%
  ggplot(aes(x = Metric, y = RE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  labs(y = "Relative Efficiency\n Supervised : Semi-Supervised") +
  xlab("") +
  scale_x_discrete(labels = scales::parse_format())+
  theme(legend.position = "bottom") +
  scale_fill_npg()

ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/re_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), "poly.png"))

re_df %>% 
  filter(!Metric %in% c("TNR", "FNR")) %>%
  filter(Method %in% c("ridge(S + X)", "Spline(S) + X")) %>%
  mutate(Group = factor(Group, levels = c("Group0", "Group1", "Delta"))) %>%
  filter(Group == "Delta") %>%
  mutate(Method = factor(Method, levels = c("ridge(S + X)", "Spline(S) + X"), labels = c("Proposed Method", "Cubic Spline Basis"))) %>%
  mutate(Metric = case_when(
    Metric == "ACC" ~ "Delta[ACC]",
    Metric == "BS"  ~ "Delta[BS]",
    Metric == "F1"  ~ "Delta[F1]",
    Metric == "FNR" ~ "Delta[FNR]",
    Metric == "FPR" ~ "Delta[FPR]",
    Metric == "NPV" ~ "Delta[NPV]",
    Metric == "PPV" ~ "Delta[PPV]",
    Metric == "TNR" ~ "Delta[TNR]",
    Metric == "TPR" ~ "Delta[TPR]")) %>%
  ggplot(aes(x = Metric, y = RE, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  labs(y = "Relative Efficiency\n Supervised : Semi-Supervised") +
  xlab("") +
  scale_x_discrete(labels = scales::parse_format())+
  theme(legend.position = "bottom") +
  scale_fill_npg()

ggsave(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/Infairness_simulation/Plots/re_", sub("^[^_]*_[^_]*_(.*?)\\.rds$", "\\1", f), "spline.png"))

}