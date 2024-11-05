library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(ggpubr)

## read in all the fit results
fit_result_rds_paths <- list.files(
  path = "results", pattern = "*fit_results.rds*", full.names = TRUE,
  recursive = TRUE
)
fit_tbl_ls <- lapply(fit_result_rds_paths, readRDS)

## read in all the evaluation results
eval_result_rds_paths <- list.files(
  path = "results", pattern = "*eval_results.rds", full.names = TRUE,
  recursive = TRUE
)
eval_tbl_ls <- lapply(eval_result_rds_paths, readRDS)

## read in the monte-carlo-based oracle ITR mean outcomes
pop_mean_itr_df <- readRDS("results/population-mean-itr.Rds")

## define some labelling functions
unihtee_filtering_labeller <- function(x) ifelse(x == TRUE, "Yes", "No")
facet_labeller <- function(x) {
  case_when(
    x == "lm_sparse_identity" ~ "Sparse LM, Identity Cov.",
    x == "lm_non_sparse_identity" ~ "Non-Sparse LM, Identity Cov.",
    x == "nl_sparse_identity" ~ "Sparse NLM, Identity Cov.",
    x == "nl_non_sparse_identity" ~ "Non-Sparse NLM, Identity Cov.",
    x == "lm_sparse_block" ~ "Sparse LM, Block Cov.",
    x == "lm_non_sparse_block" ~ "Non-Sparse LM, Block Cov.",
    x == "nl_sparse_block" ~ "Sparse NLM, Block Cov.",
    x == "nl_non_sparse_block" ~ "Non-Sparse NLM, Block Cov.",
    .default = as.character(x)
  )
}
method_labeller <- function(x) {
  case_when(
    x == "aug_mod_cov_lasso" ~ "AMC LASSO",
    x == "aug_mod_cov_xgboost" ~ "AMC XGBoost",
    x == "mod_cov_lasso" ~ "MC LASSO",
    x == "mod_cov_xgboost" ~ "MC XGBoost",
    x == "crf" ~ "Causal Random Forest",
    x == "plugin_lasso" ~ "Plug-In LASSO",
    x == "plugin_xgboost" ~ "Plug-In XGBoost",
    x == "np_aipw_lm" ~ "AIPW-Based LASSO",
    x == "np_aipw_sl" ~ "AIPW-Based Super Learner",
  )
}

method_labeller <- function(x) {
  case_when(
    x == "aug_mod_cov_lasso" ~ "Augmented Modified Covariates LASSO",
    x == "aug_mod_cov_xgboost" ~ "Augmented Modified Covariates XGBoost",
    x == "mod_cov_lasso" ~ "Modified Covariates LASSO",
    x == "mod_cov_xgboost" ~ "Modified Covariates XGBoost",
    x == "crf" ~ "Causal Random Forests",
    x == "plugin_lasso" ~ "Plug-In LASSO",
    x == "plugin_xgboost" ~ "Plug-In XGBoost",
    x == "np_aipw_lm" ~ "AIPW-Based LASSO",
    x == "np_aipw_sl" ~ "AIPW-Based Super Learner",
  )
}

method_order <- c(
  "Plug-In LASSO", "Plug-In XGBoost", "Modified Covariates LASSO",
  "Modified Covariates XGBoost", "Augmented Modified Covariates LASSO",
  "Augmented Modified Covariates XGBoost", "AIPW-Based LASSO",
  "AIPW-Based Super Learner", "Causal Random Forests"
)


## create the individual metric tibbles
source("R/eval-functions/interpretability-funs.R")
fdr_tbl <- lapply(
  fit_tbl_ls,
  function(fit_tbl) {
    if (str_detect(fit_tbl[1, ".dgp_name"], "non_sparse")) {
      fit_tbl %>% fdr_fun(gamma = c(rep(0.5, 50), rep(0, 450)))
    } else {
      fit_tbl %>% fdr_fun(gamma = c(rep(2, 10), rep(0, 490)))
    }
  }) %>%
  bind_rows() %>%
  mutate(
    type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Obs. Study"),
    .dgp_name = str_remove(.dgp_name, "rct_"),
    .dgp_name = str_remove(.dgp_name, "obs_"),
    .dgp_name = str_remove(.dgp_name, "_cov"),
    .dgp_name = facet_labeller(.dgp_name),
    .method_name = method_labeller(.method_name),
    .method_name = factor(.method_name, level = method_order)
  )
tpr_tbl <- lapply(
  fit_tbl_ls,
  function(fit_tbl) {
    if (str_detect(fit_tbl[1, ".dgp_name"], "non_sparse")) {
      fit_tbl %>% tpr_fun(gamma = c(rep(0.5, 50), rep(0, 450)))
    } else {
      fit_tbl %>% tpr_fun(gamma = c(rep(2, 10), rep(0, 490)))
    }
  }) %>%
  bind_rows() %>%
  mutate(
    type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Obs. Study"),
    .dgp_name = str_remove(.dgp_name, "rct_"),
    .dgp_name = str_remove(.dgp_name, "obs_"),
    .dgp_name = str_remove(.dgp_name, "_cov"),
    .dgp_name = facet_labeller(.dgp_name),
    .method_name = method_labeller(.method_name),
    .method_name = factor(.method_name, level = method_order)
  )
tnr_tbl <- lapply(
  fit_tbl_ls,
  function(fit_tbl) {
    if (str_detect(fit_tbl[1, ".dgp_name"], "non_sparse")) {
      fit_tbl %>% tnr_fun(gamma = c(rep(0.5, 50), rep(0, 450)))
    } else {
      fit_tbl %>% tnr_fun(gamma = c(rep(2, 10), rep(0, 490)))
    }
  }) %>%
  bind_rows() %>%
  mutate(
    type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Obs. Study"),
    .dgp_name = str_remove(.dgp_name, "rct_"),
    .dgp_name = str_remove(.dgp_name, "obs_"),
    .dgp_name = str_remove(.dgp_name, "_cov"),
    .dgp_name = facet_labeller(.dgp_name),
    .method_name = method_labeller(.method_name),
    .method_name = factor(.method_name, level = method_order)
  )
mean_outcome_tbl <- bind_rows(
  lapply(eval_tbl_ls, function(res_ls) res_ls$mean_outcome_itr)
) %>%
  mutate(
    type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Obs. Study"),
    .dgp_name = str_remove(.dgp_name, "rct_"),
    .dgp_name = str_remove(.dgp_name, "obs_"),
    .dgp_name = str_remove(.dgp_name, "_cov"),
    .method_name = method_labeller(.method_name),
    .method_name = factor(.method_name, level = method_order)
  ) %>%
  left_join(pop_mean_itr_df, by = c(".dgp_name" = "dgp")) %>%
  mutate(.dgp_name = facet_labeller(.dgp_name)) %>% 
  mutate(rel_accuracy = mean_outcome_itr / pop_mean_itr)
mean_fit_time_tbl <- bind_rows(
  lapply(eval_tbl_ls, function(res_ls) res_ls$mean_fit_time)
) %>%
  mutate(
    type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Obs. Study"),
    .dgp_name = str_remove(.dgp_name, "rct_"),
    .dgp_name = str_remove(.dgp_name, "obs_"),
    .dgp_name = str_remove(.dgp_name, "_cov"),
    .dgp_name = facet_labeller(.dgp_name),
    .method_name = method_labeller(.method_name),
    .method_name = factor(.method_name, level = method_order)
  )

## fdr plot
fdr_tbl %>%
  ggplot(aes(y = fdr, x = n, shape = unihtee_filtering, color = .method_name)) +
  facet_grid(
    rows = vars(.dgp_name), cols = vars(type),
    labeller = label_wrap_gen(width = 16, multi_line = TRUE)
  ) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.3) +
  xlab("Sample Size") +
  ylab("Empirical FDR (100 Replicates)") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()
ggsave(
  filename = "fdr.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 10,
  width = 10,
  scale = 1
)

## tpr plot
tpr_tbl %>%
  ggplot(aes(y = tpr, x = n, shape = unihtee_filtering, color = .method_name)) +
  facet_grid(
    rows = vars(.dgp_name), cols = vars(type),
    labeller = label_wrap_gen(width = 16, multi_line = TRUE)
  ) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.3) +
  xlab("Sample Size") +
  ylab("Empirical TPR (100 Replicates)") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()
ggsave(
  filename = "tpr.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 10,
  width = 10,
  scale = 1
)

## tnr plot
tnr_tbl %>%
  ggplot(aes(y = tnr, x = n, shape = unihtee_filtering, color = .method_name)) +
  facet_grid(
    rows = vars(.dgp_name), cols = vars(type),
    labeller = label_wrap_gen(width = 16, multi_line = TRUE)
  ) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.3) +
  xlab("Sample Size") +
  ylab("Empirical TNR (100 Replicates)") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()
ggsave(
  filename = "tnr.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 10,
  width = 10,
  scale = 1
)

## mean fit time plot
mean_fit_time_tbl %>%
  ggplot(aes(
    y = mean_fit_time, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  facet_grid(
    rows = vars(.dgp_name), cols = vars(type), scales = "free_y",
    labeller = label_wrap_gen(width = 16, multi_line = TRUE)
  ) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  stat_summary(
    aes(y = mean_fit_time, x = 250, yintercept = ..y..), inherit.aes = FALSE,
    fun = "min", geom = "hline", linetype = 2, alpha = 0.3
  ) +
  xlab("Sample Size") +
  ylab("Mean Fit Time (Seconds, 100 Replicates)") +
  scale_y_log10() +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()
ggsave(
  filename = "mean-fit-time.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 10,
  width = 10,
  scale = 1
)

## mean outcome plots
mean_outcome_tbl %>%
  ggplot(aes(
    y = rel_accuracy, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.3) +
  facet_grid(
    rows = vars(.dgp_name), cols = vars(type), scales = "free_y",
    labeller = label_wrap_gen(width = 16, multi_line = TRUE)
  ) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  xlab("Sample Size") +
  ylab("Relative Accuracy (100 Replicates)") +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()
ggsave(
  filename = "expected-outcome.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 10,
  width = 10,
  scale = 1
)

## summary plot

### create the summary table
selected_estimators <- c(
  "Plug-In LASSO", "Augmented Modified Covariates XGBoost",
  "AIPW-Based Super Learner", "Causal Random Forests"
)
selected_dgp <- "Sparse NLM, Block Cov."
selected_type <- "Obs. Study"
point_size <- 2

summary_rel_accuracy_plot <- mean_outcome_tbl %>% 
  filter(
    .dgp_name == selected_dgp,
    type == selected_type,
    .method_name %in% selected_estimators
  ) %>% 
  ggplot(aes(
    y = rel_accuracy, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  geom_hline(yintercept = 1.00, linetype = 2, alpha = 0.3) +
  geom_point(alpha = 0.5, size = point_size) +
  geom_line(alpha = 0.5) +
  xlab("Sample Size") +
  ylab("Relative Accuracy") +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()

summary_fdr_plot <- fdr_tbl %>% 
  filter(
    .dgp_name == selected_dgp,
    type == selected_type,
    .method_name %in% selected_estimators
  ) %>% 
  ggplot(aes(
    y = fdr, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.3) +
  geom_point(alpha = 0.5, size = point_size) +
  geom_line(alpha = 0.5) +
  xlab("Sample Size") +
  ylab("Empirical FDR") +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()

summary_tpr_plot <- tpr_tbl %>% 
  filter(
    .dgp_name == selected_dgp,
    type == selected_type,
    .method_name %in% selected_estimators
  ) %>% 
  ggplot(aes(
    y = tpr, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  geom_hline(yintercept = 1.00, linetype = 2, alpha = 0.3) +
  geom_point(alpha = 0.5, size = point_size) +
  geom_line(alpha = 0.5) +
  xlab("Sample Size") +
  ylab("Empirical TPR") +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_discrete(name = "CATE Estimator") +
  theme_bw()

summary_comp_time_plot <- mean_fit_time_tbl %>% 
  filter(
    .dgp_name == selected_dgp,
    type == selected_type,
    .method_name %in% selected_estimators
  ) %>% 
  ggplot(aes(
    y = mean_fit_time, x = n, shape = unihtee_filtering, color = .method_name
  )) +
  geom_point(alpha = 0.5, size = point_size) +
  geom_line(alpha = 0.5) +
  xlab("Sample Size") +
  ylab("Mean Fit Time (Seconds)") +
  scale_shape_discrete(
    name = "TEM-VIP-Based Filtering", labels = unihtee_filtering_labeller
  ) +
  scale_color_discrete(name = "CATE Estimator") +
  scale_y_log10() +
  theme_bw()

  
summary_plot <- ggarrange(
  summary_rel_accuracy_plot, summary_comp_time_plot,
  summary_fdr_plot, summary_tpr_plot,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = "AUTO"
)

summary_plot_annotated <- annotate_figure(
  summary_plot,
  fig.lab = "Observational Study with Sparse Non-Linear Outcome Model and Block Covariance Matrix",
  top = " "
)

ggsave(
  summary_plot_annotated,
  filename = "summary-plot.jpeg",
  path = "results/plots",
  dpi = "retina",
  height = 8,
  width = 15,
  scale = 0.7
)
