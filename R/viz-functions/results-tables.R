library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(xtable)

# setup ----

## read in all the fit results ----
fit_result_rds_paths <- list.files(
  path = "results", pattern = "*fit_results.rds*", full.names = TRUE,
  recursive = TRUE
)
fit_tbl_ls <- lapply(fit_result_rds_paths, readRDS)

## read in all the evaluation results ----
eval_result_rds_paths <- list.files(
  path = "results", pattern = "*eval_results.rds", full.names = TRUE,
  recursive = TRUE
)
eval_tbl_ls <- lapply(eval_result_rds_paths, readRDS)

## read in the monte-carlo-based oracle ITR mean outcomes ----
pop_mean_itr_df <- readRDS("results/population-mean-itr.Rds")

## function for classifying type of DGP
dgp_type_fun <- function(tbl) {
  tbl %>%
    mutate(
      type = ifelse(str_detect(.dgp_name, "rct"), "RCT", "Observational Study"),
      .dgp_name = str_remove(.dgp_name, "rct_"),
      .dgp_name = str_remove(.dgp_name, "obs_"),
      .dgp_name = str_remove(.dgp_name, "_cov")
    )
}

## define some labelling functions ----
dgp_labeller <- function(x) {
  case_when(
    x == "lm_sparse_identity" ~ "Sparse Linear Model, Identity Covariance",
    x == "lm_non_sparse_identity" ~ "Non-Sparse Linear Model, Identity Covariance",
    x == "nl_sparse_identity" ~ "Sparse Non-Linear Model, Identity Covariance",
    x == "nl_non_sparse_identity" ~ "Non-Sparse Non-Linear Model, Identity Covariance",
    x == "lm_sparse_block" ~ "Sparse Linear Model, Block Covariance",
    x == "lm_non_sparse_block" ~ "Non-Sparse Linear Model, Block Covariance",
    x == "nl_sparse_block" ~ "Sparse Non-Linear Model, Block Covariance",
    x == "nl_non_sparse_block" ~ "Non-Sparse Non-Linear Model, Block Covariance",
    .default = as.character(x)
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

## create the individual metric tibbles ----
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
  dgp_type_fun()

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
  dgp_type_fun()

tnr_tbl <- lapply(
  fit_tbl_ls,
  function(fit_tbl) {
    if (str_detect(fit_tbl[1, ".dgp_name"], "non_sparse")) {
      fit_tbl %>% tnr_fun(gamma = c(rep(0.5, 50), rep(0, 450)))
    } else {
      fit_tbl %>% tnr_fun(gamma = c(rep(2, 10), rep(0, 490)))
    }
  }) %>%
  bind_rows()  %>%
  dgp_type_fun()

mean_outcome_tbl <- bind_rows(
  lapply(eval_tbl_ls, function(res_ls) res_ls$mean_outcome_itr)
)  %>%
  dgp_type_fun() %>% 
  left_join(pop_mean_itr_df, by = c(".dgp_name" = "dgp")) %>%
  mutate(rel_mean_itr = mean_outcome_itr / pop_mean_itr)

mean_fit_time_tbl <- bind_rows(
  lapply(eval_tbl_ls, function(res_ls) res_ls$mean_fit_time)
) %>%
  dgp_type_fun()

# combined table ----
results_tbl <- mean_outcome_tbl %>% 
  left_join(
    mean_fit_time_tbl,
    by = c(".dgp_name", '.method_name', "n", "unihtee_filtering", "type")
  ) %>% 
  left_join(
    fdr_tbl,
    by = c(".dgp_name", '.method_name', "n", "unihtee_filtering", "type")
  ) %>% 
  left_join(
    tnr_tbl,
    by = c(".dgp_name", '.method_name', "n", "unihtee_filtering", "type")
  ) %>% 
  left_join(
    tpr_tbl,
    by = c(".dgp_name", '.method_name', "n", "unihtee_filtering", "type")
  ) %>% 
  mutate(
    dgp = dgp_labeller(.dgp_name),
    method = method_labeller(.method_name),
    method = factor(method, level = method_order),
    "Relative rule quality" = format(round(rel_mean_itr, 2), nsmall = 2),
    "Mean fit time (sec.)" = format(round(mean_fit_time, 2), nsmall = 2),
    "Empirical FDR (%)" = format(round(fdr * 100, 2), nsmall = 2),
    "Empirical TPR (%)" = format(round(tpr * 100, 2), nsmall = 2),
    "Empirical TNR (%)" = format(round(tnr * 100, 2), nsmall = 2)
  ) %>% 
  select(
    dgp, type, method, n, unihtee_filtering, `Relative rule quality`,
    `Mean fit time (sec.)`, `Empirical FDR (%)`, `Empirical TPR (%)`,
    `Empirical TNR (%)`
  ) %>%
  pivot_longer(
    cols = c(`Relative rule quality`, `Mean fit time (sec.)`, `Empirical FDR (%)`,
             `Empirical TPR (%)`, `Empirical TNR (%)`),
    names_to = "metric"
  ) %>% 
  pivot_wider(
    names_from = n,
    names_prefix = "n = ",
    values_from = c(value)
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Relative rule quality", "Empirical FDR (%)", "Empirical TPR (%)",
                 "Empirical TNR (%)", "Mean fit time (sec.)")
    )
  ) %>% 
  arrange(method, metric)

# dplyr-based approach ----

extract_latex_table_fun <- function(results_tbl, dgp_name, type_name) {
  
  no_filtering_tbl <- results_tbl %>%
    filter(
      dgp == dgp_name,
      type == type_name,
      unihtee_filtering == FALSE
    ) %>% 
    select(-type, -unihtee_filtering)
  
  filtering_tbl <- results_tbl %>%
    filter(
      dgp == dgp_name,
      type == type_name,
      unihtee_filtering == TRUE
    ) %>% 
    select(-type, -unihtee_filtering)
  
  dgp_tbl <- no_filtering_tbl %>% 
    left_join(
      filtering_tbl,
      by = c("dgp", "method", "metric"),
      suffix = c("", " ")
    ) %>%
    select(-dgp) %>% 
    mutate(method = if_else(metric == "Relative rule quality", method, "")) %>% 
    rename(`CATE Estimator` = method, Metric = metric)
  
  dgp_xtable <- xtable(
    dgp_tbl,
    caption = paste(
      c("Simulation results: ", dgp_name, ", ", type_name),
      collapse = ""),
    align = c("r", "|l||", "l|", "r", "r", "r|", "r", "r", "r|")
  )
  
  print(dgp_xtable, include.rownames = FALSE)
}

# print tables ----

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Linear Model, Identity Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Linear Model, Identity Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Linear Model, Block Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Linear Model, Block Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Non-Linear Model, Identity Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Non-Linear Model, Identity Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Non-Linear Model, Block Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Sparse Non-Linear Model, Block Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Linear Model, Identity Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Linear Model, Identity Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Linear Model, Block Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Linear Model, Block Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Non-Linear Model, Identity Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Non-Linear Model, Identity Covariance",
    "Observational Study"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Non-Linear Model, Block Covariance",
    "RCT"
  )

results_tbl %>%
  extract_latex_table_fun(
    "Non-Sparse Non-Linear Model, Block Covariance",
    "Observational Study"
  )
