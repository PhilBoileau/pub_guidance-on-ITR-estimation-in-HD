## load required packages
library(simChef)
library(sl3)
library(glmnet)
library(xgboost)
library(ranger)
library(unihtee)
library(personalized)
library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(future)
library(microbenchmark)
library(grf)

## set the future plan
plan(multisession, workers = 28L)

## source the R scripts
sim_functions_files = list.files(
  c(
    "R/dgp-functions", "R/method-functions", "R/eval-functions",
    "R/viz-functions/"
  ),
  pattern = "*.R$", full.names = TRUE, ignore.case = TRUE
)
sapply(sim_functions_files, source)

## generate the dgp
gamma <- c(rep(0.5, 50), rep(0, 450))
cont_out_dgp <- create_dgp(
  cont_outcome_dgp_fun,
  rct = TRUE,
  gamma = gamma,
  cov_mat = diag(1, 500),
  linear_mod = TRUE
)

## generate methods
plugin_lasso_meth <- create_method(plugin_lasso_fun)
plugin_xgboost_meth <- create_method(plugin_xgboost_fun)
mod_cov_lasso_meth <- create_method(mod_cov_lasso_fun)
mod_cov_xgboost_meth <- create_method(mod_cov_xgboost_fun)
aug_mod_cov_lasso_meth <- create_method(aug_mod_cov_lasso_fun)
aug_mod_cov_xgboost_meth <- create_method(aug_mod_cov_xgboost_fun)
crf_meth <- create_method(crf_fun)
np_aipw_lm_meth <- create_method(np_aipw_lm_fun)
np_aipw_sl_meth <- create_method(np_aipw_sl_fun)

## generate the evaluators
fdr_eval <- create_evaluator(fdr_fun, gamma = gamma)
tpr_eval <- create_evaluator(tpr_fun, gamma = gamma)
tnr_eval <- create_evaluator(tnr_fun, gamma = gamma)
mean_fit_time_eval <- create_evaluator(mean_fit_time_fun)
mean_outcome_itr_eval <- create_evaluator(mean_outcome_itr_fun)

## define the experiment
experiment <- create_experiment(name = "lm-non-sparse-rct-identity-cov") %>%
  add_dgp(cont_out_dgp, name = "lm_non_sparse_rct_identity_cov") %>%
  add_vary_across(
    .dgp = "lm_non_sparse_rct_identity_cov",
    n = c(250, 500, 1000)
  ) %>%
  add_method(plugin_lasso_meth, name = "plugin_lasso") %>%
  add_method(plugin_xgboost_meth, name = "plugin_xgboost") %>%
  add_method(mod_cov_lasso_meth, name = "mod_cov_lasso") %>%
  add_method(mod_cov_xgboost_meth, name = "mod_cov_xgboost") %>%
  add_method(aug_mod_cov_lasso_meth, name = "aug_mod_cov_lasso") %>%
  add_method(aug_mod_cov_xgboost_meth, name = "aug_mod_cov_xgboost") %>%
  add_method(crf_meth, name = "crf") %>%
  add_method(np_aipw_lm_meth, name = "np_aipw_lm") %>%
  add_method(np_aipw_sl_meth, name = "np_aipw_sl") %>%
  add_vary_across(
    .method = c(
      "plugin_lasso", "plugin_xgboost", "mod_cov_lasso", "mod_cov_xgboost",
      "aug_mod_cov_lasso", "aug_mod_cov_xgboost", "crf", "np_aipw_lm",
      "np_aipw_sl"
    ),
    unihtee_filtering = c(TRUE, FALSE)
  ) %>%
  add_evaluator(fdr_eval, name = "fdr") %>%
  add_evaluator(tpr_eval, name = "tpr") %>%
  add_evaluator(tnr_eval, name = "tnr") %>%
  add_evaluator(mean_fit_time_eval, name = "mean_fit_time") %>%
  add_evaluator(mean_outcome_itr_eval, name = "mean_outcome_itr")

## run the simulation study
set.seed(514)
results <- experiment$run(
  n_reps = 100,
  future.globals = c(
    "gamma", "unihtee_rct_fun", "unihtee_obs_study_fun",
    "compute_aipw_transform"
  ),
  future.packages = c(
    "sl3", "glmnet", "xgboost", "ranger", "personalized", "unihtee", "dplyr",
    "stringr", "purrr", "data.table", "microbenchmark", "grf"
  ),
  save = TRUE,
  checkpoint_n_reps = 20
)
