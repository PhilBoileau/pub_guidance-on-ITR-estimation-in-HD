################################################################################
## unihtee wrappers
################################################################################

library(unihtee)
library(glmnet)
library(xgboost)
library(ranger)
library(data.table)
library(sl3)

## For observational study, continuous outcome simulations
unihtee_obs_study_fun <- function(dt) {

  ## name the covariates
  covariates <- colnames(dt)[str_detect(colnames(dt), "W")]

  ## define the super learner to use## define base learners
  ## define the interactions
  interactions <- lapply(covariates, function(w) c(w, "A"))
  lrnr_interactions <- Lrnr_define_interactions$new(interactions)
  lrnr_lasso_co <- make_learner(
    Pipeline, lrnr_interactions, Lrnr_glmnet$new()
  )
  lrnr_enet_co <- make_learner(
    Pipeline, lrnr_interactions, Lrnr_glmnet$new(alpha = 0.5)
  )
  lrnr_ridge_co <- make_learner(
    Pipeline, lrnr_interactions, Lrnr_glmnet$new(alpha = 0)
  )
  lrnr_lasso_ps <- Lrnr_glmnet$new()
  lrnr_enet_ps <- Lrnr_glmnet$new(alpha = 0.5)
  lrnr_ridge_ps <- Lrnr_glmnet$new(alpha = 0)
  lrnr_xgboost <- Lrnr_xgboost$new(
    nrounds = 1000, early_stopping_rounds = 50
  )
  lrnr_ranger <- Lrnr_ranger$new()

  ## define the super learners
  lrnr_sl_ps <- Lrnr_sl$new(
    learners = list(
      lrnr_lasso_ps, lrnr_enet_ps, lrnr_ridge_ps, lrnr_xgboost, lrnr_ranger
    ),
    metalearner = make_learner(
      Lrnr_solnp, metalearner_logistic_binomial,
      loss_squared_error
    )
  )
  lrnr_sl_co <- Lrnr_sl$new(
    learners = list(
      lrnr_lasso_co, lrnr_enet_co, lrnr_ridge_co, lrnr_xgboost, lrnr_ranger
    ),
    metalearner = make_learner(
      Lrnr_solnp, metalearner_linear, loss_squared_error
    )
  )

  ## apply unicate for RD TEM VIP to continuous outcome
  results <- unihtee::unihtee(
    data = dt,
    confounders = covariates,
    modifiers = covariates,
    exposure = "A",
    outcome = "Y",
    outcome_type = "continuous",
    effect = "absolute",
    cross_fit = FALSE,
    estimator = "onestep",
    prop_score_estimator = lrnr_sl_ps,
    cond_outcome_estimator = lrnr_sl_co
  )

  return(results)

}

## For RCT continuous outcome simulations
unihtee_rct_fun <- function(dt) {

  ## name the covariates
  covariates <- colnames(dt)[str_detect(colnames(dt), "W")]

  ## add RCT propensity scores
  dt$prop_scores <- 0.5

  ## define LASSO regression for the expected conditional outcome
  interactions <- lapply(covariates, function(w) c(w, "A"))
  lrnr_interactions <- Lrnr_define_interactions$new(interactions)
  lrnr_lasso_co <- make_learner(
    Pipeline, lrnr_interactions, Lrnr_glmnet$new()
  )

  ## apply unicate for RD TEM VIP to continuous outcome
  results <- unihtee::unihtee(
    data = dt,
    confounders = covariates,
    modifiers = covariates,
    exposure = "A",
    outcome = "Y",
    outcome_type = "continuous",
    effect = "absolute",
    cross_fit = FALSE,
    estimator = "onestep",
    prop_score_values = "prop_scores",
    cond_outcome_estimator = lrnr_lasso_co
  )

  return(results)

}
