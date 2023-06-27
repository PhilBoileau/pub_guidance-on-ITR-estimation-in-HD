library(sl3)
library(ranger)
library(xgboost)
library(glmnet)
library(stringr)

np_aipw_lm_fun <- function(
  Y, Y_0, Y_1, A, W, rct, linear_mod, training,
  unihtee_filtering
) {

  ## organize the training set
  dt <- data.table::data.table(Y, A, W, training)
  covariates <- paste0("W", seq_len(ncol(W)))
  colnames(dt) <- c("Y", "A", covariates, "training")
  dt_train <- dt[training == 1, ]
  dt_valid <- dt[training == 0, ]
  dt_train$training <- NULL
  dt_valid$training <- NULL

  ## fit model
  capture.output({
    mb <- microbenchmark(
    {
      ## specify the covariates to use in the model
      if (unihtee_filtering) {
        if (rct) {
          unihtee_results <- unihtee_rct_fun(dt_train)
        } else {
          unihtee_results <- unihtee_obs_study_fun(dt_train)
        }
        selected_covariates <- as.character(
          unihtee_results[p_value_fdr < 0.05, modifier]
        )
      } else {
        selected_covariates <- covariates
      }

      ## compute the predicted differences in potential outcomes
      dt_train <- compute_aipw_transform(dt_train, selected_covariates, rct)

      ## fit a LASSO regression

      if (length(selected_covariates) > 0) {
        ## create the model formula
        formula_string <- paste(
          "pred_diff ~ ", paste(selected_covariates, collapse = " + ")
        )
        cv_glmnet_formula <- as.formula(formula_string)

        ## define the data matrix
        x <- model.matrix(cv_glmnet_formula, data = dt_train)
        pred_diff <- dt_train$pred_diff

        ## apply the CV LASSO method
        cv_lasso <- cv.glmnet(x = x, y = pred_diff, alpha = 1)

      } else {
        cv_lasso <- mean(dt_train$pred_diff)
      }

    }, times = 1
    )
  })

  ## extract the predicted TEM
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    non_zero_coefs <- rownames(coef(cv_lasso))[
      which(coef(cv_lasso, s = "lambda.min") > 0)]
    predicted_tems <- non_zero_coefs[str_detect(non_zero_coefs, "W")]
  }

  ## predict difference in potential outcomes
  if (length(selected_covariates) > 0) {
    dt_valid$pred_diff <- 0
    x_valid <- model.matrix(cv_glmnet_formula, data = dt_valid)
    diff_pred <- predict(cv_lasso, newx = x_valid, s = "lambda.min")
  } else {
    diff_pred <- rep(cv_lasso, nrow(dt_valid))
  }
  itr_classification <- as.vector(diff_pred > 0)
  Y_1_valid <- Y_1[training == 0]
  Y_0_valid <- Y_0[training == 0]
  mean_outcome_itr <- 1 / length(itr_classification) * (
    sum(Y_1_valid[itr_classification]) + sum(Y_0_valid[!itr_classification])
  )

  # return the object containing the fit and possible feature selection
  return(list(
    "predicted_tems" = predicted_tems,
    "fit_time" = mb$time / 10^9,
    "mean_outcome_itr" = mean_outcome_itr
  ))

}

np_aipw_sl_fun <- function(
  Y, Y_0, Y_1, A, W, rct, linear_mod, training,
  unihtee_filtering
) {

  ## organize the training set
  dt <- data.table::data.table(Y, A, W, training)
  covariates <- paste0("W", seq_len(ncol(W)))
  colnames(dt) <- c("Y", "A", covariates, "training")
  dt_train <- dt[training == 1, ]
  dt_valid <- dt[training == 0, ]
  dt_train$training <- NULL
  dt_valid$training <- NULL

  ## fit model
  capture.output({
    mb <- microbenchmark(
    {
      ## specify the covariates to use in the model
      if (unihtee_filtering) {
        if (rct) {
          unihtee_results <- unihtee_rct_fun(dt_train)
        } else {
          unihtee_results <- unihtee_obs_study_fun(dt_train)
        }
        selected_covariates <- as.character(
          unihtee_results[p_value_fdr < 0.05, modifier]
        )
      } else {
        selected_covariates <- covariates
      }

      ## compute the predicted differences in potential outcomes
      dt_train <- compute_aipw_transform(dt_train, selected_covariates, rct)

      ## fit a blip regression

      if (length(selected_covariates) > 0) {
        ## define the sl3 task
        blip_task <- sl3_Task$new(
          data = dt_train,
          covariates = selected_covariates,
          outcome = "pred_diff",
          outcome_type = "continuous",
          folds = 5L
        )

        ## define the sl
        lrnr_enet <- Lrnr_glmnet$new(alpha = 0.5)
        lrnr_ridge <- Lrnr_glmnet$new(alpha = 0)
        lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)
        lrnr_ranger <- Lrnr_ranger$new()
        lrnr_xgboost <- Lrnr_xgboost$new(
          nrounds = 1000, early_stopping_rounds = 50, verbose = 0
        )
        lrnr_sl <- Lrnr_sl$new(
          learners = list(
            lrnr_lasso, lrnr_enet, lrnr_ridge, lrnr_xgboost, lrnr_ranger
          ),
          metalearner = make_learner(
            Lrnr_solnp, metalearner_linear, loss_squared_error
          )
        )

        ## fit the model
        blip_fit <- lrnr_sl$train(blip_task)
      } else {
        blip_fit <- mean(dt_train$pred_diff)
      }

    }, times = 1
    )
  })

  ## extract predicted TEMs if possible
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    predicted_tems <- NA
  }

  ## predict difference in potential outcomes
  if (length(selected_covariates) > 0) {
    dt_valid$pred_diff <- 0
    blip_valid_task <- sl3_Task$new(
      data = dt_valid,
      covariates = selected_covariates,
      outcome_type = "continuous",
      outcome = "pred_diff"
    )
    pred_diff <- blip_fit$predict(blip_valid_task)
  } else {
    pred_diff <- rep(blip_fit, nrow(dt_valid))
  }
  itr_classification <- as.vector(pred_diff > 0)
  Y_1_valid <- Y_1[training == 0]
  Y_0_valid <- Y_0[training == 0]
  mean_outcome_itr <- 1 / length(itr_classification) * (
    sum(Y_1_valid[itr_classification]) + sum(Y_0_valid[!itr_classification])
  )

  # return the object containing the fit and possible feature selection
  return(list(
    "predicted_tems" = predicted_tems,
    "fit_time" = mb$time / 10^9,
    "mean_outcome_itr" = mean_outcome_itr
  ))

}

compute_aipw_transform <- function(dt, selected_covariates, rct) {

  ## construct the expected conditional outcome task
  eco_task <- sl3_Task$new(
    data = dt,
    covariates = c(selected_covariates, "A"),
    outcome = "Y",
    outcome_type = "continuous",
    folds = 5L
  )

  ## define the super learner
  if (length(selected_covariates) > 0) {
    interactions <- lapply(selected_covariates, function(w) c(w, "A"))
    lrnr_interactions <- Lrnr_define_interactions$new(interactions)
    lrnr_enet_eco <- make_learner(
      Pipeline, lrnr_interactions, Lrnr_glmnet$new(alpha = 0.5)
    )
    lrnr_ridge_eco <- make_learner(
      Pipeline, lrnr_interactions, Lrnr_glmnet$new(alpha = 0)
    )
    lrnr_lasso_eco <- make_learner(
      Pipeline, lrnr_interactions, Lrnr_glmnet$new(alpha = 1)
    )
  } else {
    lrnr_enet_eco <- Lrnr_mean$new()
    lrnr_ridge_eco <- Lrnr_mean$new()
    lrnr_lasso_eco <- Lrnr_mean$new()
  }
  lrnr_xgboost <- Lrnr_xgboost$new(
    nrounds = 1000, early_stopping_rounds = 50, verbose = 0
  )
  lrnr_ranger <- Lrnr_ranger$new()
  lrnr_sl_eco <- Lrnr_sl$new(
    learners = list(
      lrnr_lasso_eco, lrnr_enet_eco, lrnr_ridge_eco, lrnr_xgboost, lrnr_ranger
    ),
    metalearner = make_learner(
      Lrnr_solnp, metalearner_linear, loss_squared_error
    )
  )

  ## train the expected outcome super learner
  eco_fit <- lrnr_sl_eco$train(eco_task)

  ## get the predicted outcome
  dt$Y_pred <- eco_fit$predict()

  ## if in an observational setting, estimate the propensity score
  if (!rct) {

    if (length(selected_covariates) > 0) {
      ## construct the expected conditional outcome task
      ps_task <- sl3_Task$new(
        data = dt,
        covariates = selected_covariates,
        outcome = "A",
        outcome_type = "binary",
        folds = 5L
      )

      ## define the super learner
      lrnr_enet_ps <- Lrnr_glmnet$new(alpha = 0.5)
      lrnr_ridge_ps <- Lrnr_glmnet$new(alpha = 0)
      lrnr_lasso_ps <- Lrnr_glmnet$new(alpha = 1)
      lrnr_xgboost <- Lrnr_xgboost$new(
        nrounds = 1000, early_stopping_rounds = 50, verbose = 0
      )
      lrnr_ranger <- Lrnr_ranger$new()
      lrnr_sl_ps <- Lrnr_sl$new(
        learners = list(
          lrnr_lasso_ps, lrnr_enet_ps, lrnr_ridge_ps, lrnr_xgboost, lrnr_ranger
        ),
        metalearner = make_learner(
          Lrnr_solnp, metalearner_logistic_binomial, loss_squared_error
        )
      )

      ## train the expected outcome super learner
      ps_fit <- lrnr_sl_ps$train(ps_task)

      ## get predictions
      ps_pred <- ps_fit$predict()

    } else {
      ps_pred <- rep(mean(dt$A), nrow(dt))
    }

  } else {
    ps_pred <- rep(0.5, nrow(dt))
  }
  dt$ps_pred <- ps_pred

  ## compute potential outcome predictions
  exp_dt <- dt
  exp_dt$A <- 1
  exp_task <- sl3_Task$new(
    data = exp_dt,
    covariates = c(selected_covariates, "A"),
    outcome = "Y",
    outcome_type = "continuous"
  )
  dt$Y_1_pred <- eco_fit$predict(exp_task)

  no_exp_dt <- dt
  no_exp_dt$A <- 0
  no_exp_task <- sl3_Task$new(
    data = no_exp_dt,
    covariates = c(selected_covariates, "A"),
    outcome = "Y",
    outcome_type = "continuous"
  )
  dt$Y_0_pred <- eco_fit$predict(no_exp_task)

  ## compute the predicted diffrence in potential outcomes
  dt$pred_diff <- ((2 * dt$A - 1) /
    (dt$A * dt$ps_pred + (1 - dt$A) * (1 - dt$ps_pred))) * (dt$Y - dt$Y_pred) +
    dt$Y_1_pred - dt$Y_0_pred

  return(dt)

}
