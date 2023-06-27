library(personalized)
library(data.table)

mod_cov_lasso_fun <- function(
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

  ## train the model
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

      ## estimate the propensity score
      if (rct) {
        prop_func <- function(x, trt) 0.5
      } else {
        prop_func <- personalized::create.propensity.function(
          crossfit = TRUE,
          nfolds.crossfit = 5,
          cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
        )
      }

      ## estimate the conditional average treatment effect
      ## NOTE: fit.subgroup fails when x has one column
      if (length(selected_covariates) > 1) {
        x_train <- as.matrix(dt_train[, ..selected_covariates])
        x_valid <- as.matrix(dt_valid[, ..selected_covariates])
      } else {
        x_train <- as.matrix(dt_train[, ..covariates])
        x_valid <- as.matrix(dt_valid[, ..covariates])
      }
      mod_cov <- fit.subgroup(
        x = x_train,
        trt = dt_train$A,
        y = dt_train$Y,
        propensity.func = prop_func,
        loss = "sq_loss_lasso",
        nfolds = 5
      )
    }, times = 1
    )
  })

  ## prepare predicted TEMs
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    ## extract the predicted TEM
    coefs <- mod_cov$coefficients[colnames(W), ]
    predicted_tems <- names(coefs[which(coefs != 0)])
  }

  ## predict difference in potential outcomes
  itr_classification <- as.vector(predict(mod_cov, newx = x_valid) > 0)
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

mod_cov_xgboost_fun <- function(
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

  ## train the model
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

      ## estimate the propensity score
      if (rct) {
        prop_func <- function(x, trt) 0.5
      } else {
        prop_func <- personalized::create.propensity.function(
          crossfit = TRUE,
          nfolds.crossfit = 5,
          cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
        )
      }

      ## estimate the conditional average treatment effect
      if (length(selected_covariates) > 1) {
        x_train <- as.matrix(dt_train[, ..selected_covariates])
        x_valid <- as.matrix(dt_valid[, ..selected_covariates])
      } else {
        x_train <- as.matrix(dt_train[, ..covariates])
        x_valid <- as.matrix(dt_valid[, ..covariates])
      }
      mod_cov <- fit.subgroup(
        x = x_train,
        trt = dt_train$A,
        y = dt_train$Y,
        propensity.func = prop_func,
        loss = "sq_loss_xgboost",
        nfold = 5,
        nrounds = 1000,
        early_stopping_rounds = 50,
        verbose = 0
      )
    }, times = 1
    )
  })

  ## prepare predicted TEMs
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    predicted_tems <- NA
  }

  ## predict difference in potential outcomes
  itr_classification <- as.vector(predict(mod_cov, newx = x_valid) > 0)
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


aug_mod_cov_lasso_fun <- function(
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

  ## train the model
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

      ## estimate the propensity score
      if (rct) {
        prop_func <- function(x, trt) 0.5
      } else {
        prop_func <- personalized::create.propensity.function(
          crossfit = TRUE,
          nfolds.crossfit = 5,
          cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
        )
      }

      ## specify the augmentation function
      aug_func <- create.augmentation.function(
        family = "gaussian",
        crossfit = TRUE,
        nfolds.crossfit = 5,
        cv.glmnet.args = list(type.measure = "mse", nfolds = 5)
      )

      ## estimate the conditional average treatment effect
      if (length(selected_covariates) > 1) {
        x_train <- as.matrix(dt_train[, ..selected_covariates])
        x_valid <- as.matrix(dt_valid[, ..selected_covariates])
      } else {
        x_train <- as.matrix(dt_train[, ..covariates])
        x_valid <- as.matrix(dt_valid[, ..covariates])
      }
      aug_mod_cov <- fit.subgroup(
        x = x_train,
        trt = dt_train$A,
        y = dt_train$Y,
        propensity.func = prop_func,
        augment.func = aug_func,
        loss = "sq_loss_lasso",
        nfolds = 5
      )
    }, times = 1
    )
  })

  ## prepare predicted TEMs
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    ## extract the predicted TEM
    coefs <- aug_mod_cov$coefficients[colnames(W), ]
    predicted_tems <- names(coefs[which(coefs != 0)])
  }

  ## predict difference in potential outcomes
  itr_classification <- as.vector(predict(aug_mod_cov, newx = x_valid) > 0)
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

aug_mod_cov_xgboost_fun <- function(
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

  ## train the model
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

      ## estimate the propensity score
      if (rct) {
        prop_func <- function(x, trt) 0.5
      } else {
        prop_func <- personalized::create.propensity.function(
          crossfit = TRUE,
          nfolds.crossfit = 5,
          cv.glmnet.args = list(type.measure = "auc", nfolds = 10)
        )
      }

      ## specify the augmentation function
      aug_func <- create.augmentation.function(
        family = "gaussian",
        crossfit = TRUE,
        nfolds.crossfit = 5,
        cv.glmnet.args = list(type.measure = "mse", nfolds = 5)
      )

      ## estimate the conditional average treatment effect
      if (length(selected_covariates) > 1) {
        x_train <- as.matrix(dt_train[, ..selected_covariates])
        x_valid <- as.matrix(dt_valid[, ..selected_covariates])
      } else {
        x_train <- as.matrix(dt_train[, ..covariates])
        x_valid <- as.matrix(dt_valid[, ..covariates])
      }
      aug_mod_cov <- fit.subgroup(
        x = x_train,
        trt = dt_train$A,
        y = dt_train$Y,
        propensity.func = prop_func,
        augment.func = aug_func,
        loss = "sq_loss_xgboost",
        nfold = 5,
        nrounds = 1000,
        early_stopping_rounds = 50,
        verbose = 0
      )
    }, times = 1
    )
  })

  ## prepare predicted TEMs
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    predicted_tems <- NA
  }

  ## predict difference in potential outcomes
  itr_classification <- as.vector(predict(aug_mod_cov, newx = x_valid) > 0)
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
