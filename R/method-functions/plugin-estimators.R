library(glmnet)
library(xgboost)
library(microbenchmark)
library(stringr)

plugin_lasso_fun <- function(
  Y, Y_0, Y_1, A, W, rct, linear_mod, training,
  unihtee_filtering
) {

  ## organize the training and validation set
  dt <- data.table::data.table(Y, A, W, training)
  covariates <- paste0("W", seq_len(ncol(W)))
  colnames(dt) <- c("Y", "A", covariates, "training")
  dt_train <- dt[training == 1, ]
  dt_valid <- dt[training == 0, ]

  ## fit to the training data
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

      ## create the model formula
      if (length(selected_covariates) > 0) {
        formula_string <- paste(
          "Y ~ A +",
          paste(covariates, collapse = " + "),
          "+",
          paste(paste0(selected_covariates, ":A"), collapse = " + ")
        )
      } else {
        formula_string <- paste("Y ~ A +", paste(covariates, collapse = " + "))
      }

      cv_glmnet_formula <- as.formula(formula_string)

      ## define the data matrix
      x <- model.matrix(cv_glmnet_formula, data = dt_train)
      y <- dt_train$Y

      # apply the CV LASSO method
      cv_lasso <- cv.glmnet(x = x, y = y, alpha = 1)

    }, times = 1
    )
  })

  ## extract the predicted TEMs
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    non_zero_coefs <- rownames(coef(cv_lasso))[
      which(coef(cv_lasso, s = "lambda.min") > 0)
    ]
    predicted_tems <- str_remove(
      non_zero_coefs[str_detect(non_zero_coefs, "A:W")], "A:"
    )
  }

  ## predict difference in potential outcomes
  dt_valid$A <- 1
  x_valid <- model.matrix(cv_glmnet_formula, data = dt_valid)
  Y_1_pred <- predict(cv_lasso, newx = x_valid, s = "lambda.min")
  dt_valid$A <- 0
  x_valid <- model.matrix(cv_glmnet_formula, data = dt_valid)
  Y_0_pred <- predict(cv_lasso, newx = x_valid, s = "lambda.min")
  itr_classification <- as.vector((Y_1_pred - Y_0_pred) > 0)
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

plugin_xgboost_fun <- function(
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

  ## fit the xgboost regression
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

      ## prepare the data for fitting
      y_train <- dt_train$Y
      dt_train <- dt_train[, c("A", ..selected_covariates)]

      # fit the model
      xgboost_fit <- xgboost(
        data = as.matrix(dt_train), label = y_train, nrounds = 1000,
        nthread = 1, verbose = 0, early_stopping_rounds = 50
      )

    }, times = 1
    )
  })

  ## extract the TEMs if possible
  if (unihtee_filtering) {
    predicted_tems <- selected_covariates
  } else {
    predicted_tems <- NA
  }

  ## predict difference in potential outcomes
  dt_valid <- dt_valid[, c("A", ..selected_covariates)]
  dt_valid$A <- 1
  Y_1_pred <- predict(xgboost_fit, newdata = as.matrix(dt_valid))
  dt_valid$A <- 0
  Y_0_pred <- predict(xgboost_fit, newdata = as.matrix(dt_valid))
  itr_classification <- as.vector((Y_1_pred - Y_0_pred) > 0)
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
