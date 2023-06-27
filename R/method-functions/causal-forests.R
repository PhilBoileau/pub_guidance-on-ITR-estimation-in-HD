library(grf)

crf_fun <- function(
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

      # fit the model
      if (length(selected_covariates) > 0) {
        x_train <- as.matrix(dt_train[, ..selected_covariates])
        x_valid <- as.matrix(dt_valid[, ..selected_covariates])
      } else {
        x_train <- matrix(1, ncol = 1, nrow = nrow(dt_train))
        x_valid <- matrix(1, ncol = 1, nrow = nrow(dt_valid))
      }
      if (rct) {
        crf_fit <- causal_forest(
          X = x_train,
          Y = dt_train$Y,
          W = dt_train$A,
          W.hat = rep(0.5, nrow(dt_train)),
          tune.parameters = "all",
          seed = 510,
          compute.oob.predictions = FALSE,
          num.threads = 1
        )
      } else {
        crf_fit <- causal_forest(
          X = x_train,
          Y = dt_train$Y,
          W = dt_train$A,
          tune.parameters = "all",
          seed = 510,
          compute.oob.predictions = FALSE,
          num.threads = 1
        )
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

  ## predict mean outcome under ITR in validation set
  diff_pred <- predict(
    crf_fit,
    newdata = x_valid,
    num.threads = 1,
    estimate.variance = FALSE
  )
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
