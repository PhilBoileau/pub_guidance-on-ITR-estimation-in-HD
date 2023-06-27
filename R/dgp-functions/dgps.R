################################################################################
## Nonlinear Model with Treatment Effect Modification, Continuous Outcome
################################################################################

## Input:
##   n, an integer representing the number of observations to simulate
##   gamma, a p-length vector of treatment-covariate interaction coefficients
##   cov_mat, a p x p covariance matrix of the covariates
##   rct, a logical indicating whether the treatment is assigned according to an
##     RCT
##   linear_mod, a logical specifying whether the expected conditional outcome
##     is linear
## Output:
##   A list object containing p covariates, a binary treatment indicator, and the
##   potential outcomes of n independently simulated observations.
cont_outcome_dgp_fun <- function(
  n = 125,
  gamma = c(rep(2, 10), rep(0, 490)),
  cov_mat = diag(1, nrow = 500),
  rct = TRUE,
  linear_mod = TRUE
) {

  ## the number of covariates
  p <- length(gamma)

  ## increase the number of observations by 100. these are the validation set
  ## observations
  n <- n + 100

  ## generate the biomarkers
  W <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cov_mat)
  colnames(W) <- paste0("W", seq_len(p))

  ## simulate the binary treatment assignment
  if (rct) {
    prop_score <- 0.5
  } else {
    prop_score <- plogis(0.2 * (W[, 1] + W[, 2] + W[, 3] + W[, 4]))
  }
  A <- rbinom(n, 1, prop_score)

  ## simulate the outcomes under treatment and controls
  beta <- c(rep(2, 5), rep(0, p - 5))
  epsilon <- rnorm(n = n, mean = 0, sd = 1)
  W_t <- t(W)
  if (linear_mod) {
    Y_1 <- as.vector(1 + crossprod(W_t, beta) + crossprod(W_t, gamma) + epsilon)
    Y_0 <- as.vector(crossprod(W_t, beta) + epsilon)
  } else {
    Y_1 <- as.vector(
      crossprod(W_t, beta) + 2 * atan(crossprod(W_t, gamma)) + epsilon
    )
    Y_0 <- as.vector(crossprod(W_t, beta) + epsilon)
  }
  Y <- ifelse(A == 0, Y_0, Y_1)

  ## assembled into a list
  sample_ls <- list(
    "Y" = Y,
    "Y_0" = Y_0,
    "Y_1" = Y_1,
    "A" = A,
    "W" = W,
    "rct" = rct,
    "linear_mod" = linear_mod,
    "training" = c(rep(1, n - 100), rep(0, 100))
  )

  return(sample_ls)
}

## compute the true TEM-VIPs values using Monte Carlo
get_tem_vips <- function(
  n = 1e5,
  gamma = c(rep(2, 10), rep(0, 490)),
  cov_mat = diag(1, nrow = 500),
  linear_mod = TRUE
) {

  ## generate the full data
  full_data_ls <- cont_outcome_dgp_fun(n, gamma, cov_mat, TRUE, linear_mod)

  ## compute the risk-difference-scale TEM-VIP parameters
  pot_out_diff <- full_data_ls$Y_1 - full_data_ls$Y_0
  res_vec <- apply(
    full_data_ls$W, 2, function(mod) cov(pot_out_diff, mod) / var(mod)
  )

  ## return a tibble
  tibble(
    "modifier" = names(res_vec),
    "truth" = res_vec
  )
}

## block diagonal covariance matrix generator
block_diag_cov_mat <- function(p) {

  ## generate p %% 10 blocks of 10x10
  sub_mat_ls <- lapply(
    seq_len(round(p /  10)),
    function(idx) {
      cov2cor(
        clusterGeneration::genPositiveDefMat(
        10, "onion", rangeVar = c(1, 1.1)
        )$Sigma
      )
    }
  )

  return(as.matrix(Matrix::bdiag(sub_mat_ls)))
}
