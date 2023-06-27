 ###############################################################################
 ## Compute Population Mean ITR Value
 ###############################################################################

## load the DGP-functions
source("R/dgp-functions/dgps.R")

## define the conditional outcome functions for each DGP
mu_1 <- function(a, w) {
  beta <- c(rep(2, 5), rep(0, 495))
  gamma <- c(rep(2, 10), rep(0, 490))
  w_t <- t(w)
  a + crossprod(w_t, beta) + a * crossprod(w_t, gamma)
}
mu_2 <- function(a, w) {
  beta <- c(rep(2, 5), rep(0, 495))
  gamma <- c(rep(0.5, 50), rep(0, 450))
  w_t <- t(w)
  a + crossprod(w_t, beta) + a * crossprod(w_t, gamma)
}
mu_3 <- function(a, w) {
  beta <- c(rep(2, 5), rep(0, 495))
  gamma <- c(rep(2, 10), rep(0, 490))
  w_t <- t(w)
  crossprod(w_t, beta) + 2 * a * atan(crossprod(w_t, gamma))
}
mu_4 <- function(a, w) {
  beta <- c(rep(2, 5), rep(0, 495))
  gamma <- c(rep(0.5, 50), rep(0, 450))
  w_t <- t(w)
  crossprod(w_t, beta) + 2 * a * atan(crossprod(w_t, gamma))
}

## compute the monte-carlo population ITR means
population_size <- 1000000

## sparse linear model with block covariance matrix
set.seed(514)
gamma <- c(rep(2, 10), rep(0, 490))
cov_mat <- block_diag_cov_mat(p = 500)
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = cov_mat,
  rct = TRUE,
  linear_mod = TRUE
)
itr <- (mu_1(1, dat$W) - mu_1(0, dat$W)) > 0
sparse_lm_block_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## sparse linear model with identity covariance matrix
set.seed(514)
gamma <- c(rep(2, 10), rep(0, 490))
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = diag(500),
  rct = TRUE,
  linear_mod = TRUE
)
itr <- (mu_1(1, dat$W) - mu_1(0, dat$W)) > 0
sparse_lm_id_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## non-sparse linear model with block covariance matrix
set.seed(514)
gamma <- c(rep(0.5, 50), rep(0, 450))
cov_mat <- block_diag_cov_mat(p = 500)
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = cov_mat,
  rct = TRUE,
  linear_mod = TRUE
)
itr <- (mu_2(1, dat$W) - mu_2(0, dat$W)) > 0
ns_lm_block_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## non-sparse linear model with identity covariance matrix
set.seed(514)
gamma <- c(rep(0.5, 50), rep(0, 450))
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = diag(500),
  rct = TRUE,
  linear_mod = TRUE
)
itr <- (mu_2(1, dat$W) - mu_2(0, dat$W)) > 0
ns_lm_id_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## sparse nonlinear model with block covariance matrix
set.seed(514)
gamma <- c(rep(2, 10), rep(0, 490))
cov_mat <- block_diag_cov_mat(p = 500)
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = cov_mat,
  rct = TRUE,
  linear_mod = FALSE
)
itr <- (mu_3(1, dat$W) - mu_3(0, dat$W)) > 0
sparse_nlm_block_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## sparse nonlinear model with identity covariance matrix
set.seed(514)
gamma <- c(rep(2, 10), rep(0, 490))
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = diag(500),
  rct = TRUE,
  linear_mod = FALSE
)
itr <- (mu_3(1, dat$W) - mu_3(0, dat$W)) > 0
sparse_nlm_id_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## non-sparse nonlinear model with block covariance matrix
set.seed(514)
gamma <- c(rep(0.5, 50), rep(0, 450))
cov_mat <- block_diag_cov_mat(p = 500)
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = cov_mat,
  rct = TRUE,
  linear_mod = FALSE
)
itr <- (mu_4(1, dat$W) - mu_4(0, dat$W)) > 0
ns_nlm_block_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## non-sparse nonlinear model with identity covariance matrix
set.seed(514)
gamma <- c(rep(0.5, 50), rep(0, 450))
dat <- cont_outcome_dgp_fun(
  n = population_size,
  gamma = gamma,
  cov_mat = diag(500),
  rct = TRUE,
  linear_mod = FALSE
)
itr <- (mu_4(1, dat$W) - mu_4(0, dat$W)) > 0
ns_nlm_id_value <- mean(c(dat$Y_1[itr == 1], dat$Y_0[itr == 0]))

## save the parameters in a data.frame
parameter_df <- data.frame(
  dgp = c("lm_sparse_identity", "lm_non_sparse_identity", "nl_sparse_identity",
          "nl_non_sparse_identity", "lm_sparse_block", "lm_non_sparse_block",
          "nl_sparse_block", "nl_non_sparse_block"),
  pop_mean_itr= c(sparse_lm_id_value, ns_lm_id_value,
                         sparse_nlm_id_value, ns_nlm_id_value,
                         sparse_lm_block_value, ns_lm_block_value,
                         sparse_nlm_block_value, ns_nlm_block_value)
)
saveRDS(parameter_df, "results/population-mean-itr.Rds")
