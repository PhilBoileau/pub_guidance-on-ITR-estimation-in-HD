fdp_fun <- function(predicted_tems, non_tems) {

  num_false_discoveries <- sum(predicted_tems %in% non_tems)
  num_discoveries <- length(predicted_tems)
  if (num_discoveries == 0) {
    fdp <- 0
  } else if (all(is.na(predicted_tems))) {
    fdp <- NA
  } else {
    fdp <- num_false_discoveries / num_discoveries
  }

  return(fdp)

}

fdr_fun <- function(fit_results, gamma) {

  ## get non-tems
  names(gamma) <- paste0("W", seq_len(length(gamma)))
  non_tems <- names(gamma)[which(gamma == 0)]

  ## compute fdr
  group_vars <- c(
    ".dgp_name", ".method_name", "n", "unihtee_filtering"
  )
  eval_out <- fit_results %>%
    mutate(fdp = map_dbl(predicted_tems, fdp_fun, non_tems)) %>%
    group_by(across({{group_vars}})) %>%
    summarize(fdr = mean(fdp), .groups = "drop")

  return(eval_out)
}

tpp_fun <- function(predicted_tems, true_tems) {

  num_discoveries <- length(predicted_tems)
  if (num_discoveries == 0) {
    tpp <- 0
  } else if (all(is.na(predicted_tems))) {
    tpp <- NA
  } else {
    num_true_preds <- sum(predicted_tems %in% true_tems)
    num_true_tems <- length(true_tems)
    tpp <- num_true_preds / num_true_tems
  }

  return(tpp)

}

tpr_fun <- function(fit_results, gamma) {

  ## get true tems
  names(gamma) <- paste0("W", seq_len(length(gamma)))
  true_tems <- names(gamma)[which(gamma != 0)]

  ## compute tpr
  group_vars <- c(
    ".dgp_name", ".method_name", "n", "unihtee_filtering"
  )
  eval_out <- fit_results %>%
    mutate(tpp = map_dbl(predicted_tems, tpp_fun, true_tems)) %>%
    group_by(across({{group_vars}})) %>%
    summarize(tpr = mean(tpp), .groups = "drop")

  return(eval_out)
}

tnp_fun <- function(predicted_tems, non_tems) {

  num_discoveries <- length(predicted_tems)
  if (num_discoveries == 0) {
    tnp <- 1
  } else if (all(is.na(predicted_tems))) {
    tnp <- NA
  } else {
    num_false_preds <- sum(predicted_tems %in% non_tems)
    num_non_tems <- length(non_tems)
    tnp <- 1 - num_false_preds / num_non_tems
  }

  return(tnp)

}

tnr_fun <- function(fit_results, gamma) {

  ## get non tems
  names(gamma) <- paste0("W", seq_len(length(gamma)))
  non_tems <- names(gamma)[which(gamma == 0)]

  ## compute tpr
  group_vars <- c(
    ".dgp_name", ".method_name", "n", "unihtee_filtering"
  )
  eval_out <- fit_results %>%
    mutate(tnp = map_dbl(predicted_tems, tnp_fun, non_tems)) %>%
    group_by(across({{group_vars}})) %>%
    summarize(tnr = mean(tnp), .groups = "drop")

  return(eval_out)
}
