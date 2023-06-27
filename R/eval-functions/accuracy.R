## mean of mean outcome unde ITR
mean_outcome_itr_fun <- function(fit_results) {

  group_vars <- c(
    ".dgp_name", ".method_name", "n", "unihtee_filtering"
  )
  eval_out <- fit_results %>%
    group_by(across({{group_vars}})) %>%
    summarize(mean_outcome_itr = mean(mean_outcome_itr), .groups = "drop")

  return(eval_out)
}
