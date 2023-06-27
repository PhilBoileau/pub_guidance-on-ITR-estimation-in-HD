## mean fit time in seconds
mean_fit_time_fun <- function(fit_results) {

  group_vars <- c(
    ".dgp_name", ".method_name", "n", "unihtee_filtering"
  )
  eval_out <- fit_results %>%
    group_by(across({{group_vars}})) %>%
    summarize(mean_fit_time = mean(fit_time), .groups = "drop")

  return(eval_out)
}
