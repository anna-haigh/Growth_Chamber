library(dplyr)
library(tidyr)

log_linear_recovery <- function(data,
                                group_vars = c("Seedling_ID", "Needle_age", "Micro", "Elev", "Genotype"),
                                time_var = "time",
                                response_var = "Fv.Fm",
                                control_target = NULL,      # data frame with join keys + control_FvFm column, OR NULL
                                control_join_vars = c("Elev", "Micro", "Needle_age"),
                                fixed_target = NULL,         # numeric constant, used instead of control_target if provided
                                time_seq = seq(0, 7, by = 0.1),
                                start_vals = list(A = 0.2, B = 0.1)) {
  
  data <- data %>% rename(.time = all_of(time_var), .response = all_of(response_var))
  
  # Fit log-linear recovery curve per group
  fit_tab <- data %>%
    group_by(across(all_of(group_vars))) %>%
    reframe({
      fit <- tryCatch({
        m <- nls(.response ~ A + B * log(.time + 1),
                 data = pick(everything()),
                 start = start_vals,
                 control = nls.control(maxiter = 100))
        data.frame(A = coef(m)["A"], B = coef(m)["B"], converged = TRUE)
      }, error = function(e) {
        data.frame(A = NA_real_, B = NA_real_, converged = FALSE)
      })
      fit
    })
  
  message(sprintf("Converged: %d / %d", sum(fit_tab$converged), nrow(fit_tab)))
  
  # Generate smooth fitted curves
  fit_curves <- fit_tab %>%
    rowwise() %>%
    mutate(time = list(time_seq), fitted = list(A + B * log(time_seq + 1))) %>%
    tidyr::unnest(c(time, fitted)) %>%
    ungroup()
  
  # Determine target and solve for days-to-recovery
  if (!is.null(fixed_target)) {
    # constant target (of 0.8) applied to every group
    target_tab <- fit_tab %>% mutate(target = fixed_target)
  } else if (!is.null(control_target)) {
    # group-specific control target, joined in
    target_tab <- fit_tab %>%
      left_join(control_target, by = control_join_vars) %>%
      rename(target = control_FvFm)
  } else {
    stop("Provide either `control_target` (a data frame) or `fixed_target` (a number).")
  }
  
  target_tab <- target_tab %>%
    mutate(
      days_to_recovery = case_when(
        !converged      ~ NA_real_,
        target <= A     ~ 0,       # already at/above target at time 0
        B <= 0          ~ Inf,     # flat/declining curve never reaches target
        TRUE            ~ exp((target - A) / B) - 1
      )
    )
  
  list(
    fit_tab = fit_tab,          # per-seedling A, B, converged
    fit_curves = fit_curves,    # smooth time/fitted values for plotting
    target_tab = target_tab     # A, B, target, days_to_recovery
  )
}
