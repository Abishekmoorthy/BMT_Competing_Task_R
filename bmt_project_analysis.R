# =============================================================================
# PROJECT: When One Event Prevents Another:
#          Competing Risks in Bone Marrow Transplantation
#
# Purpose:
#   Reproducible analysis of the KMsurv::bmt dataset, competing-risks models,
#   diagnostics, and a simulation study comparing cumulative-incidence
#   estimation with the naive 1 - Kaplan-Meier approach.
#
# IMPORTANT:
#   - This script is the analytical pipeline, not the Shiny app.
#   - Run from RStudio or with: source("bmt_project_analysis.R")
#   - Results are written to ./bmt_project_output/
#   - The BMT data are historical and small; results are educational and
#     should not be interpreted as contemporary clinical recommendations.
# =============================================================================

# ------------------------------- 0. Setup ------------------------------------

required_packages <- c(
  "KMsurv", "survival", "cmprsk", "ggplot2", "dplyr", "tidyr",
  "purrr", "readr", "tibble", "stringr", "broom", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(KMsurv)
  library(survival)
  library(cmprsk)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(broom)
  library(scales)
})

set.seed(20261014)


# Set the project directory explicitly
project_dir <- "C:/Users/Abishek/OneDrive/Documents"

# Save all analysis results inside the project folder
output_dir <- file.path(project_dir, "bmt_project_output")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "models"), showWarnings = FALSE)
dir.create(file.path(output_dir, "simulation"), showWarnings = FALSE)


write_csv_safe <- function(x, path) {
  readr::write_csv(as.data.frame(x), path, na = "")
}

save_plot <- function(plot, filename, width = 9, height = 6) {
  ggsave(
    filename = file.path(output_dir, "figures", filename),
    plot = plot, width = width, height = height, dpi = 300,
    bg = "white"
  )
}

analysis_data <- analysis_data |>
  mutate(
    relapse_event = as.integer(event_status == 1),
    death_event = as.integer(event_status == 2),
    age10 = patient_age / 10
  )

cox_relapse <- survival::coxph(Surv(time, relapse_event) ~ disease_group + age10 + patient_sex, data = analysis_data, x = TRUE, model = TRUE)
cox_death <- survival::coxph(Surv(time, death_event) ~ disease_group + age10 + patient_sex, data = analysis_data, x = TRUE, model = TRUE)

cox_relapse_tidy <- broom::tidy(cox_relapse, exponentiate = TRUE, conf.int = TRUE)
cox_death_tidy <- broom::tidy(cox_death, exponentiate = TRUE, conf.int = TRUE)

write_csv_safe(cox_relapse_tidy, file.path(output_dir, "models", "cox_cause_specific_relapse.csv"))
write_csv_safe(cox_death_tidy, file.path(output_dir, "models", "cox_cause_specific_death.csv"))

ph_relapse <- survival::cox.zph(cox_relapse)
ph_death <- survival::cox.zph(cox_death)

write_csv_safe(as.data.frame(ph_relapse$table) |> rownames_to_column("term"), file.path(output_dir, "models", "cox_zph_relapse.csv"))
write_csv_safe(as.data.frame(ph_death$table) |> rownames_to_column("term"), file.path(output_dir, "models", "cox_zph_death.csv"))

png(file.path(output_dir, "figures", "05_cox_zph_relapse.png"), width = 1200, height = 900, res = 140)
plot(ph_relapse)
dev.off()

png(file.path(output_dir, "figures", "06_cox_zph_death.png"), width = 1200, height = 900, res = 140)
plot(ph_death)
dev.off()

# --- Section 7: Fine-Gray models ---
fg_covariates <- model.matrix(~ disease_group + age10 + patient_sex, data = analysis_data)[, -1, drop = FALSE]

fit_fine_gray <- function(failcode, covariates, data) {
  tryCatch(cmprsk::crr(ftime = data$time, fstatus = data$event_status, cov1 = covariates, failcode = failcode, cencode = 0), error = function(e) e)
}

fg_relapse <- fit_fine_gray(1, fg_covariates, analysis_data)
fg_death <- fit_fine_gray(2, fg_covariates, analysis_data)

tidy_crr <- function(fit, term_names) {
  if (inherits(fit, "error")) {
    return(tibble(term = term_names, subdistribution_HR = NA_real_, model_status = conditionMessage(fit)))
  }
  beta <- fit$coef
  se <- sqrt(diag(fit$var))
  z <- beta / se
  tibble(
    term = names(beta),
    estimate = unname(beta),
    subdistribution_HR = exp(unname(beta)),
    conf_low = exp(unname(beta) - 1.96 * unname(se)),
    conf_high = exp(unname(beta) + 1.96 * unname(se)),
    p_value = unname(2 * pnorm(abs(z), lower.tail = FALSE)),
    model_status = "Fitted"
  )
}

fg_relapse_table <- tidy_crr(fg_relapse, colnames(fg_covariates))
fg_death_table <- tidy_crr(fg_death, colnames(fg_covariates))

write_csv_safe(fg_relapse_table, file.path(output_dir, "models", "fine_gray_relapse.csv"))
write_csv_safe(fg_death_table, file.path(output_dir, "models", "fine_gray_death.csv"))

# --- Section 8: Model comparison ---
model_comparison <- full_join(
  cox_relapse_tidy |> transmute(term, cause_specific_HR = estimate, cs_low = conf.low, cs_high = conf.high, cs_p = p.value),
  fg_relapse_table |> transmute(term, subdistribution_HR, fg_low = conf_low, fg_high = conf_high, fg_p = p_value),
  by = "term"
)
write_csv_safe(model_comparison, file.path(output_dir, "models", "model_comparison_relapse.csv"))

# --- Section 9 & 10: Simulation Study ---
aj_cif_at <- function(time, status, tau, cause = 1L) {
  ok <- is.finite(time) & !is.na(status) & time >= 0
  time <- time[ok]
  status <- status[ok]
  event_times <- sort(unique(time[status != 0 & time <= tau]))
  if (length(event_times) == 0) return(0)
  surv_prev <- 1
  cif <- 0
  for (tt in event_times) {
    at_risk <- sum(time >= tt)
    if (at_risk <= 0) next
    d_all <- sum(time == tt & status != 0)
    d_cause <- sum(time == tt & status == cause)
    cif <- cif + surv_prev * d_cause / at_risk
    surv_prev <- surv_prev * (1 - d_all / at_risk)
  }
  cif
}

naive_km_risk_at <- function(time, status, tau, cause = 1L) {
  fit <- survival::survfit(Surv(time, status == cause) ~ 1)
  ss <- summary(fit, times = tau, extend = TRUE)
  if (length(ss$surv) == 0 || is.na(ss$surv[1])) return(NA_real_)
  1 - ss$surv[1]
}

simulate_one <- function(n, lambda1, lambda2, lambda_c, tau) {
  t_interest <- rexp(n, rate = lambda1)
  t_competing <- rexp(n, rate = lambda2)
  t_censor <- rexp(n, rate = lambda_c)
  time <- pmin(t_interest, t_competing, t_censor)
  status <- ifelse(t_interest <= t_competing & t_interest <= t_censor, 1L, ifelse(t_competing < t_interest & t_competing <= t_censor, 2L, 0L))
  true_cif <- lambda1 / (lambda1 + lambda2) * (1 - exp(-(lambda1 + lambda2) * tau))
  tibble(
    n = n, lambda1 = lambda1, lambda2 = lambda2, lambda_c = lambda_c, tau = tau, true_cif = true_cif,
    aj_estimate = aj_cif_at(time, status, tau, cause = 1L),
    naive_km_estimate = naive_km_risk_at(time, status, tau, cause = 1L)
  )
}

simulation_design <- tidyr::crossing(n = c(100L, 300L), lambda2 = c(0.02, 0.08, 0.16), censoring_level = c("Low", "High")) |>
  mutate(lambda1 = 0.08, lambda_c = if_else(censoring_level == "Low", 0.02, 0.12), tau = 5)

n_simulations <- 500L
simulation_results <- pmap_dfr(simulation_design, function(n, lambda2, censoring_level, lambda1, lambda_c, tau) {
  map_dfr(seq_len(n_simulations), function(rep) {
    simulate_one(n = n, lambda1 = lambda1, lambda2 = lambda2, lambda_c = lambda_c, tau = tau) |>
      mutate(replication = rep, censoring_level = censoring_level)
  })
})

write_csv_safe(simulation_results, file.path(output_dir, "simulation", "simulation_replicates.csv"))

simulation_long <- simulation_results |>
  pivot_longer(cols = c(aj_estimate, naive_km_estimate), names_to = "method", values_to = "estimate") |>
  mutate(method = recode(method, aj_estimate = "Aalen-Johansen", naive_km_estimate = "Naive Kaplan-Meier"), error = estimate - true_cif)

simulation_summary <- simulation_long |>
  group_by(n, lambda1, lambda2, lambda_c, tau, censoring_level, true_cif, method) |>
  summarise(
    replications = sum(is.finite(estimate)),
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias = mean(error, na.rm = TRUE),
    empirical_sd = sd(estimate, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv_safe(simulation_summary, file.path(output_dir, "simulation", "simulation_summary_uncertainty.csv"))

# --- Section 11 & 12: Metadata & Finish ---
project_metadata <- tibble(
  item = c("Dataset", "Sample size", "Primary event", "Competing event", "Random seed", "Simulations per design"),
  value = c("KMsurv::bmt", as.character(nrow(analysis_data)), "Relapse before death", "Death before relapse", "20261014", as.character(n_simulations))
)
write_csv_safe(project_metadata, file.path(output_dir, "tables", "project_metadata.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Complete reproducible pipeline executed successfully across all team sections.")