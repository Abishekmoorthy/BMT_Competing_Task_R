# -----------------------------------------------------------------------------
# 06_fine_gray_model.R
#
# Purpose : Fine-Gray subdistribution hazard model for relapse.
#           Unlike the cause-specific Cox model, patients who die before
#           relapse stay in the risk set, so the subdistribution hazard ratio
#           (SHR) describes the effect of each covariate on the cumulative
#           incidence (overall probability) of relapse.
# Covariates : disease group, patient age (per 10 years), patient sex
# Output  : output/tables/fine_gray_relapse.csv
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")


# Covariate matrix ------------------------------------------------------------
# crr() requires a numeric matrix, so factors are converted to indicator
# (dummy) columns. The first column (intercept) is removed.

fg_covariates <- model.matrix(
  ~ disease_group + age10 + patient_sex,
  data = analysis_data
)[, -1]


# Fit the model ---------------------------------------------------------------
# failcode = 1 : relapse is the event of interest
# cencode  = 0 : censored observations
# (event_status = 2, death before relapse, is the competing event)

fg_relapse <- crr(
  ftime    = analysis_data$time,
  fstatus  = analysis_data$event_status,
  cov1     = fg_covariates,
  failcode = 1,
  cencode  = 0
)
fg_summary <- summary(fg_relapse)
print(fg_summary)

# To model death before relapse instead, set failcode = 2.


# Results table: subdistribution hazard ratios with 95% CI --------------------
# summary.crr columns  - coef: coef, exp(coef), se(coef), z, p-value
#                      - conf.int: exp(coef), exp(-coef), 2.5%, 97.5%

fg_results <- data.frame(
  term     = colnames(fg_covariates),
  variable = term_labels[colnames(fg_covariates)],
  SHR      = round(fg_summary$conf.int[, 1], 2),
  lower_95 = round(fg_summary$conf.int[, 3], 2),
  upper_95 = round(fg_summary$conf.int[, 4], 2),
  p_value  = round(fg_summary$coef[, 5], 3),
  row.names = NULL
)
print(fg_results)
write.csv(fg_results, "output/tables/fine_gray_relapse.csv", row.names = FALSE)
