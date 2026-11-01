# -----------------------------------------------------------------------------
# 05_cox_model.R
#
# Purpose : Cause-specific Cox proportional hazards model for relapse.
#           Deaths before relapse are treated as censored, so the model
#           estimates the effect of each covariate on the relapse hazard
#           among patients who are still alive and relapse-free.
# Covariates : disease group, patient age (per 10 years), patient sex
# Output  : output/tables/cox_relapse.csv, output/tables/cox_ph_test.csv
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")


# Fit the model ---------------------------------------------------------------

cox_relapse <- coxph(
  Surv(time, relapse_event) ~ disease_group + age10 + patient_sex,
  data = analysis_data
)
print(summary(cox_relapse))

# To model death before relapse instead, replace relapse_event with
# death_event in the formula above.


# Results table: hazard ratios with 95% confidence intervals -----------------

cox_summary <- summary(cox_relapse)

cox_results <- data.frame(
  term     = rownames(cox_summary$conf.int),
  variable = term_labels[rownames(cox_summary$conf.int)],
  HR       = round(cox_summary$conf.int[, "exp(coef)"], 2),
  lower_95 = round(cox_summary$conf.int[, "lower .95"], 2),
  upper_95 = round(cox_summary$conf.int[, "upper .95"], 2),
  p_value  = round(cox_summary$coefficients[, "Pr(>|z|)"], 3),
  row.names = NULL
)
print(cox_results)
write.csv(cox_results, "output/tables/cox_relapse.csv", row.names = FALSE)


# Proportional hazards assumption ---------------------------------------------
# Test based on scaled Schoenfeld residuals. A p-value below 0.05 suggests
# the effect of that covariate changes over time.

ph_test <- cox.zph(cox_relapse)

ph_results <- data.frame(
  term       = rownames(ph_test$table),
  chi_square = round(ph_test$table[, "chisq"], 2),
  df         = ph_test$table[, "df"],
  p_value    = round(ph_test$table[, "p"], 3),
  row.names  = NULL
)
print(ph_results)
write.csv(ph_results, "output/tables/cox_ph_test.csv", row.names = FALSE)
