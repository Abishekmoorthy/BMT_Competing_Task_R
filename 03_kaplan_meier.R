# -----------------------------------------------------------------------------
# 03_kaplan_meier.R
#
# Purpose : Kaplan-Meier estimate of overall survival by disease group, with a
#           log-rank test for differences between groups.
# Note    : Kaplan-Meier is suitable for overall survival because death is the
#           only event. It is not used for relapse on its own, because death
#           is a competing event for relapse (see 04_cumulative_incidence.R).
# Output  : output/figures/km_overall_survival.png,
#           output/tables/km_survival_1_year.csv
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")


# Kaplan-Meier estimate -------------------------------------------------------

km_fit <- survfit(Surv(os_days, died) ~ disease_group, data = analysis_data)
print(km_fit)   # deaths and median survival time in each group


# Log-rank test ---------------------------------------------------------------
# H0: overall survival is the same in the three disease groups

logrank_test <- survdiff(Surv(os_days, died) ~ disease_group, data = analysis_data)
print(logrank_test)


# Survival curves -------------------------------------------------------------

km_plot <- ggsurvplot(
  km_fit,
  data         = analysis_data,
  pval         = TRUE,        # log-rank p-value
  conf.int     = FALSE,
  risk.table   = TRUE,        # number still at risk below the plot
  censor       = TRUE,        # + marks show censored patients
  legend.title = "Disease group",
  legend.labs  = levels(analysis_data$disease_group),
  title        = "Overall survival by disease group",
  xlab         = "Days since transplant",
  ylab         = "Survival probability",
  ggtheme      = theme_minimal()
)
print(km_plot)

png("output/figures/km_overall_survival.png", width = 1400, height = 1100, res = 150)
print(km_plot)
dev.off()


# Survival probability at one year --------------------------------------------

one_year <- summary(km_fit, times = 365)

km_one_year <- data.frame(
  disease_group = levels(analysis_data$disease_group),
  at_risk       = one_year$n.risk,
  survival      = round(one_year$surv, 3),
  lower_95      = round(one_year$lower, 3),
  upper_95      = round(one_year$upper, 3)
)
print(km_one_year)
write.csv(km_one_year, "output/tables/km_survival_1_year.csv", row.names = FALSE)
