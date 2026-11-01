# -----------------------------------------------------------------------------
# 04_cumulative_incidence.R
#
# Purpose : Estimate the cumulative incidence of relapse and of death before
#           relapse using the Aalen-Johansen estimator (cmprsk::cuminc), and
#           compare disease groups with Gray's test.
#           Also shows why 1 - Kaplan-Meier overestimates the relapse risk
#           when deaths are treated as censored observations.
# Output  : output/figures/cif_*.png, output/tables/cif_timepoints.csv
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")


# First event observed for each patient ---------------------------------------

print(table(analysis_data$disease_group, analysis_data$event_label))


# Cumulative incidence - all patients -----------------------------------------
# cuminc() returns one curve per event type:
#   "1 1" = relapse, "1 2" = death before relapse

cif_overall <- cuminc(
  ftime   = analysis_data$time,
  fstatus = analysis_data$event_status,
  cencode = 0
)

# Estimated probabilities at 100 days, 1 year and 2 years
cif_times <- timepoints(cif_overall, times = c(100, 365, 730))$est
rownames(cif_times) <- c("Relapse", "Death before relapse")
print(round(cif_times, 3))
write.csv(data.frame(event = rownames(cif_times), round(cif_times, 3), check.names = FALSE),
          "output/tables/cif_timepoints.csv", row.names = FALSE)


# Cumulative incidence by disease group, with Gray's test ---------------------

cif_group <- cuminc(
  ftime   = analysis_data$time,
  fstatus = analysis_data$event_status,
  group   = analysis_data$disease_group,
  cencode = 0
)

# Gray's test: H0 = cumulative incidence is the same across groups
# Row 1 tests relapse, row 2 tests death before relapse
print(cif_group$Tests)

# Collect the curves into one data frame for plotting.
# Curve names have the form "<group> <cause>", e.g. "AML low risk 1".
cif_curves <- data.frame()
for (curve_name in names(cif_group)) {
  if (curve_name == "Tests") next
  curve <- cif_group[[curve_name]]
  cif_curves <- rbind(cif_curves, data.frame(
    group    = sub(" [12]$", "", curve_name),
    cause    = ifelse(grepl("1$", curve_name), "Relapse", "Death before relapse"),
    time     = curve$time,
    estimate = curve$est
  ))
}
cif_curves$group <- factor(cif_curves$group, levels = levels(analysis_data$disease_group))
cif_curves$cause <- factor(cif_curves$cause, levels = c("Relapse", "Death before relapse"))

cif_group_plot <- ggplot(cif_curves, aes(x = time, y = estimate, colour = group)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ cause) +
  ylim(0, 1) +
  labs(title = "Cumulative incidence of relapse and death before relapse",
       subtitle = "Aalen-Johansen estimates by disease group",
       x = "Days since transplant", y = "Cumulative incidence",
       colour = "Disease group") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
print(cif_group_plot)
ggsave("output/figures/cif_by_group.png", cif_group_plot,
       width = 10, height = 5, dpi = 150, bg = "white")


# Naive 1 - Kaplan-Meier versus Aalen-Johansen for relapse --------------------
# The naive approach censors patients at death, which implicitly assumes they
# could still relapse later. This inflates the estimated relapse risk.

naive_km <- survfit(Surv(time, relapse_event) ~ 1, data = analysis_data)

relapse_comparison <- rbind(
  data.frame(method   = "Aalen-Johansen (competing risks)",
             time     = cif_overall[["1 1"]]$time,
             estimate = cif_overall[["1 1"]]$est),
  data.frame(method   = "1 - Kaplan-Meier (deaths censored)",
             time     = c(0, naive_km$time),
             estimate = c(0, 1 - naive_km$surv))
)

naive_plot <- ggplot(relapse_comparison, aes(x = time, y = estimate, colour = method)) +
  geom_step(linewidth = 1) +
  ylim(0, 1) +
  labs(title = "Estimated probability of relapse",
       subtitle = "1 - Kaplan-Meier overestimates relapse when death is a competing event",
       x = "Days since transplant", y = "Probability of relapse", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
print(naive_plot)
ggsave("output/figures/cif_naive_vs_aalen_johansen.png", naive_plot,
       width = 9, height = 5, dpi = 150, bg = "white")

# The same comparison at one year
relapse_1y_aj    <- timepoints(cif_overall, times = 365)$est[1, 1]
relapse_1y_naive <- 1 - summary(naive_km, times = 365)$surv
cat(sprintf("Relapse by 1 year: Aalen-Johansen = %.3f, 1 - KM = %.3f\n",
            relapse_1y_aj, relapse_1y_naive))
