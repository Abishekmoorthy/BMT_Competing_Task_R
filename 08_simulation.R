# -----------------------------------------------------------------------------
# 08_simulation.R
#
# Purpose : Simulation study comparing the Aalen-Johansen estimator with the
#           naive 1 - Kaplan-Meier estimator of the probability of relapse.
#
# Design  : Times to relapse, death and censoring are drawn from independent
#           exponential distributions, so the true cumulative incidence of
#           relapse at time tau is known exactly:
#               F1(tau) = l1 / (l1 + l2) * (1 - exp(-(l1 + l2) * tau))
#           where l1 = relapse rate and l2 = death rate.
#
#           Sample size      : 137 (same as the BMT data)
#           Relapse rate     : 0.08
#           Death rate       : 0.02, 0.08, 0.16  (low, medium, high competition)
#           Censoring rate   : 0.02 (low), 0.12 (high)
#           Replications     : 500 for each of the 6 scenarios
#           Time point (tau) : 5
#
# Measures: bias and root mean squared error (RMSE) of each estimator.
# Output  : output/tables/simulation_summary.csv,
#           output/figures/simulation_boxplots.png
# Runtime : about one minute.
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")

set.seed(20261014)

n_patients     <- 137
n_replications <- 500
relapse_rate   <- 0.08
tau            <- 5

scenarios <- expand.grid(
  death_rate = c(0.02, 0.08, 0.16),
  censoring  = c("Low", "High"),
  stringsAsFactors = FALSE
)
scenarios$censor_rate <- ifelse(scenarios$censoring == "Low", 0.02, 0.12)
print(scenarios)


# Run the simulation ----------------------------------------------------------

simulation_rows <- list()

for (s in seq_len(nrow(scenarios))) {

  death_rate  <- scenarios$death_rate[s]
  censor_rate <- scenarios$censor_rate[s]
  true_cif    <- relapse_rate / (relapse_rate + death_rate) *
    (1 - exp(-(relapse_rate + death_rate) * tau))

  aj_estimates    <- numeric(n_replications)
  naive_estimates <- numeric(n_replications)

  for (r in seq_len(n_replications)) {

    # One simulated dataset: the earliest of the three times is observed
    time_relapse <- rexp(n_patients, rate = relapse_rate)
    time_death   <- rexp(n_patients, rate = death_rate)
    time_censor  <- rexp(n_patients, rate = censor_rate)

    observed_time <- pmin(time_relapse, time_death, time_censor)
    status <- ifelse(observed_time == time_relapse, 1,
                     ifelse(observed_time == time_death, 2, 0))

    # Aalen-Johansen estimate at tau
    aj_fit <- cuminc(ftime = observed_time, fstatus = status, cencode = 0)
    aj_estimates[r] <- timepoints(aj_fit, times = tau)$est[1, 1]

    # Naive estimate: deaths treated as censored
    km_fit <- survfit(Surv(observed_time, status == 1) ~ 1)
    naive_estimates[r] <- 1 - summary(km_fit, times = tau, extend = TRUE)$surv
  }

  simulation_rows[[s]] <- data.frame(
    death_rate = death_rate,
    censoring  = scenarios$censoring[s],
    true_cif   = true_cif,
    method     = rep(c("Aalen-Johansen", "1 - Kaplan-Meier"), each = n_replications),
    estimate   = c(aj_estimates, naive_estimates)
  )
  cat("Scenario", s, "of", nrow(scenarios), "complete\n")
}

simulation_results <- do.call(rbind, simulation_rows)
simulation_results$censoring <- factor(simulation_results$censoring, levels = c("Low", "High"))


# Summarise performance -------------------------------------------------------
# bias = average difference from the true value (0 is ideal)
# rmse = typical size of the error (smaller is better)

simulation_summary <- simulation_results %>%
  group_by(death_rate, censoring, method) %>%
  summarise(
    true_cif      = mean(true_cif),
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias          = mean(estimate - true_cif, na.rm = TRUE),
    rmse          = sqrt(mean((estimate - true_cif)^2, na.rm = TRUE)),
    .groups = "drop"
  )
print(as.data.frame(simulation_summary), digits = 3)
write.csv(simulation_summary, "output/tables/simulation_summary.csv", row.names = FALSE)


# Plot ------------------------------------------------------------------------

true_values <- distinct(simulation_results, death_rate, true_cif)

simulation_plot <- ggplot(simulation_results, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_cif),
             colour = "red", linewidth = 0.8) +
  facet_grid(censoring ~ death_rate, labeller = label_both) +
  labs(title = "Estimated probability of relapse at time 5 (500 simulated datasets each)",
       subtitle = "Red line: true value",
       x = NULL, y = "Estimate") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
print(simulation_plot)
ggsave("output/figures/simulation_boxplots.png", simulation_plot,
       width = 10, height = 6, dpi = 150, bg = "white")
