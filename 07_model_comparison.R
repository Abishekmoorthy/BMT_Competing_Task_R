# -----------------------------------------------------------------------------
# 07_model_comparison.R
#
# Purpose : Compare the cause-specific Cox model (05) with the Fine-Gray
#           model (06) for relapse, as a table and a forest plot.
# Reading : HR  - effect on the relapse rate among patients still at risk
#           SHR - effect on the overall probability of relapse
#           The two can differ when a covariate also affects the competing
#           event (death before relapse).
# Output  : output/tables/model_comparison.csv,
#           output/figures/model_forest_plot.png
# -----------------------------------------------------------------------------

if (!exists("cox_results")) source("05_cox_model.R")
if (!exists("fg_results"))  source("06_fine_gray_model.R")


# Side-by-side table ----------------------------------------------------------
# Both models use the same covariates, so the rows match by term

model_comparison <- merge(
  cox_results[, c("term", "variable", "HR", "lower_95", "upper_95", "p_value")],
  fg_results[,  c("term", "SHR", "lower_95", "upper_95", "p_value")],
  by = "term",
  suffixes = c("_cox", "_fine_gray"),
  sort = FALSE
)
print(model_comparison)
write.csv(model_comparison, "output/tables/model_comparison.csv", row.names = FALSE)


# Forest plot -----------------------------------------------------------------

forest_data <- rbind(
  data.frame(model    = "Cause-specific Cox (HR)",
             variable = cox_results$variable,
             ratio    = cox_results$HR,
             lower    = cox_results$lower_95,
             upper    = cox_results$upper_95),
  data.frame(model    = "Fine-Gray (SHR)",
             variable = fg_results$variable,
             ratio    = fg_results$SHR,
             lower    = fg_results$lower_95,
             upper    = fg_results$upper_95)
)

forest_plot <- ggplot(forest_data, aes(x = variable, y = ratio, colour = model)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5)) +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Covariate effects on relapse",
       subtitle = "Estimates with 95% confidence intervals (log scale)",
       x = NULL, y = "Hazard ratio / subdistribution hazard ratio", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
print(forest_plot)
ggsave("output/figures/model_forest_plot.png", forest_plot,
       width = 9, height = 5, dpi = 150, bg = "white")
