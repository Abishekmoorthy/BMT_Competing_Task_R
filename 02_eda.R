# -----------------------------------------------------------------------------
# 02_eda.R
#
# Purpose : Exploratory data analysis of the BMT cohort.
#           1. Summary table of baseline characteristics by disease group
#           2. Disease groups and patient/donor characteristics
#           3. Age and waiting time
#           4. Observed events and their timing
# Output  : output/tables/baseline_summary.csv, output/figures/eda_*.png
# -----------------------------------------------------------------------------

if (!exists("analysis_data")) source("01_preprocessing.R")


# 1. Baseline summary table ---------------------------------------------------
# One row per disease group: number of patients, ages, waiting time and the
# proportion of patients in each category.

baseline_summary <- analysis_data %>%
  group_by(disease_group) %>%
  summarise(
    patients          = n(),
    median_age        = median(patient_age),
    median_donor_age  = median(donor_age),
    median_wait_days  = median(waiting_days),
    pct_male          = round(100 * mean(patient_sex == "Male"), 1),
    pct_cmv_positive  = round(100 * mean(patient_cmv == "Positive"), 1),
    pct_mtx           = round(100 * mean(mtx == "Yes"), 1),
    pct_relapse       = round(100 * mean(event_status == 1), 1),
    pct_death_first   = round(100 * mean(event_status == 2), 1)
  )
print(as.data.frame(baseline_summary))
write.csv(baseline_summary, "output/tables/baseline_summary.csv", row.names = FALSE)


# 2. Disease groups and characteristics ---------------------------------------

group_plot <- ggplot(bmt_clean, aes(x = disease_group, fill = disease_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.4) +
  labs(title = "Patients by disease group", x = NULL, y = "Number of patients") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
print(group_plot)
ggsave("output/figures/eda_1_disease_groups.png", group_plot,
       width = 7, height = 5, dpi = 150, bg = "white")

# Categorical characteristics shown together in one faceted plot
characteristics <- bmt_clean %>%
  transmute(
    `Patient sex` = as.character(patient_sex),
    `Donor sex`   = as.character(donor_sex),
    `Patient CMV` = as.character(patient_cmv),
    `Donor CMV`   = as.character(donor_cmv),
    `Hospital`    = as.character(hospital),
    `MTX given`   = as.character(mtx)
  ) %>%
  pivot_longer(everything(), names_to = "characteristic", values_to = "value")

characteristics_plot <- ggplot(characteristics, aes(x = value)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~ characteristic, scales = "free_x") +
  labs(title = "Patient and donor characteristics", x = NULL, y = "Number of patients") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
print(characteristics_plot)
ggsave("output/figures/eda_2_characteristics.png", characteristics_plot,
       width = 10, height = 6, dpi = 150, bg = "white")


# 3. Age and waiting time -----------------------------------------------------

age_long <- bmt_clean %>%
  select(disease_group, Patient = patient_age, Donor = donor_age) %>%
  pivot_longer(c(Patient, Donor), names_to = "person", values_to = "age")

age_plot <- ggplot(age_long, aes(x = disease_group, y = age, fill = person)) +
  geom_boxplot() +
  labs(title = "Patient and donor age by disease group",
       x = NULL, y = "Age (years)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
print(age_plot)
ggsave("output/figures/eda_3_age.png", age_plot,
       width = 8, height = 5, dpi = 150, bg = "white")

wait_plot <- ggplot(bmt_clean, aes(x = disease_group, y = waiting_days, fill = disease_group)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Waiting time from diagnosis to transplant",
       x = NULL, y = "Days (log scale)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
print(wait_plot)
ggsave("output/figures/eda_4_waiting_time.png", wait_plot,
       width = 7, height = 5, dpi = 150, bg = "white")


# 4. Observed events ----------------------------------------------------------
# First event per patient: relapse, death before relapse, or censored

event_plot <- ggplot(analysis_data, aes(x = disease_group, fill = event_label)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(title = "First observed event by disease group",
       x = NULL, y = "Share of patients", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
print(event_plot)
ggsave("output/figures/eda_5_first_event.png", event_plot,
       width = 8, height = 5, dpi = 150, bg = "white")

# Timing of relapses and deaths before relapse
timing_plot <- ggplot(filter(analysis_data, event_status > 0),
                      aes(x = time, fill = event_label)) +
  geom_histogram(bins = 20, colour = "white") +
  facet_wrap(~ event_label, scales = "free_y") +
  labs(title = "Timing of relapses and deaths before relapse",
       x = "Days since transplant", y = "Number of patients") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
print(timing_plot)
ggsave("output/figures/eda_6_event_timing.png", timing_plot,
       width = 9, height = 5, dpi = 150, bg = "white")

