# BMT DATASET: EXPLORATORY DATA ANALYSIS (EDA)
# Keep this file in the same folder as bmt_preprocessing.R.
# This line runs preprocessing and creates bmt_clean before the EDA begins.
source("bmt_preprocessing.R")

# Patient characteristics and disease-group distribution
table(bmt_clean$group) # count patients in each disease group 
prop.table(table(bmt_clean$group)) * 100 # calculate  percentages


ggplot(bmt_clean, aes(x = group, fill = group)) +
  geom_bar() +
  labs(
    title = "Distribution of Patients by Disease Group",
    x = "Disease Group",
    y = "Number of Patients"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Patient sex
ggplot(bmt_clean, aes(x = z3, fill = z3)) +
  geom_bar() +
  labs(title = "Patient Sex Distribution",
       x = "Patient Sex", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Donor sex
ggplot(bmt_clean, aes(x = z4, fill = z4)) +
  geom_bar() +
  labs(title = "Donor Sex Distribution",
       x = "Donor Sex", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Patient CMV status
ggplot(bmt_clean, aes(x = z5, fill = z5)) +
  geom_bar() +
  labs(title = "Patient CMV Status",
       x = "CMV Status", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Donor CMV status
ggplot(bmt_clean, aes(x = z6, fill = z6)) +
  geom_bar() +
  labs(title = "Donor CMV Status",
       x = "CMV Status", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Hospital distribution
ggplot(bmt_clean, aes(x = z9)) +
  geom_bar(fill = "mediumpurple") +
  labs(title = "Patients by Hospital",
       x = "Hospital Code", y = "Count") +
  theme_minimal()

# MTX usage
ggplot(bmt_clean, aes(x = z10, fill = z10)) +
  geom_bar() +
  labs(title = "MTX Prophylaxis Distribution",
       x = "MTX Usage", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")



# Age and baseline-characteristic exploration
summary(bmt_clean$z1)
summary(bmt_clean$z2)

# Mean ages
mean(bmt_clean$z1)
mean(bmt_clean$z2)

# Standard deviations
sd(bmt_clean$z1)
sd(bmt_clean$z2)



#Patient age distribution with density curve
ggplot(bmt_clean, aes(x = z1)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 15,
    fill = "steelblue",
    color = "white"
  ) +
  geom_density(color = "darkred", linewidth = 1) +
  labs(
    title = "Distribution of Patient Age",
    subtitle = "Histogram with density curve",
    x = "Patient age (years)",
    y = "Density"
  ) +
  theme_minimal()


# Comparing patient and donor age distributions
age_data <- bmt_clean %>%
  select(z1, z2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "person",
    values_to = "age"
  ) %>%
  mutate(
    person = recode(
      person,
      z1 = "Patient",
      z2 = "Donor"
    )
  )

ggplot(age_data, aes(x = age, fill = person, color = person)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  labs(
    title = "Patient and Donor Age Distributions",
    x = "Age (years)",
    y = "Density",
    fill = "Person",
    color = "Person"
  ) +
  theme_minimal()

# Violin plot of patient age by disease group
ggplot(bmt_clean, aes(x = group, y = z1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.5) +
  labs(
    title = "Patient Age by Disease Group",
    x = "Disease group",
    y = "Patient age (years)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#Empirical cumulative distribution of patient age
ggplot(bmt_clean, aes(x = z1)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(
    title = "Cumulative Distribution of Patient Age",
    x = "Patient age (years)",
    y = "Proportion of patients at or below this age"
  ) +
  theme_minimal()

#Patient age against overall follow-up time
ggplot(
  bmt_clean,
  aes(x = z1, y = t1, color = factor(d1))
) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Patient Age and Recorded Follow-up Time",
    subtitle = "Follow-up time shown on a logarithmic scale",
    x = "Patient age (years)",
    y = "Recorded time (days, log scale)",
    color = "Death observed",
    caption = "0 = censored; 1 = death observed"
  ) +
  theme_minimal()

# Correlation heatmap for baseline numeric characteristics
baseline_numeric <- bmt_clean %>%
  select(
    `Patient age` = z1,
    `Donor age` = z2,
    `Waiting time` = z7
  )

correlation_matrix <- cor(
  baseline_numeric,
  use = "complete.obs",
  method = "spearman"
)

correlation_data <- as.data.frame(
  as.table(correlation_matrix)
)

ggplot(
  correlation_data,
  aes(x = Var1, y = Var2, fill = Freq)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 2)), size = 5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    title = "Correlation Between Baseline Characteristics",
    x = NULL,
    y = NULL,
    fill = "Spearman\ncorrelation"
  ) +
  theme_minimal()

ggplot(
  bmt_clean,
  aes(x = z7, fill = group)
) +
  geom_histogram(
    bins = 15,
    alpha = 0.6,
    position = "identity",
    color = "white"
  ) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Waiting Time Before Transplantation",
    x = "Waiting time (days)",
    y = "Number of patients"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Follow-up-time distributions and observed events
summary(bmt_clean$t1)

median(bmt_clean$t1)
mean(bmt_clean$t1)

# Plot the recorded follow-up times
ggplot(bmt_clean, aes(x = t1)) +
  geom_histogram(bins = 20, fill = "cornflowerblue", color = "white") +
  labs(
    title = "Distribution of Overall Survival Follow-up Times",
    x = "Time (days)",
    y = "Number of Patients"
  ) +
  theme_minimal()

summary(bmt_clean$t2)

ggplot(bmt_clean, aes(x = t2)) +
  geom_histogram(bins = 20, fill = "coral", color = "white") +
  labs(
    title = "Distribution of Disease-Free Survival Times",
    x = "Time (days)",
    y = "Number of Patients"
  ) +
  theme_minimal()



event_data <- bmt_clean %>%
  summarise(
    Death = sum(d1 == 1),
    Relapse = sum(d2 == 1),
    `Disease-free event` = sum(d3 == 1),
    `Acute GVHD` = sum(da == 1),
    `Chronic GVHD` = sum(dc == 1),
    `Platelet recovery` = sum(dp == 1)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Event",
    values_to = "Count"
  )

ggplot(event_data, aes(x = reorder(Event, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Observed Events in the BMT Dataset",
    x = "Event",
    y = "Number of Patients"
  ) +
  theme_minimal()

#Distribution of observed death and relapse times
event_times <- bind_rows(
  bmt_clean %>%
    filter(d1 == 1) %>%
    transmute(
      outcome = "Observed death",
      time = t1
    ),
  
  bmt_clean %>%
    filter(d2 == 1) %>%
    transmute(
      outcome = "Observed relapse",
      time = t2
    )
)

ggplot(event_times, aes(x = time, fill = outcome)) +
  geom_density(alpha = 0.35) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(
    title = "Timing of Observed Deaths and Relapses",
    x = "Time (days)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


#Time to platelet recovery among patients with observed recovery
ggplot(
  bmt_clean %>% filter(dp == 1),
  aes(x = tp)
) +
  geom_histogram(
    bins = 15,
    fill = "darkseagreen",
    color = "white"
  ) +
  geom_rug(alpha = 0.5) +
  labs(
    title = "Time to Observed Platelet Recovery",
    subtitle = "Patients with dp = 1",
    x = "Time to platelet recovery (days)",
    y = "Number of patients"
  ) +
  theme_minimal()

#Time to acute and chronic GVHD among patients with observed events
gvhd_times <- bind_rows(
  bmt_clean %>%
    filter(da == 1) %>%
    transmute(
      outcome = "Acute GVHD",
      time = ta
    ),
  
  bmt_clean %>%
    filter(dc == 1) %>%
    transmute(
      outcome = "Chronic GVHD",
      time = tc
    )
)

ggplot(gvhd_times, aes(x = time)) +
  geom_histogram(
    bins = 12,
    fill = "mediumpurple",
    color = "white"
  ) +
  facet_wrap(~ outcome, scales = "free") +
  labs(
    title = "Timing of Observed GVHD Events",
    x = "Time (days)",
    y = "Number of observed events"
  ) +
  theme_minimal()

#Compare observed event proportions by disease group
event_rates <- bmt_clean %>%
  group_by(group) %>%
  summarise(
    Death = mean(d1 == 1),
    Relapse = mean(d2 == 1),
    `Disease-free event` = mean(d3 == 1),
    `Acute GVHD` = mean(da == 1),
    `Chronic GVHD` = mean(dc == 1),
    `Platelet recovery` = mean(dp == 1),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -group,
    names_to = "outcome",
    values_to = "proportion"
  )

ggplot(
  event_rates,
  aes(x = group, y = proportion, color = group)
) +
  geom_point(size = 3) +
  geom_line(aes(group = 1), color = "grey60") +
  facet_wrap(~ outcome, scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  labs(
    title = "Observed Event Proportions by Disease Group",
    x = "Disease group",
    y = "Observed proportion",
    color = "Disease group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Survival analysis
install.packages("survival")
library(survival)
# Overall survival object
surv_obj <- Surv(
  time = bmt_clean$t1,
  event = bmt_clean$d1
)

# Estimate survival curves for each disease group
km_fit <- survfit(
  surv_obj ~ group,
  data = bmt_clean
)

install.packages("survminer")
library(survminer)

#Overall Kaplan–Meier survival curve
km_overall <- survfit(
  Surv(t1, d1) ~ 1,
  data = bmt_clean
)

ggsurvplot(
  km_overall,
  data = bmt_clean,
  conf.int = TRUE,
  risk.table = TRUE,
  censor = TRUE,
  xlab = "Time (days)",
  ylab = "Estimated probability of survival",
  title = "Overall Kaplan–Meier Survival Curve",
  risk.table.height = 0.25,
  ggtheme = theme_minimal()
)

#Compare survival across disease groups
km_group <- survfit(
  Surv(t1, d1) ~ group,
  data = bmt_clean
)

ggsurvplot(
  km_group,
  data = bmt_clean,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  censor = TRUE,
  xlab = "Time (days)",
  ylab = "Estimated survival probability",
  title = "Overall Survival by Disease Group",
  legend.title = "Disease group",
  risk.table.height = 0.25,
  ggtheme = theme_minimal()
)

#Disease-free survival curves by disease group
km_dfs <- survfit(
  Surv(t2, d3) ~ group,
  data = bmt_clean
)

ggsurvplot(
  km_dfs,
  data = bmt_clean,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  censor = TRUE,
  xlab = "Time (days)",
  ylab = "Estimated probability of remaining disease-free",
  title = "Disease-Free Survival by Disease Group",
  legend.title = "Disease group",
  risk.table.height = 0.25,
  ggtheme = theme_minimal()
)

#Cumulative hazard by disease group
ggsurvplot(
  km_group,
  data = bmt_clean,
  fun = "cumhaz",
  conf.int = TRUE,
  censor = TRUE,
  xlab = "Time (days)",
  ylab = "Estimated cumulative hazard",
  title = "Cumulative Hazard by Disease Group",
  legend.title = "Disease group",
  ggtheme = theme_minimal()
)

