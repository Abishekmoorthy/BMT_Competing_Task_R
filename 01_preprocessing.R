# -----------------------------------------------------------------------------
# 01_preprocessing.R
#
# Purpose : Load the bone marrow transplant (BMT) data from the KMsurv
#           package, recode the variables and create the competing-risks
#           outcome used in the modelling scripts.
# Creates : bmt_clean      - data with descriptive names and labelled factors
#           analysis_data  - bmt_clean plus the outcome variables for modelling
# -----------------------------------------------------------------------------


# Packages --------------------------------------------------------------------

required_packages <- c("KMsurv", "survival", "cmprsk", "survminer",
                       "dplyr", "tidyr", "ggplot2", "shiny")

not_installed <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(not_installed) > 0) {
  install.packages(not_installed)
}

suppressPackageStartupMessages({
  library(KMsurv)     # source of the bmt dataset
  library(survival)   # Kaplan-Meier and Cox models
  library(cmprsk)     # cumulative incidence and Fine-Gray model
  library(survminer)  # survival plots
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Folders for saved figures and tables
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)


# Load data -------------------------------------------------------------------

data("bmt", package = "KMsurv")

cat("Rows:", nrow(bmt), " Columns:", ncol(bmt), "\n")
cat("Missing values:", sum(is.na(bmt)), " Duplicated rows:", sum(duplicated(bmt)), "\n")


# Rename and recode -----------------------------------------------------------
# Original variable names (Klein & Moeschberger, 2003):
#   group       disease group: 1 = ALL, 2 = AML low risk, 3 = AML high risk
#   t1, d1      time to death or end of follow-up (days), death indicator
#   t2, d2      disease-free survival time (days), relapse indicator
#   d3          disease-free survival indicator (relapse or death)
#   ta, da      acute GVHD time and indicator
#   tc, dc      chronic GVHD time and indicator
#   tp, dp      platelet recovery time and indicator
#   z1, z2      patient age, donor age (years)
#   z3, z4      patient sex, donor sex (1 = male)
#   z5, z6      patient CMV status, donor CMV status (1 = positive)
#   z7          waiting time from diagnosis to transplant (days)
#   z8          FAB classification grade 4 or 5 and AML (1 = yes)
#   z9          hospital
#   z10         MTX used as GVHD prophylaxis (1 = yes)

bmt_clean <- bmt %>%
  rename(
    disease_group = group,
    os_days       = t1,  died      = d1,
    dfs_days      = t2,  relapsed  = d2,  dfs_event = d3,
    agvhd_days    = ta,  agvhd     = da,
    cgvhd_days    = tc,  cgvhd     = dc,
    platelet_days = tp,  platelet_recovered = dp,
    patient_age   = z1,  donor_age = z2,
    patient_sex   = z3,  donor_sex = z4,
    patient_cmv   = z5,  donor_cmv = z6,
    waiting_days  = z7,  fab       = z8,
    hospital      = z9,  mtx       = z10
  ) %>%
  mutate(
    disease_group = factor(disease_group, levels = 1:3,
                           labels = c("ALL", "AML low risk", "AML high risk")),
    patient_sex = factor(patient_sex, levels = 0:1, labels = c("Female", "Male")),
    donor_sex   = factor(donor_sex,   levels = 0:1, labels = c("Female", "Male")),
    patient_cmv = factor(patient_cmv, levels = 0:1, labels = c("Negative", "Positive")),
    donor_cmv   = factor(donor_cmv,   levels = 0:1, labels = c("Negative", "Positive")),
    fab         = factor(fab,         levels = 0:1, labels = c("No", "Yes")),
    mtx         = factor(mtx,         levels = 0:1, labels = c("No", "Yes")),
    hospital    = factor(hospital,    levels = 1:4,
                         labels = c("Ohio State", "Alfred", "St. Vincent", "Hahnemann"))
  )


# Competing-risks outcome -----------------------------------------------------
# The event of interest is relapse. Death without a prior relapse prevents
# relapse from ever being observed, so it is treated as a competing event.
#   event_status = 0  censored
#   event_status = 1  relapse
#   event_status = 2  death before relapse
# Time to the first of these events is the disease-free survival time (t2).

analysis_data <- bmt_clean %>%
  mutate(
    time = dfs_days,
    event_status = case_when(
      relapsed == 1 ~ 1,
      died == 1     ~ 2,
      TRUE          ~ 0
    ),
    event_label = factor(event_status, levels = 0:2,
                         labels = c("Censored", "Relapse", "Death before relapse")),
    relapse_event = as.integer(event_status == 1),
    death_event   = as.integer(event_status == 2),
    age10       = patient_age / 10,   # hazard ratios are reported per 10 years
    donor_age10 = donor_age / 10
  )

print(table(analysis_data$event_label))

# Every relapse or death should agree with the original d3 indicator
cat("Outcome coding agrees with d3:",
    all((analysis_data$event_status > 0) == (analysis_data$dfs_event == 1)), "\n")


# Readable names for model terms (used in results tables) --------------------

term_labels <- c(
  "disease_groupAML low risk"  = "AML low risk (vs ALL)",
  "disease_groupAML high risk" = "AML high risk (vs ALL)",
  "age10"                      = "Patient age, per 10 years",
  "patient_sexMale"            = "Male patient (vs female)",
  "donor_age10"                = "Donor age, per 10 years",
  "donor_sexMale"              = "Male donor (vs female)",
  "patient_cmvPositive"        = "Patient CMV positive",
  "donor_cmvPositive"          = "Donor CMV positive",
  "mtxYes"                     = "MTX given",
  "fabYes"                     = "FAB grade 4/5 (AML)"
)
