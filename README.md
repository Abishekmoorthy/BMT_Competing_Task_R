# Competing Risks in Bone Marrow Transplantation

## Dataset

This project uses the `bmt` dataset from the `KMsurv` R package. The dataset contains information on 137 patients who underwent bone marrow transplantation and belong to three disease groups:

- Acute Lymphoblastic Leukemia (ALL)
- AML Low-Risk
- AML High-Risk

The data include follow-up times and information about events such as leukemia relapse, death, and censored observations.

The main focus is to understand the probability of leukemia relapse when death can occur before relapse. Since death prevents a future relapse, it is considered a competing event. Accounting for these competing events is important for obtaining more appropriate estimates of relapse probability.

## Expected Outcomes

The project is expected to provide a clearer understanding of how competing risks affect survival analysis and the estimation of relapse probabilities.

By comparing different statistical approaches, the analysis aims to identify differences between standard survival estimates and methods that account for competing events. The results are also expected to illustrate the strengths and limitations of different estimation and regression methods.

A Monte Carlo simulation will be used to compare the performance of the estimation methods across 500 replications, considering measures such as bias and root mean squared error (RMSE).

## R Libraries

The following R libraries will be used for data analysis, statistical modelling, visualisation, and application development.

- `KMsurv` – Provides the bone marrow transplantation dataset.
- `survival` – Supports Kaplan–Meier estimation and Cox proportional hazards regression.
- `cmprsk` – Supports cumulative incidence functions and Fine–Gray regression.
- `dplyr` – Used for data filtering, transformation, and summarisation.
- `tidyr` – Helps organise and reshape data.
- `ggplot2` – Used to create statistical plots and visualisations.
- `broom` – Helps organise statistical model outputs.
- `scales` – Supports formatting of plot axes and labels.
- `shiny` – Used to develop an interactive web application for exploring the analysis results.

## Shiny Application

A Shiny application is planned to make the analysis easier to explore without requiring users to work directly with the R code.

The application will present relevant statistical results through interactive tables and plots. It is intended to help users explore survival estimates, cumulative incidence, competing-risk model results, and comparisons between the statistical methods.

The goal is to make the findings more accessible and provide a simple way to understand how competing events influence the estimated probability of relapse.
