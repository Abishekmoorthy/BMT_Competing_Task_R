# Competing Risks Analysis of Relapse after Bone Marrow Transplantation

Analysis of relapse, with death before relapse as a competing event, in 137 patients
with acute leukaemia who received a bone marrow transplant (`bmt` dataset, R package KMsurv).

A plain-language explanation of every method is in **`GUIDE.html`** (open it in any browser).

## How to run (in this order)

1. Install **R** (https://cran.r-project.org) and **RStudio** (https://posit.co/download/rstudio-desktop/).
2. Unzip this folder and double-click **`BMT_Project.Rproj`**. RStudio opens in the project folder.
3. Open **`run_all.R`** and click **Source** (top-right of the script window).
   - First run only: missing packages are installed automatically (needs internet, a few minutes).
   - Results print in the Console, plots appear in the Plots pane, and everything is saved in `output/`.
   - Total run time is about two minutes, most of it the simulation.
4. Open **`app.R`** and click **Run App** to start the interactive Shiny app.
   The app works on its own; you do not need to run step 3 first.

## Project structure

| File | Content |
|---|---|
| `BMT_Project.Rproj` | RStudio project file (sets the working folder) |
| `01_preprocessing.R` | Packages, data loading, readable variable names, competing-risks outcome |
| `02_eda.R` | Exploratory analysis: summary table and six plots |
| `03_kaplan_meier.R` | Kaplan-Meier overall survival by disease group, log-rank test |
| `04_cumulative_incidence.R` | Aalen-Johansen cumulative incidence, Gray's test, comparison with 1 - KM |
| `05_cox_model.R` | Cause-specific Cox model for relapse, proportional hazards test |
| `06_fine_gray_model.R` | Fine-Gray model for relapse |
| `07_model_comparison.R` | Cox vs Fine-Gray table and forest plot |
| `08_simulation.R` | Simulation study: Aalen-Johansen vs 1 - Kaplan-Meier |
| `run_all.R` | Runs scripts 01 to 08 in order |
| `app.R` | Interactive Shiny app (6 tabs) |
| `GUIDE.html` | Explanation of all concepts and how to read the results |

Each script can also be run on its own; it loads the data automatically if needed.

## Troubleshooting

- `cannot open file '01_preprocessing.R'`: RStudio is not in the project folder. Open the project
  with `BMT_Project.Rproj`, or use *Session > Set Working Directory > To Source File Location*.
- `there is no package called ...`: run `install.packages("<name>")` in the Console.
- Plots look squashed: make the Plots pane larger, or open the saved PNG in `output/figures/`.

## References

- Klein, J. P. and Moeschberger, M. L. (2003). *Survival Analysis: Techniques for Censored
  and Truncated Data*, 2nd ed. Springer. (Source of the `bmt` data.)
- Kaplan, E. L. and Meier, P. (1958). Nonparametric estimation from incomplete observations.
  *Journal of the American Statistical Association*, 53, 457-481.
- Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical
  Society, Series B*, 34, 187-220.
- Gray, R. J. (1988). A class of K-sample tests for comparing the cumulative incidence of a
  competing risk. *Annals of Statistics*, 16, 1141-1154.
- Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution
  of a competing risk. *Journal of the American Statistical Association*, 94, 496-509.
- Putter, H., Fiocco, M. and Geskus, R. B. (2007). Tutorial in biostatistics: competing risks
  and multi-state models. *Statistics in Medicine*, 26, 2389-2430.
- R packages: survival, cmprsk, survminer, KMsurv, dplyr, tidyr, ggplot2, shiny.
