# -----------------------------------------------------------------------------
# run_all.R
#
# Runs the full analysis in order. Open BMT_Project.Rproj first so that the
# working directory is the project folder, then source this file.
# Results are printed to the console and saved in the output/ folder.
#
# The dashboard is started separately: open dashboard.Rmd and click "Run Document".
# -----------------------------------------------------------------------------

source("01_preprocessing.R")       # data preparation
source("02_eda.R")                 # exploratory analysis
source("03_kaplan_meier.R")        # overall survival by disease group
source("04_cumulative_incidence.R")# relapse vs death (competing risks)
source("05_cox_model.R")           # cause-specific Cox model
source("06_fine_gray_model.R")     # Fine-Gray model
source("07_model_comparison.R")    # Cox vs Fine-Gray
source("08_simulation.R")          # simulation study

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("Analysis complete. See the output/ folder.\n")
