# -----------------------------------------------------------------------------
# app.R - BMT Survival Explorer
#
# Interactive version of the analysis in scripts 02 to 08. It uses the same
# prepared data (01_preprocessing.R) and the same models.
#
# How to run: open this file in RStudio and click "Run App".
#
# A Shiny app has two parts:
#   ui     - the layout of the page: tabs, menus, sliders and output areas
#   server - the R code that fills each output area. Each output is
#            recalculated automatically when an input it uses changes.
# Inputs are read with input$<name>; outputs are created as output$<name>.
# -----------------------------------------------------------------------------

source("01_preprocessing.R")
library(shiny)


# Menu choices (text shown in the app = column name in the data) -------------

explore_choices <- c(
  "Patient age (years)"         = "patient_age",
  "Donor age (years)"           = "donor_age",
  "Waiting time (days)"         = "waiting_days",
  "Disease group"               = "disease_group",
  "Patient sex"                 = "patient_sex",
  "Patient CMV"                 = "patient_cmv",
  "Donor CMV"                   = "donor_cmv",
  "MTX given"                   = "mtx",
  "First event"                 = "event_label"
)

group_choices <- c(
  "Disease group" = "disease_group",
  "Patient sex"   = "patient_sex",
  "Patient CMV"   = "patient_cmv",
  "Donor CMV"     = "donor_cmv",
  "MTX given"     = "mtx"
)


# Results that do not depend on any input are calculated once at start-up ---

# Main regression models (same as 05_cox_model.R and 06_fine_gray_model.R)
cox_relapse <- coxph(Surv(time, relapse_event) ~ disease_group + age10 + patient_sex,
                     data = analysis_data)
cox_summary <- summary(cox_relapse)

fg_covariates <- model.matrix(~ disease_group + age10 + patient_sex, data = analysis_data)[, -1]
fg_relapse <- crr(ftime = analysis_data$time, fstatus = analysis_data$event_status,
                  cov1 = fg_covariates, failcode = 1, cencode = 0)
fg_summary <- summary(fg_relapse)

model_results <- data.frame(
  Variable  = term_labels[colnames(fg_covariates)],
  HR        = cox_summary$conf.int[, "exp(coef)"],
  HR_lower  = cox_summary$conf.int[, "lower .95"],
  HR_upper  = cox_summary$conf.int[, "upper .95"],
  HR_p      = cox_summary$coefficients[, "Pr(>|z|)"],
  SHR       = fg_summary$conf.int[, 1],
  SHR_lower = fg_summary$conf.int[, 3],
  SHR_upper = fg_summary$conf.int[, 4],
  SHR_p     = fg_summary$coef[, 5],
  row.names = NULL
)

ph_test <- cox.zph(cox_relapse)$table

# Relapse curves for the naive vs Aalen-Johansen comparison
cif_all  <- cuminc(analysis_data$time, analysis_data$event_status, cencode = 0)
naive_km <- survfit(Surv(time, relapse_event) ~ 1, data = analysis_data)


# User interface --------------------------------------------------------------

ui <- navbarPage(
  title = "BMT Survival Explorer",

  # Tab 1: overview
  tabPanel("Overview",
    fluidRow(column(8, offset = 1,
      h2("Relapse after bone marrow transplantation"),
      p("137 patients with acute leukaemia who received a bone marrow transplant ",
        "(dataset bmt, R package KMsurv). For each patient the first of three ",
        "outcomes is recorded:"),
      tags$ul(
        tags$li(strong("Relapse"), " - the event of interest"),
        tags$li(strong("Death before relapse"), " - a competing event: after death, relapse cannot occur"),
        tags$li(strong("Censored"), " - follow-up ended before either event")
      ),
      tableOutput("overview_table"),
      p("Use the tabs above in order: data, survival, competing risks, models, simulation.")
    ))
  ),

  # Tab 2: data exploration
  tabPanel("Data",
    sidebarLayout(
      sidebarPanel(
        selectInput("explore_var", "Variable", choices = explore_choices),
        selectInput("explore_group", "Compare by", choices = group_choices)
      ),
      mainPanel(plotOutput("explore_plot", height = "450px"))
    )
  ),

  # Tab 3: Kaplan-Meier
  tabPanel("Survival",
    sidebarLayout(
      sidebarPanel(
        selectInput("km_group", "Compare by", choices = group_choices),
        helpText("Kaplan-Meier estimate of overall survival. The log-rank test ",
                 "checks whether the curves differ between groups.")
      ),
      mainPanel(
        plotOutput("km_plot", height = "450px"),
        verbatimTextOutput("km_test")
      )
    )
  ),

  # Tab 4: competing risks
  tabPanel("Competing risks",
    sidebarLayout(
      sidebarPanel(
        selectInput("cif_group", "Compare by", choices = group_choices),
        sliderInput("cif_day", "Estimate at day", min = 0, max = 2000, value = 365, step = 50),
        helpText("Cumulative incidence (Aalen-Johansen) of relapse and of death before ",
                 "relapse. Gray's test compares the groups.")
      ),
      mainPanel(
        plotOutput("cif_plot", height = "420px"),
        verbatimTextOutput("cif_test"),
        hr(),
        h4("Why not use 1 - Kaplan-Meier for relapse?"),
        plotOutput("naive_plot", height = "380px"),
        textOutput("naive_text")
      )
    )
  ),

  # Tab 5: regression models
  tabPanel("Models",
    fluidRow(column(10, offset = 1,
      h3("Effect of patient characteristics on relapse"),
      p("Cause-specific Cox model (HR): effect on the relapse rate among patients ",
        "still alive and relapse-free. Fine-Gray model (SHR): effect on the overall ",
        "probability of relapse. A value above 1 means higher risk."),
      plotOutput("forest_plot", height = "380px"),
      tableOutput("model_table"),
      h4("Proportional hazards check for the Cox model"),
      p("A p-value above 0.05 means no evidence that the effect changes over time."),
      tableOutput("ph_table")
    ))
  ),

  # Tab 6: simulation
  tabPanel("Simulation",
    sidebarLayout(
      sidebarPanel(
        sliderInput("sim_death", "Death rate (competing event)", 0, 0.30, 0.08, step = 0.02),
        sliderInput("sim_censor", "Censoring rate", 0.02, 0.30, 0.05, step = 0.02),
        sliderInput("sim_reps", "Number of simulated datasets", 50, 300, 100, step = 50),
        actionButton("sim_run", "Run simulation"),
        helpText("Each dataset has 137 patients and a relapse rate of 0.08. ",
                 "The red line is the true probability of relapse by time 5.")
      ),
      mainPanel(
        plotOutput("sim_plot", height = "400px"),
        tableOutput("sim_table")
      )
    )
  )
)


# Server ----------------------------------------------------------------------

server <- function(input, output, session) {

  # Tab 1: overview table of first events by disease group
  output$overview_table <- renderTable({
    as.data.frame.matrix(table(analysis_data$disease_group, analysis_data$event_label))
  }, rownames = TRUE)

  # Tab 2: histogram for numeric variables, bar chart for categories
  output$explore_plot <- renderPlot({
    var   <- input$explore_var
    group <- input$explore_group
    var_name   <- names(explore_choices)[explore_choices == var]
    group_name <- names(group_choices)[group_choices == group]

    if (is.numeric(analysis_data[[var]])) {
      ggplot(analysis_data, aes(x = .data[[group]], y = .data[[var]], fill = .data[[group]])) +
        geom_boxplot() +
        labs(x = group_name, y = var_name) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")
    } else {
      ggplot(analysis_data, aes(x = .data[[group]], fill = .data[[var]])) +
        geom_bar(position = "dodge") +
        labs(x = group_name, y = "Number of patients", fill = var_name) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")
    }
  })

  # Tab 3: Kaplan-Meier curves for the chosen grouping
  km_formula <- reactive({
    as.formula(paste("Surv(os_days, died) ~", input$km_group))
  })

  output$km_plot <- renderPlot({
    fit <- survfit(km_formula(), data = analysis_data)
    group_levels <- levels(analysis_data[[input$km_group]])
    plot(fit, col = 2:4, lwd = 2, mark.time = TRUE,
         xlab = "Days since transplant", ylab = "Survival probability",
         main = "Overall survival (Kaplan-Meier)")
    legend("topright", legend = group_levels, col = 2:4, lwd = 2, bty = "n")
  })

  output$km_test <- renderPrint({
    survdiff(km_formula(), data = analysis_data)
  })

  # Tab 4: cumulative incidence for the chosen grouping
  cif_fit <- reactive({
    cuminc(ftime   = analysis_data$time,
           fstatus = analysis_data$event_status,
           group   = analysis_data[[input$cif_group]],
           cencode = 0)
  })

  output$cif_plot <- renderPlot({
    fit <- cif_fit()
    # Curve names look like "<group> <cause>", e.g. "ALL 1" (1 = relapse, 2 = death)
    curve_names <- names(fit)[names(fit) != "Tests"]
    curve_group <- sub(" [12]$", "", curve_names)
    is_relapse  <- grepl(" 1$", curve_names)

    plot(fit,
         curvlab = paste(curve_group, ifelse(is_relapse, "- relapse", "- death")),
         color = match(curve_group, unique(curve_group)) + 1,
         lty   = ifelse(is_relapse, 1, 2),
         lwd   = 2,
         xlab  = "Days since transplant", ylab = "Cumulative incidence",
         main  = "Relapse (solid) and death before relapse (dashed)")
    abline(v = input$cif_day, lty = 3)
  })

  output$cif_test <- renderPrint({
    cat("Gray's test (row 1 = relapse, row 2 = death before relapse)\n")
    print(round(cif_fit()$Tests, 4))
    cat("\nEstimates at day", input$cif_day, "\n")
    print(round(timepoints(cif_fit(), times = input$cif_day)$est, 3))
  })

  output$naive_plot <- renderPlot({
    relapse_curve <- cif_all[["1 1"]]
    plot(relapse_curve$time, relapse_curve$est, type = "s", lwd = 2, col = "steelblue",
         ylim = c(0, 0.5), xlab = "Days since transplant", ylab = "Probability of relapse",
         main = "Relapse: Aalen-Johansen vs 1 - Kaplan-Meier (all patients)")
    lines(c(0, naive_km$time), 1 - c(1, naive_km$surv), type = "s", lwd = 2, col = "firebrick")
    abline(v = input$cif_day, lty = 3)
    legend("topleft", bty = "n", lwd = 2, col = c("steelblue", "firebrick"),
           legend = c("Aalen-Johansen (correct)", "1 - Kaplan-Meier (deaths censored)"))
  })

  output$naive_text <- renderText({
    aj    <- timepoints(cif_all, times = input$cif_day)$est[1, 1]
    naive <- 1 - summary(naive_km, times = input$cif_day, extend = TRUE)$surv
    sprintf("At day %d: Aalen-Johansen = %.3f, 1 - Kaplan-Meier = %.3f.",
            input$cif_day, aj, naive)
  })

  # Tab 5: regression models (fitted once at start-up)
  output$forest_plot <- renderPlot({
    forest <- rbind(
      data.frame(model = "Cause-specific Cox (HR)", variable = model_results$Variable,
                 ratio = model_results$HR, lower = model_results$HR_lower,
                 upper = model_results$HR_upper),
      data.frame(model = "Fine-Gray (SHR)", variable = model_results$Variable,
                 ratio = model_results$SHR, lower = model_results$SHR_lower,
                 upper = model_results$SHR_upper)
    )
    ggplot(forest, aes(x = variable, y = ratio, colour = model)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
      scale_y_log10() +
      coord_flip() +
      labs(x = NULL, y = "Ratio with 95% confidence interval (log scale)", colour = NULL) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })

  output$model_table <- renderTable(model_results, digits = 3)

  output$ph_table <- renderTable({
    data.frame(Term = rownames(ph_test), Chi_square = ph_test[, "chisq"], p_value = ph_test[, "p"])
  }, digits = 3)

  # Tab 6: simulation, run when the button is pressed
  sim_results <- eventReactive(input$sim_run, {
    n_patients   <- 137
    relapse_rate <- 0.08
    tau          <- 5
    true_cif <- relapse_rate / (relapse_rate + input$sim_death) *
      (1 - exp(-(relapse_rate + input$sim_death) * tau))

    aj_estimates    <- numeric(input$sim_reps)
    naive_estimates <- numeric(input$sim_reps)

    withProgress(message = "Simulating", value = 0, {
      for (r in seq_len(input$sim_reps)) {
        time_relapse <- rexp(n_patients, relapse_rate)
        time_death   <- rexp(n_patients, input$sim_death)
        time_censor  <- rexp(n_patients, input$sim_censor)
        observed_time <- pmin(time_relapse, time_death, time_censor)
        status <- ifelse(observed_time == time_relapse, 1,
                         ifelse(observed_time == time_death, 2, 0))

        aj_estimates[r] <- timepoints(cuminc(observed_time, status, cencode = 0), tau)$est[1, 1]
        km <- survfit(Surv(observed_time, status == 1) ~ 1)
        naive_estimates[r] <- 1 - summary(km, times = tau, extend = TRUE)$surv

        incProgress(1 / input$sim_reps)
      }
    })

    data.frame(method   = rep(c("Aalen-Johansen", "1 - Kaplan-Meier"), each = input$sim_reps),
               estimate = c(aj_estimates, naive_estimates),
               true_cif = true_cif)
  }, ignoreNULL = FALSE)

  output$sim_plot <- renderPlot({
    sims <- sim_results()
    boxplot(estimate ~ method, data = sims, col = c("firebrick", "steelblue"),
            xlab = "", ylab = "Estimated probability of relapse at time 5",
            main = "Estimates from the simulated datasets")
    abline(h = sims$true_cif[1], col = "red", lwd = 2)
  })

  output$sim_table <- renderTable({
    sim_results() %>%
      group_by(method) %>%
      summarise(true_value    = mean(true_cif),
                mean_estimate = mean(estimate, na.rm = TRUE),
                bias          = mean(estimate - true_cif, na.rm = TRUE),
                rmse          = sqrt(mean((estimate - true_cif)^2, na.rm = TRUE)))
  }, digits = 4)
}

shinyApp(ui, server)
