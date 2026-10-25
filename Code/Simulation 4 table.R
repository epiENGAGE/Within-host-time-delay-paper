library(dplyr)
library(readr)
library(knitr)

# --------------------------------------------------
# Simulation 4 results directory
# --------------------------------------------------
SIM_DIR <- "C:/Users/nj7786/Downloads/WHsim4/Results/Simulations/Simulation 4"

# Read existing simulation summary
results <- readr::read_csv(
  file.path(SIM_DIR, "paper_joint_delay_recovery_summary.csv"),
  show_col_types = FALSE
)

# --------------------------------------------------
# Format manuscript table
# --------------------------------------------------
table_sim4 <- results %>%
  mutate(
    delay_code = factor(
      delay_code,
      levels = c("lp", "ip", "gen", "si")
    ),
    moment = factor(
      moment,
      levels = c("Mean", "SD")
    )
  ) %>%
  arrange(delay_code, moment) %>%
  transmute(
    Delay = case_when(
      moment == "SD" ~ "",
      delay_code == "lp" ~ "Latent period",
      delay_code == "ip" ~ "Incubation period",
      delay_code == "gen" ~ "Generation interval",
      delay_code == "si" ~ "Serial interval"
    ),
    Quantity = as.character(moment),
    Truth = sprintf("%.2f", true_value),
    Bias = sprintf("%.3f", median_bias),
    AE = sprintf("%.3f", median_abs_error),
    RMSE = sprintf("%.3f", rmse),
    CRPS = sprintf("%.3f", mean_crps),
    Width = sprintf("%.3f", median_interval_width),
    Coverage = sprintf("%.2f", coverage_90)
  )

# --------------------------------------------------
# Produce booktabs LaTeX table
# --------------------------------------------------
latex_table <- knitr::kable(
  table_sim4,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  row.names = FALSE,
  align = c("l", "l", rep("r", 7)),
  col.names = c(
    "Delay",
    "Quantity",
    "Truth",
    "Bias",
    "AE",
    "RMSE",
    "CRPS",
    "90\\% width",
    "Coverage"
  ),
  linesep = ""
)

# Save LaTeX table
writeLines(
  as.character(latex_table),
  file.path(SIM_DIR, "table_sim4_joint_delay_recovery.tex")
)

# Also save a CSV version
readr::write_csv(
  table_sim4,
  file.path(SIM_DIR, "table_sim4_joint_delay_recovery.csv")
)

# Display table
print(latex_table)