library(posterior)
library(dplyr)
library(tibble)

OUT_DIR <- path.expand(
  "~/Within-host-time-delay-paper/Results/Covid/Cens/lognormal"
)

# Load fitted models
fit_omicron <- readRDS(
  file.path(OUT_DIR, "fit_censored_lognormal_sgtf.rds")
)

fit_delta <- readRDS(
  file.path(OUT_DIR, "fit_censored_lognormal_nsgtf.rds")
)

# Function to extract posterior median + 95% CrI for mean generation time
get_gt_summary <- function(fit, variant) {
  
  draws <- posterior::as_draws_df(
    fit$draws(variables = "mean_gen")
  )
  
  tibble(
    variant = variant,
    median = median(draws$mean_gen),
    lower_95 = quantile(draws$mean_gen, 0.025),
    upper_95 = quantile(draws$mean_gen, 0.975)
  )
}

gt_summary <- bind_rows(
  get_gt_summary(fit_omicron, "Omicron"),
  get_gt_summary(fit_delta, "Delta")
)

print(gt_summary)

# Optional: print directly in manuscript-style wording
for (i in seq_len(nrow(gt_summary))) {
  cat(
    gt_summary$variant[i], ": ",
    sprintf("%.2f", gt_summary$median[i]),
    " days (95% CrI: ",
    sprintf("%.2f", gt_summary$lower_95[i]),
    "--",
    sprintf("%.2f", gt_summary$upper_95[i]),
    ")\n",
    sep = ""
  )
}