# Get GT estimates for Covid
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

OLD_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Noncens"

NEW_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

OUT_DIR <- file.path(NEW_DIR, "generation_time_estimates_three_models")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

find_censored_fit <- function(strain_slug) {
  env_name <- if (strain_slug == "sgtf") "CENS_SGTF_FIT" else "CENS_NSGTF_FIT"
  env_path <- Sys.getenv(env_name, unset = "")
  
  if (nzchar(env_path) && file.exists(env_path)) {
    return(env_path)
  }
  
  hits <- list.files(
    NEW_DIR,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  hits <- hits[
    !grepl("_loo\\.rds$|_waic\\.rds$|loo\\.rds$|waic\\.rds$", basename(hits), ignore.case = TRUE)
  ]
  
  hits <- hits[
    grepl("fit|cens|censored", basename(hits), ignore.case = TRUE)
  ]
  
  if (strain_slug == "sgtf") {
    hits <- hits[
      grepl("sgtf", basename(hits), ignore.case = TRUE) &
        !grepl("nsgtf|non[-_ ]?sgtf", basename(hits), ignore.case = TRUE)
    ]
  } else {
    hits <- hits[
      grepl("nsgtf|non[-_ ]?sgtf", basename(hits), ignore.case = TRUE)
    ]
  }
  
  if (length(hits) != 1) {
    message("\nCandidate censored fit files for ", strain_slug, ":")
    print(hits)
    stop(
      "\nCould not uniquely identify censored ", strain_slug, " fit file.\n",
      "Set manually before running:\n",
      "Sys.setenv(CENS_SGTF_FIT='/full/path/to/censored_sgtf_fit.rds')\n",
      "Sys.setenv(CENS_NSGTF_FIT='/full/path/to/censored_nsgtf_fit.rds')\n"
    )
  }
  
  hits[1]
}

fit_registry <- tibble::tribble(
  ~model,                    ~strain,     ~path,
  "Lognormal",               "SGTF",      file.path(OLD_DIR, "fit_park_lognormal_fast_sgtf.rds"),
  "Lognormal",               "non-SGTF",  file.path(OLD_DIR, "fit_park_lognormal_fast_nsgtf.rds"),
  "WH-informed",             "SGTF",      file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds"),
  "WH-informed",             "non-SGTF",  file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds"),
  "Censored WH-informed",    "SGTF",      find_censored_fit("sgtf"),
  "Censored WH-informed",    "non-SGTF",  find_censored_fit("nsgtf")
)

missing <- fit_registry %>% filter(!file.exists(path))

print(fit_registry)

extract_gt_draws <- function(path, model, strain) {
  fit <- readRDS(path)
  
  vars <- posterior::variables(fit$draws(inc_warmup = FALSE))
  
  needed <- c("mean_gen", "sd_gen")
  missing_vars <- setdiff(needed, vars)
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variables in ", model, " / ", strain, ":\n",
      paste(missing_vars, collapse = ", "), "\n",
      "File: ", path
    )
  }
  
  posterior::as_draws_df(
    fit$draws(variables = needed, inc_warmup = FALSE)
  ) %>%
    as_tibble() %>%
    transmute(
      model = model,
      strain = strain,
      mean_gen = as.numeric(mean_gen),
      sd_gen = as.numeric(sd_gen)
    )
}

gt_draws <- fit_registry %>%
  rowwise() %>%
  do(
    extract_gt_draws(
      path = .$path,
      model = .$model,
      strain = .$strain
    )
  ) %>%
  ungroup() %>%
  mutate(
    strain_label = recode(
      strain,
      "SGTF" = "Omicron",
      "non-SGTF" = "Delta"
    )
  )

gt_summary <- gt_draws %>%
  group_by(strain_label, strain, model) %>%
  summarise(
    mean_gen_median = median(mean_gen, na.rm = TRUE),
    mean_gen_mean = mean(mean_gen, na.rm = TRUE),
    mean_gen_q025 = quantile(mean_gen, 0.025, na.rm = TRUE),
    mean_gen_q975 = quantile(mean_gen, 0.975, na.rm = TRUE),
    
    sd_gen_median = median(sd_gen, na.rm = TRUE),
    sd_gen_mean = mean(sd_gen, na.rm = TRUE),
    sd_gen_q025 = quantile(sd_gen, 0.025, na.rm = TRUE),
    sd_gen_q975 = quantile(sd_gen, 0.975, na.rm = TRUE),
    
    n_draws = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(strain_label, model)

gt_summary_rounded <- gt_summary %>%
  mutate(
    across(
      c(
        mean_gen_median, mean_gen_mean, mean_gen_q025, mean_gen_q975,
        sd_gen_median, sd_gen_mean, sd_gen_q025, sd_gen_q975
      ),
      ~ round(.x, 2)
    )
  )

print(gt_summary_rounded)

write_csv(
  gt_draws,
  file.path(OUT_DIR, "generation_time_draws_three_models.csv")
)

write_csv(
  gt_summary,
  file.path(OUT_DIR, "generation_time_summary_three_models_full_precision.csv")
)

write_csv(
  gt_summary_rounded,
  file.path(OUT_DIR, "generation_time_summary_three_models_rounded.csv")
)

# Convenient paragraph-ready table
paragraph_table <- gt_summary_rounded %>%
  select(
    Strain = strain_label,
    Model = model,
    `Mean GT median` = mean_gen_median,
    `Mean GT 95% CrI lower` = mean_gen_q025,
    `Mean GT 95% CrI upper` = mean_gen_q975,
    `SD GT median` = sd_gen_median,
    `SD GT 95% CrI lower` = sd_gen_q025,
    `SD GT 95% CrI upper` = sd_gen_q975
  )

print(paragraph_table)

write_csv(
  paragraph_table,
  file.path(OUT_DIR, "generation_time_paragraph_table.csv")
)
