# Covid moments plot
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(tibble)
})

# Paths
OLD_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Noncens"

NEW_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

FIG_DIR <- file.path(NEW_DIR, "figure3_generation_time_three_models")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Labels / settings
nice_model_levels <- c(
  "Lognormal",
  "Mechanistic",
  "Censored WH-informed"
)

model_plot_labels <- c(
  "Lognormal" = "Lognormal",
  "Mechanistic" = "WH-informed",
  "Censored WH-informed" = "Censored WH-informed"
)

nice_strain_levels <- c("SGTF", "non-SGTF")

strain_plot_labels <- c(
  "SGTF" = "Omicron",
  "non-SGTF" = "Delta"
)

set.seed(20240513)

# Helpers
theme_paper <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "plain"),
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90")
    )
}

save_both <- function(plot, filename_stem, width, height, dpi = 450) {
  pdf_file <- file.path(FIG_DIR, paste0(filename_stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename_stem, ".png"))
  
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = dpi)
  
  message("Wrote ", pdf_file)
  message("Wrote ", png_file)
}

read_fit <- function(path, model, strain) {
  if (!file.exists(path)) {
    stop("Missing fit: ", path)
  }
  
  list(
    fit = readRDS(path),
    path = path,
    model = model,
    strain = strain
  )
}

fit_variables <- function(fit_obj) {
  tryCatch(
    posterior::variables(fit_obj$draws(inc_warmup = FALSE)),
    error = function(e) character()
  )
}

variable_regex <- function(variable) {
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", variable)
  paste0("^", escaped, "(\\[|$)")
}

has_base_variable <- function(fit_obj, variable) {
  vars <- fit_variables(fit_obj)
  any(str_detect(vars, variable_regex(variable)))
}

available_base_variables <- function(fit_obj, variables) {
  variables[vapply(variables, function(v) has_base_variable(fit_obj, v), logical(1))]
}

try_draws_df <- function(fit_obj, variables) {
  have <- available_base_variables(fit_obj, variables)
  
  if (length(have) == 0) return(tibble())
  
  tryCatch(
    posterior::as_draws_df(
      fit_obj$draws(variables = have, inc_warmup = FALSE)
    ) %>%
      as_tibble(),
    error = function(e) tibble()
  )
}

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

empty_gen_draws <- function() {
  tibble(
    model = character(),
    strain = character(),
    mean_gen = numeric(),
    sd_gen = numeric()
  )
}

# Fit
fit_registry <- list(
  read_fit(
    file.path(OLD_DIR, "fit_park_lognormal_fast_sgtf.rds"),
    "Lognormal",
    "SGTF"
  ),
  read_fit(
    file.path(OLD_DIR, "fit_park_lognormal_fast_nsgtf.rds"),
    "Lognormal",
    "non-SGTF"
  ),
  read_fit(
    file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds"),
    "Mechanistic",
    "SGTF"
  ),
  read_fit(
    file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds"),
    "Mechanistic",
    "non-SGTF"
  ),
  read_fit(
    find_censored_fit("sgtf"),
    "Censored WH-informed",
    "SGTF"
  ),
  read_fit(
    find_censored_fit("nsgtf"),
    "Censored WH-informed",
    "non-SGTF"
  )
)

print(
  tibble(
    model = vapply(fit_registry, `[[`, character(1), "model"),
    strain = vapply(fit_registry, `[[`, character(1), "strain"),
    path = vapply(fit_registry, `[[`, character(1), "path")
  )
)

# Extract generation-time draws
gen_draws <- purrr::map_dfr(fit_registry, function(x) {
  d <- try_draws_df(x$fit, c("mean_gen", "sd_gen"))
  
  if (nrow(d) == 0 || !all(c("mean_gen", "sd_gen") %in% names(d))) {
    return(empty_gen_draws())
  }
  
  d %>%
    transmute(
      model = x$model,
      strain = x$strain,
      mean_gen = .data$mean_gen,
      sd_gen = .data$sd_gen
    )
})

gen_draws <- gen_draws %>%
  mutate(
    model = factor(.data$model, levels = nice_model_levels),
    strain = factor(.data$strain, levels = nice_strain_levels)
  )

fig3_dat <- gen_draws %>%
  pivot_longer(
    c("mean_gen", "sd_gen"),
    names_to = "quantity",
    values_to = "value"
  ) %>%
  mutate(
    quantity = recode(
      .data$quantity,
      mean_gen = "Mean generation time",
      sd_gen = "SD generation time"
    )
  )

# Plot
p_fig3 <- ggplot(fig3_dat, aes(x = model, y = value, fill = strain)) +
  geom_violin(
    position = position_dodge(width = 0.85),
    trim = TRUE,
    linewidth = 0.25,
    alpha = 0.65
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.16,
    outlier.shape = NA,
    linewidth = 0.25,
    alpha = 0.85
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.85),
    size = 1.4,
    colour = "white"
  ) +
  facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = model_plot_labels) +
  scale_fill_discrete(labels = strain_plot_labels) +
  labs(x = NULL, y = "Days") +
  theme_paper(10) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

save_both(
  p_fig3,
  "figure3_generation_mean_sd_lognormal_vs_wh_vs_censored",
  width = 8.8,
  height = 3.9
)

write.csv(
  gen_draws,
  file.path(FIG_DIR, "figure3_generation_mean_sd_draws_lognormal_vs_wh_vs_censored.csv"),
  row.names = FALSE
)