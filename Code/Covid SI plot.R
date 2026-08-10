# Serial interval script
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble)
})

# Paths
OLD_DIR <- path.expand(
  "~/Within-host-time-delay-framework/Results/Covid/Noncens"
)

NEW_DIR <- "~/Within-host-time-delay-framework/Results/Covid/Cens"

SERIAL_FILE <- "~/Within-host-time-delay-framework/Data/serial-netherlands.xlsx"

FIG_DIR <- file.path(NEW_DIR, "figure7_style_three_model_overlay")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

LOGNORM_SGTF  <- file.path(OLD_DIR, "fit_park_lognormal_fast_sgtf.rds")
LOGNORM_NSGTF <- file.path(OLD_DIR, "fit_park_lognormal_fast_nsgtf.rds")

WH_SGTF  <- file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_sgtf.rds")
WH_NSGTF <- file.path(OLD_DIR, "fit_mechanistic_proportional_incubation_updatedV_obspair_nsgtf.rds")

# Figure labels and colours
model_levels <- c(
  "Lognormal",
  "WH-informed",
  "Censored WH-informed"
)

model_plot_labels <- c(
  "Lognormal" = "Lognormal",
  "WH-informed" = "WH-informed",
  "Censored WH-informed" = "Censored WH-informed"
)

model_cols <- c(
  "Lognormal" = "#F8766D",
  "WH-informed" = "#00BFC4",
  "Censored WH-informed" = "#619CFF"
)

strain_plot_labels <- c(
  "SGTF" = "Omicron",
  "non-SGTF" = "Delta"
)

strain_panel_levels <- c("Omicron", "Delta")

si_grid_lognormal_default <- seq(-10, 20, by = 0.2)
si_grid_wh_default        <- seq(-10, 20, by = 0.25)

# Find censored model fits
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

CENS_SGTF  <- find_censored_fit("sgtf")
CENS_NSGTF <- find_censored_fit("nsgtf")

# Helpers
read_fit <- function(path, model, strain) {
  if (!file.exists(path)) stop("Missing fit file: ", path)
  
  list(
    fit = readRDS(path),
    model = model,
    strain = strain,
    path = path
  )
}

fit_variables <- function(fit_obj) {
  posterior::variables(fit_obj$draws(inc_warmup = FALSE))
}

variable_regex <- function(variable) {
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", variable)
  paste0("^", escaped, "(\\[|$)")
}

has_base_variable <- function(fit_obj, variable) {
  vars <- fit_variables(fit_obj)
  any(stringr::str_detect(vars, variable_regex(variable)))
}

try_draws_matrix_one <- function(fit_obj, variable) {
  if (!has_base_variable(fit_obj, variable)) return(NULL)
  
  posterior::as_draws_matrix(
    fit_obj$draws(variables = variable, inc_warmup = FALSE)
  )
}

row_log_sum_exp <- function(x) {
  apply(x, 1, function(v) {
    m <- max(v, na.rm = TRUE)
    if (!is.finite(m)) return(-Inf)
    m + log(sum(exp(v - m), na.rm = TRUE))
  })
}

expand_serial_values <- function(serial_df) {
  rep(serial_df$serial, times = serial_df$n)
}

infer_serial_grid <- function(n_cols) {
  candidates <- list(
    si_grid_lognormal_default,
    si_grid_wh_default,
    seq(-10, 20, by = 0.1),
    seq(-5, 15, by = 0.25),
    seq(-5, 15, by = 0.2),
    seq(-5, 15, by = 0.1)
  )
  
  lens <- vapply(candidates, length, integer(1))
  idx <- which(lens == n_cols)
  
  if (length(idx) == 1) return(candidates[[idx]])
  
  stop("Could not infer serial grid from matrix with ", n_cols, " columns.")
}

theme_paper_overlay <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      axis.title = element_text(face = "plain"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90")
    )
}

make_weighted_serial <- function(df, strain_value, household_value = "within") {
  df %>%
    filter(
      .data$strain == strain_value,
      .data$household == household_value
    ) %>%
    mutate(serial = as.numeric(.data$serial)) %>%
    group_by(.data$serial) %>%
    dplyr::summarize(n = sum(.data$n), .groups = "drop") %>%
    filter(.data$n > 0) %>%
    arrange(.data$serial)
}

# Observed serial data
serialdata <- readxl::read_xlsx(SERIAL_FILE)

obs_serial <- bind_rows(
  make_weighted_serial(serialdata, "SGTF", "within") %>%
    mutate(strain = "SGTF"),
  make_weighted_serial(serialdata, "non-SGTF", "within") %>%
    mutate(strain = "non-SGTF")
) %>%
  group_by(.data$strain) %>%
  mutate(freq = .data$n / sum(.data$n)) %>%
  ungroup() %>%
  mutate(
    serial = as.numeric(.data$serial),
    strain_label = recode(.data$strain, !!!strain_plot_labels),
    strain_label = factor(.data$strain_label, levels = strain_panel_levels)
  )

# Extract fitted serial density
extract_serial_density_direct <- function(fit_obj, model, strain) {
  candidates <- c(
    "serial_prob_grid",
    "serial_density_norm_grid"
  )
  
  mat <- NULL
  used_name <- NULL
  
  for (nm in candidates) {
    tmp <- try_draws_matrix_one(fit_obj, nm)
    if (!is.null(tmp)) {
      mat <- tmp
      used_name <- nm
      break
    }
  }
  
  if (is.null(mat)) {
    logu <- try_draws_matrix_one(fit_obj, "log_unnorm_si")
    logn <- try_draws_matrix_one(fit_obj, "log_norm")
    
    if (!is.null(logu) && !is.null(logn)) {
      mat <- exp(sweep(logu, 1, as.numeric(logn[, 1]), FUN = "-"))
      used_name <- "log_unnorm_si/log_norm"
    }
  }
  
  if (is.null(mat)) return(NULL)
  
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  
  si_grid <- infer_serial_grid(ncol(mat))
  
  qs <- apply(
    mat,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  tibble(
    model = model,
    strain = strain,
    serial = as.numeric(si_grid),
    q025 = as.numeric(qs[1, ]),
    q50  = as.numeric(qs[2, ]),
    q975 = as.numeric(qs[3, ]),
    source = used_name
  )
}

extract_serial_density_from_loglik <- function(fit_obj, model, strain, obs_serial) {
  ll <- try_draws_matrix_one(fit_obj, "log_lik")
  
  if (is.null(ll)) {
    stop("No serial-density variable and no log_lik saved for ", model, " / ", strain)
  }
  
  obs_df <- obs_serial %>%
    filter(.data$strain == !!strain) %>%
    arrange(.data$serial)
  
  serial_expanded <- expand_serial_values(obs_df)
  
  if (ncol(ll) != length(serial_expanded)) {
    stop(
      "log_lik columns (", ncol(ll),
      ") do not match expanded observed serial length (",
      length(serial_expanded), ") for ", model, " / ", strain
    )
  }
  
  integer_support <- sort(unique(as.numeric(obs_df$serial)))
  first_col <- match(integer_support, serial_expanded)
  present <- !is.na(first_col)
  
  logp <- ll[, first_col[present], drop = FALSE]
  lse <- row_log_sum_exp(logp)
  prob <- exp(logp - lse)
  
  qs <- apply(
    prob,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  )
  
  tibble(
    model = model,
    strain = strain,
    serial = integer_support[present],
    q025 = as.numeric(qs[1, ]),
    q50  = as.numeric(qs[2, ]),
    q975 = as.numeric(qs[3, ]),
    source = "log_lik fallback"
  )
}

extract_serial_density <- function(fit_obj, model, strain, obs_serial) {
  out <- extract_serial_density_direct(fit_obj, model, strain)
  if (!is.null(out)) return(out)
  
  extract_serial_density_from_loglik(fit_obj, model, strain, obs_serial)
}

# Load fits
fit_registry <- list(
  read_fit(LOGNORM_SGTF,  "Lognormal",             "SGTF"),
  read_fit(WH_SGTF,       "WH-informed",           "SGTF"),
  read_fit(CENS_SGTF,     "Censored WH-informed",  "SGTF"),
  
  read_fit(LOGNORM_NSGTF, "Lognormal",             "non-SGTF"),
  read_fit(WH_NSGTF,      "WH-informed",           "non-SGTF"),
  read_fit(CENS_NSGTF,    "Censored WH-informed",  "non-SGTF")
)

print(tibble(
  model = vapply(fit_registry, `[[`, character(1), "model"),
  strain = vapply(fit_registry, `[[`, character(1), "strain"),
  path = vapply(fit_registry, `[[`, character(1), "path")
))

serial_density <- purrr::map_dfr(fit_registry, function(x) {
  extract_serial_density(x$fit, x$model, x$strain, obs_serial)
}) %>%
  mutate(
    strain_label = recode(.data$strain, !!!strain_plot_labels),
    strain_label = factor(.data$strain_label, levels = strain_panel_levels),
    model = factor(.data$model, levels = model_levels)
  )

# Plot
p <- ggplot() +
  geom_col(
    data = obs_serial,
    aes(x = serial, y = freq),
    fill = "grey70",
    colour = NA,
    width = 0.9,
    alpha = 0.75
  ) +
  geom_ribbon(
    data = serial_density,
    aes(
      x = serial,
      ymin = q025,
      ymax = q975,
      fill = model,
      group = interaction(model, strain_label)
    ),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(
    data = serial_density,
    aes(
      x = serial,
      y = q50,
      colour = model,
      group = interaction(model, strain_label)
    ),
    linewidth = 0.95
  ) +
  facet_wrap(~ strain_label, nrow = 1) +
  scale_colour_manual(
    values = model_cols,
    labels = model_plot_labels,
    breaks = model_levels
  ) +
  scale_fill_manual(
    values = model_cols,
    labels = model_plot_labels,
    breaks = model_levels
  ) +
  labs(
    x = "Observed serial interval, days",
    y = "Posterior predictive frequency"
  ) +
  theme_paper_overlay(13)

ggsave(
  file.path(FIG_DIR, "figure7_style_serial_fit_overlay_three_models_omicron_delta.png"),
  p,
  width = 10,
  height = 4.8,
  dpi = 450
)

ggsave(
  file.path(FIG_DIR, "figure7_style_serial_fit_overlay_three_models_omicron_delta.pdf"),
  p,
  width = 10,
  height = 4.8,
  device = cairo_pdf
)

write.csv(
  serial_density,
  file.path(FIG_DIR, "figure7_style_serial_fit_overlay_three_models_density.csv"),
  row.names = FALSE
)

print(p)
