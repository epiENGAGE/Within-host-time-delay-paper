# Plot Simulation 4: joint recovery of epidemiological delays
# -----------------------------------------------------------
# Main manuscript-style figure for Simulation 4.
#
# Panel A:
#   Full distribution recovery for latent period, incubation period,
#   generation time, and serial interval.
#   Dashed line = exact data-generating distribution
#   Solid line = median fitted distribution across simulation replicates
#   Ribbon     = 5th-95th percentile across fitted distributions
#
# Panel B:
#   Recovery of the mean and SD of all four delays.
#   Violin/boxplot = posterior median across simulation replicates
#   Cross          = known data-generating value
#
# The styling deliberately follows the Simulation 1/2 figures:
# theme_classic(), bold unboxed facet strips, light horizontal grid lines,
# restrained black/grey geometry, and patchwork panel labels.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
SIM_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 4"
SIM_DIR <- "C:/Users/nj7786/Downloads/WHsim4/Results/Simulations/Simulation 4"

FIG_DIR <- Sys.getenv("FIG_DIR", unset = SIM_DIR)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

delay_file <- file.path(
  SIM_DIR,
  "joint_delay_recovery_by_replicate.csv"
)

density_file <- file.path(
  SIM_DIR,
  "fitted_delay_densities_by_replicate.csv"
)

truth_file <- file.path(
  SIM_DIR,
  "true_joint_delay_targets.csv"
)

summary_file <- file.path(
  SIM_DIR,
  "paper_joint_delay_recovery_summary.csv"
)

required_files <- c(
  delay_file,
  density_file,
  truth_file,
  summary_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing Simulation 4 output file(s):\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
delay_results <- readr::read_csv(
  delay_file,
  show_col_types = FALSE
)

fitted_densities <- readr::read_csv(
  density_file,
  show_col_types = FALSE
)

truth <- readr::read_csv(
  truth_file,
  show_col_types = FALSE
)

recovery_summary <- readr::read_csv(
  summary_file,
  show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# Common labels and ordering
# -----------------------------------------------------------------------------
delay_levels <- c(
  "Latent period",
  "Incubation period",
  "Generation time",
  "Serial interval"
)

delay_short_levels <- c(
  "Latent",
  "Incubation",
  "Generation",
  "Serial"
)

delay_short_lookup <- c(
  "Latent period" = "Latent",
  "Incubation period" = "Incubation",
  "Generation time" = "Generation",
  "Serial interval" = "Serial"
)

# -----------------------------------------------------------------------------
# Theme -- matched to Simulation 1/2
# -----------------------------------------------------------------------------
theme_sim <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size + 1
      ),
      axis.title = ggplot2::element_text(face = "plain"),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.2,
        colour = "grey90"
      )
    )
}

# -----------------------------------------------------------------------------
# Reconstruct the true delay densities
#
# The Simulation 4 fitting script writes the true scalar targets but does not
# write the true density curves separately. The fitted-density output contains
# the same numerical support used for each delay, so the true curve can be
# reconstructed here from the replicate fits only if a true-density file was
# saved. To keep this plotting script standalone, first look for the optional
# file. If it is absent, use the true curve written by the Simulation 4 script
# if you add the small block described in the message below.
# -----------------------------------------------------------------------------
true_density_file <- file.path(
  SIM_DIR,
  "true_joint_delay_densities.csv"
)

if (!file.exists(true_density_file)) {
  stop(
    paste0(
      "\nThe plotting script needs the exact DGM density curves for Panel A.\n",
      "Please add this ONE-TIME export to the Simulation 4 script immediately ",
      "after TRUE_DENSITIES is created:\n\n",
      "  readr::write_csv(\n",
      "    TRUE_DENSITIES,\n",
      "    file.path(OUT_DIR, \"true_joint_delay_densities.csv\")\n",
      "  )\n\n",
      "Then either rerun only the setup section through TRUE_DENSITIES, or run ",
      "that write_csv line in the existing Simulation 4 R session. ",
      "You do NOT need to rerun Stan.\n"
    )
  )
}

true_densities <- readr::read_csv(
  true_density_file,
  show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# Panel A data: full distribution recovery
# -----------------------------------------------------------------------------
density_summary <- fitted_densities %>%
  dplyr::group_by(
    delay_code,
    delay,
    time
  ) %>%
  dplyr::summarise(
    fitted_q05 = stats::quantile(
      density,
      0.05,
      na.rm = TRUE
    ),
    fitted_median = stats::median(
      density,
      na.rm = TRUE
    ),
    fitted_q95 = stats::quantile(
      density,
      0.95,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    true_densities %>%
      dplyr::select(
        delay_code,
        delay,
        time,
        true_density
      ),
    by = c(
      "delay_code",
      "delay",
      "time"
    )
  ) %>%
  dplyr::mutate(
    delay = factor(
      delay,
      levels = delay_levels
    )
  )

if (any(!is.finite(density_summary$true_density))) {
  stop(
    "Some true density values are missing after joining fitted and true curves."
  )
}

# -----------------------------------------------------------------------------
# Panel A: joint recovery of the four full delay distributions
# -----------------------------------------------------------------------------
p_density <- ggplot2::ggplot(
  density_summary,
  ggplot2::aes(x = time)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = fitted_q05,
      ymax = fitted_q95
    ),
    fill = "grey75",
    alpha = 0.55,
    linewidth = 0
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = fitted_median),
    linewidth = 0.75,
    colour = "black"
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = true_density),
    linetype = "dashed",
    linewidth = 0.65,
    colour = "black"
  ) +
  ggplot2::facet_wrap(
    ~ delay,
    scales = "free",
    nrow = 1
  ) +
  ggplot2::labs(
    x = "Time, days",
    y = "Density"
  ) +
  theme_sim(11) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x = ggplot2::element_blank()
  )

# -----------------------------------------------------------------------------
# Panel B data: recovery of mean and SD
# -----------------------------------------------------------------------------
summary_plot <- delay_results %>%
  dplyr::mutate(
    delay = factor(
      delay,
      levels = delay_levels
    ),
    delay_short = dplyr::recode(
      as.character(delay),
      !!!delay_short_lookup
    ),
    delay_short = factor(
      delay_short,
      levels = delay_short_levels
    ),
    moment = factor(
      moment,
      levels = c("Mean", "SD")
    ),
    moment_label = dplyr::recode(
      as.character(moment),
      "Mean" = "Delay mean",
      "SD" = "Delay SD"
    ),
    moment_label = factor(
      moment_label,
      levels = c(
        "Delay mean",
        "Delay SD"
      )
    )
  )

truth_plot <- truth %>%
  dplyr::mutate(
    delay = factor(
      delay,
      levels = delay_levels
    ),
    delay_short = dplyr::recode(
      as.character(delay),
      !!!delay_short_lookup
    ),
    delay_short = factor(
      delay_short,
      levels = delay_short_levels
    ),
    moment = factor(
      moment,
      levels = c("Mean", "SD")
    ),
    moment_label = dplyr::recode(
      as.character(moment),
      "Mean" = "Delay mean",
      "SD" = "Delay SD"
    ),
    moment_label = factor(
      moment_label,
      levels = c(
        "Delay mean",
        "Delay SD"
      )
    )
  )

# -----------------------------------------------------------------------------
# Panel B: posterior-median recovery of the delay summaries
#
# This deliberately mirrors Simulation 1:
# violin + narrow boxplot + jittered replicate points + truth marker.
# -----------------------------------------------------------------------------
p_summary <- ggplot2::ggplot(
  summary_plot,
  ggplot2::aes(
    x = delay_short,
    y = median
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.65,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.35,
    alpha = 0.90
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0,
      seed = 44004
    ),
    alpha = 0.18,
    size = 0.8
  ) +
  ggplot2::stat_summary(
    fun = stats::median,
    geom = "point",
    shape = 21,
    size = 1.8,
    fill = "white",
    colour = "black",
    stroke = 0.3
  ) +
  ggplot2::geom_point(
    data = truth_plot,
    ggplot2::aes(
      x = delay_short,
      y = true_value
    ),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.0,
    size = 2.7
  ) +
  ggplot2::facet_wrap(
    ~ moment_label,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::labs(
    x = "Epidemiological delay",
    y = "Posterior median estimate, days"
  ) +
  theme_sim(11) +
  ggplot2::theme(
    legend.position = "none"
  )

# -----------------------------------------------------------------------------
# Main combined figure
# -----------------------------------------------------------------------------
p_main <- (p_density / p_summary) +
  patchwork::plot_layout(
    heights = c(1.05, 1.0)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 16
    ),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_main)

ggplot2::ggsave(
  file.path(
    FIG_DIR,
    "figure_sim4_joint_delay_recovery.png"
  ),
  p_main,
  width = 11.5,
  height = 7.4,
  dpi = 450
)

ggplot2::ggsave(
  file.path(
    FIG_DIR,
    "figure_sim4_joint_delay_recovery.pdf"
  ),
  p_main,
  width = 11.5,
  height = 7.4,
  device = cairo_pdf
)

# -----------------------------------------------------------------------------
# Optional secondary figure: probabilistic performance
#
# Kept separate from the main figure so the manuscript message remains:
# "the model jointly recovers all four epidemiological delays."
# -----------------------------------------------------------------------------
crps_plot_dat <- delay_results %>%
  dplyr::mutate(
    delay_short = dplyr::recode(
      delay,
      !!!delay_short_lookup
    ),
    delay_short = factor(
      delay_short,
      levels = delay_short_levels
    ),
    moment_label = dplyr::recode(
      moment,
      "Mean" = "Delay mean",
      "SD" = "Delay SD"
    ),
    moment_label = factor(
      moment_label,
      levels = c(
        "Delay mean",
        "Delay SD"
      )
    )
  )

p_crps <- ggplot2::ggplot(
  crps_plot_dat,
  ggplot2::aes(
    x = delay_short,
    y = crps
  )
) +
  ggplot2::geom_boxplot(
    outlier.alpha = 0.20,
    width = 0.55,
    linewidth = 0.35
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.12,
      height = 0,
      seed = 44005
    ),
    alpha = 0.25,
    size = 0.8
  ) +
  ggplot2::facet_wrap(
    ~ moment_label,
    scales = "free_y",
    nrow = 1
  ) +
  ggplot2::labs(
    x = "Epidemiological delay",
    y = "CRPS, days"
  ) +
  theme_sim(11)

ggplot2::ggsave(
  file.path(
    FIG_DIR,
    "figure_sim4_joint_delay_crps_secondary.png"
  ),
  p_crps,
  width = 9.5,
  height = 4.5,
  dpi = 450
)

ggplot2::ggsave(
  file.path(
    FIG_DIR,
    "figure_sim4_joint_delay_crps_secondary.pdf"
  ),
  p_crps,
  width = 9.5,
  height = 4.5,
  device = cairo_pdf
)

# -----------------------------------------------------------------------------
# Paper-ready quantitative table
# -----------------------------------------------------------------------------
paper_table <- recovery_summary %>%
  dplyr::mutate(
    delay = factor(
      delay,
      levels = delay_levels
    ),
    moment = factor(
      moment,
      levels = c("Mean", "SD")
    )
  ) %>%
  dplyr::arrange(
    delay,
    moment
  ) %>%
  dplyr::select(
    delay,
    moment,
    n_successful_fits,
    true_value,
    median_bias,
    median_abs_error,
    rmse,
    mean_crps,
    median_interval_width,
    coverage_90
  )

readr::write_csv(
  paper_table,
  file.path(
    FIG_DIR,
    "table_sim4_joint_delay_recovery_for_paper.csv"
  )
)

cat(
  "\nSimulation 4 figure written to:\n  ",
  file.path(
    FIG_DIR,
    "figure_sim4_joint_delay_recovery.png"
  ),
  "\n\n",
  sep = ""
)

print(
  paper_table,
  n = Inf
)
