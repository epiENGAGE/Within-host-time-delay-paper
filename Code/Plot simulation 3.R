# Plot simulation 3
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Paths and data loading
OUT_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 3"

fit_summary_by_prior <- read_csv(
  file.path(OUT_DIR, "paper_fit_summary_by_prior.csv"),
  show_col_types = FALSE
)

param_results <- read_csv(
  file.path(OUT_DIR, "parameter_recovery_results.csv"),
  show_col_types = FALSE
)

# Parameters targeted by the prior-sensitivity analysis
PRIOR_TARGET_SET <- c(
  "log_r",
  "log_d",
  "m",
  "log_kappa",
  "log_sigma_z"
)

plot_heat_data <- fit_summary_by_prior %>%
  filter(prior_scenario != "manuscript_default") %>%
  mutate(
    accuracy = factor(
      accuracy,
      levels = c("accurate", "mildly_wrong", "badly_wrong"),
      labels = c("Accurate prior", "Mildly wrong prior", "Badly wrong prior")
    ),
    strength = factor(
      strength,
      levels = c("weak", "moderate", "strong"),
      labels = c("Weak", "Moderate", "Strong")
    )
  ) %>%
  select(
    accuracy,
    strength,
    median_abs_error_mean_gen,
    median_abs_error_sd_gen,
    median_kl_true_to_postmean,
    median_test_lpd_per_pair
  ) %>%
  pivot_longer(
    cols = c(
      median_abs_error_mean_gen,
      median_abs_error_sd_gen,
      median_kl_true_to_postmean,
      median_test_lpd_per_pair
    ),
    names_to = "metric_raw",
    values_to = "value"
  ) %>%
  group_by(accuracy, metric_raw) %>%
  mutate(
    weak_value = value[strength == "Weak"][1],
    improvement_vs_weak = case_when(
      metric_raw == "median_test_lpd_per_pair" ~ value - weak_value,
      TRUE ~ weak_value - value
    ),
    percent_improvement_vs_weak = 100 * improvement_vs_weak / abs(weak_value)
  ) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric_raw,
      median_abs_error_mean_gen = "GT mean absolute error",
      median_abs_error_sd_gen = "GT SD absolute error",
      median_kl_true_to_postmean = "KL distance to true SI",
      median_test_lpd_per_pair = "Held-out LPD per pair"
    ),
    metric = factor(
      metric,
      levels = c(
        "GT mean absolute error",
        "GT SD absolute error",
        "KL distance to true SI",
        "Held-out LPD per pair"
      )
    ),
    label = case_when(
      metric == "Held-out LPD per pair" ~ sprintf(
        "%.5f\n%+.3f%%",
        value,
        percent_improvement_vs_weak
      ),
      TRUE ~ sprintf(
        "%.4f\n%+.1f%%",
        value,
        percent_improvement_vs_weak
      )
    ),
    panel = "A",
    parameter = NA_character_,
    moved_toward_truth_rate = NA_real_
  ) %>%
  select(
    panel,
    accuracy,
    strength,
    metric,
    value,
    weak_value,
    improvement_vs_weak,
    percent_improvement_vs_weak,
    label,
    parameter,
    moved_toward_truth_rate
  )

# Panel B data: posterior correction rate for wrong priors only
plot_move_data <- param_results %>%
  filter(
    variable %in% PRIOR_TARGET_SET,
    prior_scenario != "manuscript_default",
    accuracy %in% c("mildly_wrong", "badly_wrong")
  ) %>%
  group_by(accuracy, strength, variable) %>%
  summarise(
    moved_toward_truth_rate = mean(posterior_moved_toward_truth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    accuracy = factor(
      accuracy,
      levels = c("mildly_wrong", "badly_wrong"),
      labels = c("Mildly wrong prior", "Badly wrong prior")
    ),
    strength = factor(
      strength,
      levels = c("weak", "moderate", "strong"),
      labels = c("Weak", "Moderate", "Strong")
    ),
    parameter = recode(
      variable,
      log_r = "r",
      log_d = "d",
      m = "m",
      log_kappa = "κ",
      log_sigma_z = "σ_z"
    ),
    parameter = factor(
      parameter,
      levels = c("r", "d", "m", "κ", "σ_z")
    ),
    panel = "B",
    metric = NA_character_,
    value = NA_real_,
    weak_value = NA_real_,
    improvement_vs_weak = NA_real_,
    percent_improvement_vs_weak = NA_real_,
    label = NA_character_
  ) %>%
  select(
    panel,
    accuracy,
    strength,
    metric,
    value,
    weak_value,
    improvement_vs_weak,
    percent_improvement_vs_weak,
    label,
    parameter,
    moved_toward_truth_rate
  )

# Save one CSV for moving to another RStudio instance
sim3_panel_AB_plot_data <- bind_rows(
  plot_heat_data,
  plot_move_data
)

write_csv(
  sim3_panel_AB_plot_data,
  file.path(OUT_DIR, "sim3_panel_AB_plot_data.csv"),
  na = ""
)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

OUT_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 3"
FIG_DIR <- OUT_DIR
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

sim3_panel_AB_plot_data <- read_csv(
  file.path(OUT_DIR, "sim3_panel_AB_plot_data.csv"),
  show_col_types = FALSE
)

plot_heat <- sim3_panel_AB_plot_data %>%
  filter(panel == "A")

plot_move <- sim3_panel_AB_plot_data %>%
  filter(panel == "B") %>%
  mutate(
    parameter_plotmath = recode(
      as.character(parameter),
      "r" = "r",
      "d" = "d",
      "m" = "m",
      "κ" = "kappa",
      "σ_z" = "sigma[z]"
    )
  )

# Panel A
p_heat <- ggplot(
  plot_heat,
  aes(x = strength, y = accuracy, fill = percent_improvement_vs_weak)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 4.0, lineheight = 0.9) +
  facet_wrap(~ metric, nrow = 2) +
  scale_fill_gradient2(
    low = "firebrick3",
    mid = "white",
    high = "steelblue4",
    midpoint = 0,
    name = "Improvement\nvs weak prior (%)"
  ) +
  labs(
    x = "Prior strength",
    y = "Prior accuracy"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.margin = margin(8, 10, 8, 8)
  )

# Panel B
p_move_clean <- ggplot(
  plot_move,
  aes(x = strength, y = moved_toward_truth_rate)
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5) +
  geom_col(width = 0.65, fill = "grey35") +
  facet_grid(
    accuracy ~ parameter_plotmath,
    labeller = labeller(parameter_plotmath = label_parsed)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Prior strength",
    y = "Posterior closer to truth than prior mean"
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

# Combined Panel A + Panel B figure
p_AB <- p_heat / p_move_clean +
  plot_layout(heights = c(1.15, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_AB)

ggsave(
  file.path(FIG_DIR, "figure_sim3_panel_A_B_prior_sensitivity_large_text.png"),
  p_AB,
  width = 14,
  height = 14,
  dpi = 400
)

ggsave(
  file.path(FIG_DIR, "figure_sim3_panel_A_B_prior_sensitivity_large_text.pdf"),
  p_AB,
  width = 14,
  height = 14,
  device = cairo_pdf
)