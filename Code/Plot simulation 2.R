# Plot simulation 3
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# Paths
SIM_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 2"

FIG_DIR <- SIM_DIR

fit_file <- file.path(SIM_DIR, "simulation_fit_level_results.csv")
summary_file <- file.path(SIM_DIR, "paper_recovery_summary_by_regime.csv")
truth_file <- file.path(SIM_DIR, "true_regime_targets.csv")

fit_results <- read_csv(fit_file, show_col_types = FALSE)
recovery_summary <- read_csv(summary_file, show_col_types = FALSE)
truth <- read_csv(truth_file, show_col_types = FALSE)

# Clean data
regime_lookup <- truth %>%
  arrange(target_mean_gen) %>%
  transmute(
    regime,
    target_mean_gen,
    regime_label = paste0(round(target_mean_gen), "-day regime"),
    regime_label = factor(regime_label, levels = paste0(round(target_mean_gen), "-day regime")),
    true_mean_gen,
    true_sd_gen
  )

plot_results <- fit_results %>%
  left_join(regime_lookup %>% select(regime, regime_label), by = "regime") %>%
  mutate(
    regime_label = factor(regime_label, levels = levels(regime_lookup$regime_label)),
    mean_gen_covered = as.logical(mean_gen_covered),
    sd_gen_covered = as.logical(sd_gen_covered),
    convergence_ok = max_rhat < 1.05 & n_divergent == 0
  )

estimate_long <- plot_results %>%
  select(
    regime, regime_label, target_mean_gen, rep_id,
    mean_gen_median, mean_gen_q05, mean_gen_q95, mean_gen_covered,
    sd_gen_median, sd_gen_q05, sd_gen_q95, sd_gen_covered,
    true_mean_gen, true_sd_gen
  ) %>%
  pivot_longer(
    cols = c(
      mean_gen_median, mean_gen_q05, mean_gen_q95, mean_gen_covered,
      sd_gen_median, sd_gen_q05, sd_gen_q95, sd_gen_covered
    ),
    names_to = c("quantity", ".value"),
    names_pattern = "(mean_gen|sd_gen)_(median|q05|q95|covered)"
  ) %>%
  mutate(
    quantity = recode(
      quantity,
      mean_gen = "Generation-time mean",
      sd_gen = "Generation-time SD"
    ),
    quantity = factor(quantity, levels = c("Generation-time mean", "Generation-time SD")),
    true_value = ifelse(quantity == "Generation-time mean", true_mean_gen, true_sd_gen),
    error = median - true_value,
    abs_error = abs(error),
    abs_percent_error = 100 * abs_error / true_value,
    interval_width = q95 - q05,
    within_0.5_days = abs_error <= 0.5,
    within_1_day = abs_error <= 1.0,
    covered = as.logical(covered)
  )

# Offsets prevent points/intervals from sitting directly on top of one another.
offsets <- estimate_long %>%
  distinct(quantity, regime_label, rep_id) %>%
  group_by(quantity, regime_label) %>%
  mutate(
    rep_rank = row_number(),
    n_rep_regime = n(),
    offset_index = rep_rank - (n_rep_regime + 1) / 2
  ) %>%
  ungroup() %>%
  mutate(offset = offset_index * 0.018)

estimate_long <- estimate_long %>%
  left_join(
    offsets %>% select(quantity, regime_label, rep_id, offset),
    by = c("quantity", "regime_label", "rep_id")
  ) %>%
  mutate(
    true_value_plot = true_value * (1 + offset)
  )

median_error <- estimate_long %>%
  group_by(regime_label, target_mean_gen, quantity) %>%
  summarise(
    median_abs_percent_error = median(abs_percent_error, na.rm = TRUE),
    .groups = "drop"
  )

# Theme
theme_sim <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(colour = "black"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey94"),
      plot.title = element_blank()
    )
}

# Panel A: recovery against truth
p_identity <- ggplot(
  estimate_long,
  aes(x = true_value_plot, y = median, ymin = q05, ymax = q95)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_errorbar(
    width = 0,
    alpha = 0.45,
    linewidth = 0.35
  ) +
  geom_point(
    alpha = 0.80,
    size = 1.8
  ) +
  facet_wrap(~ quantity, scales = "free", nrow = 1) +
  labs(
    x = "True data-generating value, days",
    y = "Posterior estimate, days"
  ) +
  theme_sim(11)

# Panel B: percentage error
p_percent_error <- ggplot(
  estimate_long,
  aes(x = regime_label, y = abs_percent_error)
) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 52002),
    alpha = 0.70,
    size = 1.9
  ) +
  geom_crossbar(
    data = median_error,
    aes(
      x = regime_label,
      y = median_abs_percent_error,
      ymin = median_abs_percent_error,
      ymax = median_abs_percent_error
    ),
    inherit.aes = FALSE,
    width = 0.45,
    linewidth = 0.50
  ) +
  facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
  labs(
    x = "Data-generating regime",
    y = "Absolute percentage error, %"
  ) +
  theme_sim(11) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

# Combine and save
p_combined <- (p_identity / p_percent_error) +
  plot_layout(heights = c(1.1, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  file.path(FIG_DIR, "figure_sim_regime_recovery_identity.png"),
  p_combined,
  width = 10.5,
  height = 7.0,
  dpi = 450
)

ggsave(
  file.path(FIG_DIR, "figure_sim_regime_recovery_identity.pdf"),
  p_combined,
  width = 10.5,
  height = 7.0,
  device = cairo_pdf
)

print(p_combined)

# Quantitative summaries
quant_summary <- estimate_long %>%
  group_by(target_mean_gen, regime_label, quantity) %>%
  summarise(
    n_replicates = n(),
    true_value = first(true_value),
    
    median_posterior_median = median(median, na.rm = TRUE),
    q25_posterior_median = quantile(median, 0.25, na.rm = TRUE),
    q75_posterior_median = quantile(median, 0.75, na.rm = TRUE),
    
    median_bias = median(error, na.rm = TRUE),
    mean_bias = mean(error, na.rm = TRUE),
    
    median_abs_error = median(abs_error, na.rm = TRUE),
    q75_abs_error = quantile(abs_error, 0.75, na.rm = TRUE),
    q90_abs_error = quantile(abs_error, 0.90, na.rm = TRUE),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    
    median_abs_percent_error = median(abs_percent_error, na.rm = TRUE),
    q90_abs_percent_error = quantile(abs_percent_error, 0.90, na.rm = TRUE),
    
    proportion_within_0.5_days = mean(within_0.5_days, na.rm = TRUE),
    proportion_within_1_day = mean(within_1_day, na.rm = TRUE),
    median_interval_width = median(interval_width, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(quantity, target_mean_gen)

cat("\n--- Quantitative summary by regime ---\n\n")
print(quant_summary, n = Inf)

write_csv(
  quant_summary,
  file.path(FIG_DIR, "section_3_2_2_quantitative_summary.csv")
)

text_numbers <- quant_summary %>%
  mutate(
    quantity_short = recode(
      quantity,
      "Generation-time mean" = "mean",
      "Generation-time SD" = "sd"
    )
  ) %>%
  select(
    target_mean_gen,
    quantity_short,
    true_value,
    median_posterior_median,
    median_abs_error,
    median_abs_percent_error,
    q90_abs_error,
    q90_abs_percent_error,
    proportion_within_0.5_days,
    proportion_within_1_day,
    median_interval_width,
    coverage
  ) %>%
  arrange(quantity_short, target_mean_gen)

cat("\n--- Text numbers ---\n\n")
print(text_numbers, n = Inf)

write_csv(
  text_numbers,
  file.path(FIG_DIR, "section_3_2_2_text_numbers.csv")
)
