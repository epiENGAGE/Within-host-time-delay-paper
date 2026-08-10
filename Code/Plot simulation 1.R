# Plot sim 1
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# Paths
SIM_DIR <- "~/Within-host-time-delay-framework/Results/Simulations/Simulation 1"

FIG_DIR <- Sys.getenv("FIG_DIR", unset = SIM_DIR)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

fit_file <- file.path(SIM_DIR, "simulation_fit_level_results.csv")
summary_file <- file.path(SIM_DIR, "paper_recovery_summary_by_sample_size.csv")
truth_file <- file.path(SIM_DIR, "true_generation_time_targets.csv")

fit_results <- read_csv(fit_file, show_col_types = FALSE)
recovery_summary <- read_csv(summary_file, show_col_types = FALSE)
truth <- read_csv(truth_file, show_col_types = FALSE)

# Clean data
sample_levels <- sort(unique(fit_results$n_pairs))

plot_results <- fit_results %>%
  mutate(
    n_pairs_f = factor(n_pairs, levels = sample_levels),
    mean_gen_covered = as.logical(mean_gen_covered),
    sd_gen_covered = as.logical(sd_gen_covered),
    convergence_ok = max_rhat < 1.05 & n_divergent == 0
  )

estimate_long <- plot_results %>%
  select(
    n_pairs, n_pairs_f, rep_id,
    mean_gen_median, mean_gen_q05, mean_gen_q95, mean_gen_covered,
    sd_gen_median, sd_gen_q05, sd_gen_q95, sd_gen_covered
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
    quantity = factor(
      quantity,
      levels = c("Generation-time mean", "Generation-time SD")
    )
  )

truth_long <- tibble(
  quantity = factor(
    c("Generation-time mean", "Generation-time SD"),
    levels = c("Generation-time mean", "Generation-time SD")
  ),
  true_value = c(truth$true_mean_gen[1], truth$true_sd_gen[1])
)

error_long <- plot_results %>%
  select(n_pairs, n_pairs_f, rep_id, mean_gen_abs_error, sd_gen_abs_error) %>%
  pivot_longer(
    cols = c(mean_gen_abs_error, sd_gen_abs_error),
    names_to = "quantity",
    values_to = "absolute_error_days"
  ) %>%
  mutate(
    quantity = recode(
      quantity,
      mean_gen_abs_error = "Generation-time mean",
      sd_gen_abs_error = "Generation-time SD"
    ),
    quantity = factor(
      quantity,
      levels = c("Generation-time mean", "Generation-time SD")
    )
  )

# Theme
theme_sim <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1),
      axis.title = element_text(face = "plain"),
      plot.title = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90")
    )
}

# Row A: posterior-median recovery of target summaries
p_recovery <- ggplot(
  estimate_long,
  aes(x = n_pairs_f, y = median)
) +
  geom_hline(
    data = truth_long,
    aes(yintercept = true_value),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.65,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.35,
    alpha = 0.9
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 11),
    alpha = 0.18,
    size = 0.8
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.8,
    fill = "white",
    colour = "black",
    stroke = 0.3
  ) +
  facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = "Posterior median estimate, days"
  ) +
  theme_sim(11)

# Row B: absolute error distribution
p_error <- ggplot(
  error_long,
  aes(x = n_pairs_f, y = absolute_error_days)
) +
  geom_boxplot(
    outlier.alpha = 0.20,
    width = 0.55,
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 12),
    alpha = 0.25,
    size = 0.8
  ) +
  facet_wrap(~ quantity, scales = "free_y", nrow = 1) +
  labs(
    x = "Number of simulated transmission pairs",
    y = "Absolute error, days"
  ) +
  theme_sim(11) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

# Combine and save
p_combined <- (p_recovery / p_error) +
  plot_layout(heights = c(1.05, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  file.path(FIG_DIR, "figure_sim_small_sample_recovery.png"),
  p_combined,
  width = 10.5,
  height = 7.2,
  dpi = 450
)

ggsave(
  file.path(FIG_DIR, "figure_sim_small_sample_recovery.pdf"),
  p_combined,
  width = 10.5,
  height = 7.2,
  device = cairo_pdf
)

# Table of results
paper_table <- recovery_summary %>%
  transmute(
    n_pairs,
    n_successful_fits,
    median_abs_error_mean_gen,
    rmse_mean_gen,
    median_ci_width_mean_gen,
    coverage_mean_gen,
    median_abs_error_sd_gen,
    rmse_sd_gen,
    median_ci_width_sd_gen,
    coverage_sd_gen,
    convergence_ok_rate
  )

write_csv(
  paper_table,
  file.path(FIG_DIR, "table_sim_small_sample_recovery_for_paper.csv")
)

print(p_combined)

# More summaries
library(dplyr)
library(tidyr)
library(readr)

# Compact quantitative summary by sample size and quantity
quant_summary <- estimate_long %>%
  left_join(truth_long, by = "quantity") %>%
  mutate(
    error = median - true_value,
    abs_error = abs(error),
    interval_width = q95 - q05,
    within_0.5_days = abs_error <= 0.5,
    within_1_day = abs_error <= 1.0,
    covered = as.logical(covered)
  ) %>%
  group_by(n_pairs, quantity) %>%
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
    
    proportion_within_0.5_days = mean(within_0.5_days, na.rm = TRUE),
    proportion_within_1_day = mean(within_1_day, na.rm = TRUE),
    
    median_interval_width = median(interval_width, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    
    .groups = "drop"
  )

print(quant_summary, n = Inf)

# More useful values
text_numbers <- quant_summary %>%
  mutate(
    quantity_short = recode(
      quantity,
      "Generation-time mean" = "mean",
      "Generation-time SD" = "sd"
    )
  ) %>%
  select(
    n_pairs,
    quantity_short,
    true_value,
    median_posterior_median,
    median_abs_error,
    q90_abs_error,
    proportion_within_0.5_days,
    proportion_within_1_day,
    median_interval_width,
    coverage
  ) %>%
  arrange(quantity_short, n_pairs)

print(text_numbers, n = Inf)

write_csv(
  text_numbers,
  file.path(FIG_DIR, "section_3_2_1_text_numbers.csv")
)
