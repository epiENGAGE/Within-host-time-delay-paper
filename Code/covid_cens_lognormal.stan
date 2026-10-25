functions {
  real lognormal_lpdf_from_logx(real log_x, real mu, real sigma) {
    return normal_lpdf(log_x | mu, sigma) - log_x;
  }
}

data {
  // Reported integer serial-interval categories.
  int<lower=1> N_cat;
  array[N_cat] int<lower=0> counts;
  vector[N_cat] serial_cat;

  // Expanded observation indices for pointwise LOO/WAIC.
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N_cat> obs_cat_pt;

  // Fine latent serial-interval grid used before date coarsening.
  int<lower=2> N_si_latent;
  vector[N_si_latent] si_latent_grid;
  real<lower=0> dx_si;

  // Sparse quadrature terms for the Park-style SI convolution.
  int<lower=1> N_terms;
  array[N_terms] int<lower=1, upper=N_si_latent> s_id;
  vector[N_terms] const_lp;
  vector[N_terms] log_inc1;
  vector[N_terms] log_gen;

  // Observation kernel: probability weight of reported category k
  // conditional on a latent continuous SI grid point j. For nearest-day
  // recording of both onset dates this is triangular: max(0,1-|s-y|).
  matrix<lower=0>[N_cat, N_si_latent] coarse_w;

  // Parameters of the fixed forward incubation-period lognormal approximation.
  real logmean_inc;
  real<lower=0> logsd_inc;
  real<lower=-0.999, upper=0.999> rho;

  // Weakly informative priors for the lognormal generation interval.
  real prior_logmean_gen_mean;
  real<lower=0> prior_logmean_gen_sd;
  real<lower=0> prior_logsd_gen_mean;
  real<lower=0> prior_logsd_gen_sd;
}

parameters {
  real logmean_gen;
  real<lower=0.02, upper=2.5> logsd_gen;
}

transformed parameters {
  vector[N_si_latent] log_unnorm_latent;
  real log_norm_latent;
  vector[N_si_latent] latent_density;
  vector[N_cat] p_cat_raw;
  vector[N_cat] p_cat;

  {
    vector[N_terms] term_lp;
    real cond_sd = logsd_gen * sqrt(1 - square(rho));

    for (m in 1:N_terms) {
      real cond_mean =
        logmean_gen
        + logsd_gen * rho * (log_inc1[m] - logmean_inc) / logsd_inc;

      term_lp[m] =
        const_lp[m]
        + lognormal_lpdf_from_logx(log_gen[m], cond_mean, cond_sd);
    }

    log_unnorm_latent = rep_vector(negative_infinity(), N_si_latent);

    for (m in 1:N_terms) {
      int k = s_id[m];
      log_unnorm_latent[k] = log_sum_exp(log_unnorm_latent[k], term_lp[m]);
    }

    // Normalize the latent SI density over the numerical support.
    log_norm_latent = log_sum_exp(log_unnorm_latent) + log(dx_si);
    for (j in 1:N_si_latent) {
      latent_density[j] = exp(log_unnorm_latent[j] - log_norm_latent);
    }

    // Integrate the latent density through the day-level observation kernel.
    p_cat_raw = coarse_w * latent_density * dx_si;

    // Match the censored WH model: condition on the reported-category support.
    p_cat = p_cat_raw / fmax(sum(p_cat_raw), 1e-300);
  }
}

model {
  logmean_gen ~ normal(prior_logmean_gen_mean, prior_logmean_gen_sd);
  logsd_gen ~ normal(prior_logsd_gen_mean, prior_logsd_gen_sd);

  for (k in 1:N_cat) {
    if (counts[k] > 0) {
      target += counts[k] * log(fmax(p_cat[k], 1e-300));
    }
  }
}

generated quantities {
  real mean_gen;
  real sd_gen;
  vector[N_obs] log_lik;
  vector[N_cat] log_lik_grouped;
  vector[N_cat] serial_prob_category;
  vector[N_si_latent] serial_density_latent_grid;

  mean_gen = exp(logmean_gen + 0.5 * square(logsd_gen));
  sd_gen = sqrt(
    (exp(square(logsd_gen)) - 1)
    * exp(2 * logmean_gen + square(logsd_gen))
  );

  for (n in 1:N_obs) {
    log_lik[n] = log(fmax(p_cat[obs_cat_pt[n]], 1e-300));
  }

  for (k in 1:N_cat) {
    log_lik_grouped[k] = counts[k] * log(fmax(p_cat[k], 1e-300));
    serial_prob_category[k] = p_cat[k];
  }

  serial_density_latent_grid = latent_density;
}
