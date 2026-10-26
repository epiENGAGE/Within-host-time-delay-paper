functions {
  real lognormal_lpdf_from_x(real x, real mu, real sigma) {
    return lognormal_lpdf(x | mu, sigma);
  }
}

data {
  int<lower=1> N_norm;
  vector[N_norm] si_norm_grid;
  real<lower=0> si_norm_step;

  int<lower=1> N_si;
  array[N_si] int<lower=1, upper=N_norm> obs_index;
  array[N_si] int<lower=0> counts;

  // Expanded pointwise observations, one entry per transmission pair.
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N_norm> obs_index_obs;

  // Precomputed incubation-difference convolution grid.
  // delta = IP2 - IP1, so SI = GT + delta and GT = SI - delta.
  int<lower=1> N_delta;
  vector[N_delta] delta_grid;
  real<lower=0> delta_step;
  vector<lower=0>[N_delta] C_delta_density;

  real prior_logmean_gen_mean;
  real<lower=0> prior_logmean_gen_sd;
  real prior_logsd_gen_mean;
  real<lower=0> prior_logsd_gen_sd;
}

parameters {
  real logmean_gen;
  real<lower=0> logsd_gen;
}

transformed parameters {
  vector[N_norm] log_unnorm_norm;
  real log_norm;
  vector[N_si] log_p_si;

  for (k in 1:N_norm) {
    vector[N_delta] term_lp;

    for (dd in 1:N_delta) {
      real gen = si_norm_grid[k] - delta_grid[dd];

      if (gen > 0 && C_delta_density[dd] > 0) {
        term_lp[dd] = log(C_delta_density[dd])
                      + log(delta_step)
                      + lognormal_lpdf_from_x(gen, logmean_gen, logsd_gen);
      } else {
        term_lp[dd] = negative_infinity();
      }
    }

    log_unnorm_norm[k] = log_sum_exp(term_lp);
  }

  log_norm = log_sum_exp(log_unnorm_norm) + log(si_norm_step);

  for (n in 1:N_si) {
    log_p_si[n] = log_unnorm_norm[obs_index[n]] - log_norm;
  }
}

model {
  logmean_gen ~ normal(prior_logmean_gen_mean, prior_logmean_gen_sd);
  logsd_gen ~ normal(prior_logsd_gen_mean, prior_logsd_gen_sd);

  for (n in 1:N_si) {
    target += counts[n] * log_p_si[n];
  }
}

generated quantities {
  real mean_gen;
  real sd_gen;
  vector[N_obs] log_lik;
  vector[N_si] log_lik_grouped;
  vector[N_norm] serial_density_norm_grid;

  mean_gen = exp(logmean_gen + 0.5 * square(logsd_gen));

  sd_gen = sqrt(
    (exp(square(logsd_gen)) - 1)
    * exp(2 * logmean_gen + square(logsd_gen))
  );

  for (n in 1:N_obs) {
    log_lik[n] = log_unnorm_norm[obs_index_obs[n]] - log_norm;
  }

  for (n in 1:N_si) {
    log_lik_grouped[n] = counts[n] * log_p_si[n];
  }

  for (k in 1:N_norm) {
    serial_density_norm_grid[k] = exp(log_unnorm_norm[k] - log_norm);
  }
}
