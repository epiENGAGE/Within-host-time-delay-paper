functions {
  real lognormal_lpdf_from_logx(real log_x, real mu, real sigma) {
    return normal_lpdf(log_x | mu, sigma) - log_x;
  }
}

data {
  int<lower=1> N_norm;
  vector[N_norm] si_norm_grid;
  real<lower=0> si_norm_step;

  int<lower=1> N_si;
  array[N_si] int<lower=1, upper=N_norm> obs_index;
  array[N_si] int<lower=0> counts;

  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N_norm> obs_index_obs;

  int<lower=1> N_terms;
  array[N_terms] int<lower=1, upper=N_norm> s_id;

  vector[N_terms] const_lp;
  vector[N_terms] log_inc1;
  vector[N_terms] log_gen;

  real logmean_inc;
  real<lower=0> logsd_inc;
  real<lower=-1, upper=1> rho;

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

    for (k in 1:N_norm) {
      log_unnorm_norm[k] = negative_infinity();
    }

    // Accumulate log_sum_exp by serial-grid point.
    for (m in 1:N_terms) {
      int k = s_id[m];
      log_unnorm_norm[k] = log_sum_exp(log_unnorm_norm[k], term_lp[m]);
    }
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

  sd_gen =
    sqrt(
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
