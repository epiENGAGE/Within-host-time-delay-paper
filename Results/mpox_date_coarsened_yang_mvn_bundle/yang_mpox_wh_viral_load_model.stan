functions {
  real log_viral_shape_one(real t, real log_r, real log_d, real m, real log_kappa) {
    real rr = exp(log_r);
    real dd = exp(log_d);
    real kappa = exp(log_kappa);
    real a = log(rr / dd) + kappa * (t - m);
    real log1pexp_a;

    if (a > 40) log1pexp_a = a;
    else log1pexp_a = log1p_exp(a);

    return rr * t - ((rr + dd) / kappa) * log1pexp_a;
  }
}

data {
  int<lower=1> N;
  vector[N] t;
  vector[N] y_log10;
  array[N] int<lower=0, upper=1> is_left_censored;
  real y_lod_log10;
  int<lower=1> K_site;
  array[N] int<lower=1, upper=K_site> site_id;
  array[N] int<lower=0, upper=1> hiv;

  real prior_log_r_mean;
  real<lower=0> prior_log_r_sd;
  real prior_log_d_mean;
  real<lower=0> prior_log_d_sd;
  real prior_m_mean;
  real<lower=0> prior_m_sd;
  real prior_log_kappa_mean;
  real<lower=0> prior_log_kappa_sd;
}

parameters {
  real<lower=-5, upper=2> log_r;
  real<lower=-5, upper=2> log_d;
  real<lower=-20, upper=40> m;
  real<lower=-5, upper=3> log_kappa;

  vector[K_site] alpha_site;
  real beta_hiv;
  real<lower=0> sigma_y;
}

transformed parameters {
  real<lower=0> r = exp(log_r);
  real<lower=0> d = exp(log_d);
  real<lower=0> kappa = exp(log_kappa);
}

model {
  log_r ~ normal(prior_log_r_mean, prior_log_r_sd);
  log_d ~ normal(prior_log_d_mean, prior_log_d_sd);
  m ~ normal(prior_m_mean, prior_m_sd);
  log_kappa ~ normal(prior_log_kappa_mean, prior_log_kappa_sd);

  alpha_site ~ normal(6, 3);
  beta_hiv ~ normal(0, 2);
  sigma_y ~ normal(0, 1.5);

  for (n in 1:N) {
    real mu = alpha_site[site_id[n]] + beta_hiv * hiv[n]
              + log_viral_shape_one(t[n], log_r, log_d, m, log_kappa) / log(10);
    if (is_left_censored[n] == 1) {
      target += normal_lcdf(y_lod_log10 | mu, sigma_y);
    } else {
      y_log10[n] ~ normal(mu, sigma_y);
    }
  }
}

generated quantities {
  vector[N] log_lik;

  for (n in 1:N) {
    real mu = alpha_site[site_id[n]] + beta_hiv * hiv[n]
              + log_viral_shape_one(t[n], log_r, log_d, m, log_kappa) / log(10);
    if (is_left_censored[n] == 1) {
      log_lik[n] = normal_lcdf(y_lod_log10 | mu, sigma_y);
    } else {
      log_lik[n] = normal_lpdf(y_log10[n] | mu, sigma_y);
    }
  }
}
