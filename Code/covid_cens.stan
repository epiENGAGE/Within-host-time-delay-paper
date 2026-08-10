functions {
  real viral_shape_one(real t, real log_r, real log_d, real m, real log_kappa) {
    real rr = exp(log_r);
    real dd = exp(log_d);
    real kappa = exp(log_kappa);
    real a = log(rr / dd) + kappa * (t - m);
    real log1pexp_a;

    if (a > 40) log1pexp_a = a;
    else log1pexp_a = log1p_exp(a);

    return exp(fmin(rr * t - ((rr + dd) / kappa) * log1pexp_a, 700));
  }

  vector normalize_density(vector raw, real dx) {
    int N = num_elements(raw);
    vector[N] out;
    real z = sum(raw) * dx;
    if (z > 0) out = raw / z;
    else out = rep_vector(1.0 / (N * dx), N);
    return out;
  }

  vector normalize_mass(vector raw) {
    int N = num_elements(raw);
    vector[N] out;
    real z = sum(raw);
    if (z > 0) out = raw / z;
    else out = rep_vector(1.0 / N, N);
    return out;
  }

  vector cumulative_trapezoid(vector y, real dx) {
    int N = num_elements(y);
    vector[N] out;
    out[1] = 0;
    for (i in 2:N) {
      out[i] = out[i - 1] + dx * 0.5 * (y[i - 1] + y[i]);
    }
    return out;
  }
}

data {
  int<lower=1> N_cat;
  array[N_cat] int counts;
  vector[N_cat] serial_cat;

  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N_cat> obs_cat_pt;

  int<lower=1> N_l;
  vector[N_l] l_grid;
  real<lower=0> dx_l;

  int<lower=1> N_u;
  vector[N_u] u_grid;
  real<lower=0> dx_u;

  int<lower=1> N_gt;
  vector[N_gt] gt_grid;

  int<lower=1> N_ip;
  vector[N_ip] ip_grid;
  real<lower=0> dx_ip;

  int<lower=1> N_t;
  vector[N_t] t_grid;

  int<lower=1> N_delta;
  vector[N_delta] delta_grid;

  array[N_l, N_u] int<lower=1, upper=N_t> lu_t_idx;
  array[N_l, N_u] int<lower=1, upper=N_gt> gt_idx;

  int<lower=1> N_ip_pair;
  array[N_ip_pair] int<lower=1, upper=N_delta> pair_delta;
  array[N_ip_pair] int<lower=1, upper=N_ip> pair_i1;
  array[N_ip_pair] int<lower=1, upper=N_ip> pair_i2;

  int<lower=1> N_coarse;
  array[N_coarse] int<lower=1, upper=N_cat> coarse_cat;
  array[N_coarse] int<lower=1, upper=N_gt> coarse_g;
  array[N_coarse] int<lower=1, upper=N_delta> coarse_delta;
  vector<lower=0>[N_coarse] coarse_w;

  vector<lower=0>[N_ip] park_weight_ip;

  int<lower=1> N_z;
  vector[N_z] z_std_grid;
  vector<lower=0>[N_z] z_std_w;

  real prior_log_lambda1_mean;
  real<lower=0> prior_log_lambda1_sd;

  real prior_log_lambda_I_mean;
  real<lower=0> prior_log_lambda_I_sd;

  real prior_log_r_mean;
  real<lower=0> prior_log_r_sd;

  real prior_log_d_mean;
  real<lower=0> prior_log_d_sd;

  real prior_m_mean;
  real<lower=0> prior_m_sd;

  real prior_log_kappa_mean;
  real<lower=0> prior_log_kappa_sd;

  real prior_log_sigma_z_mean;
  real<lower=0> prior_log_sigma_z_sd;
}

parameters {
  real<lower=-8, upper=3> log_lambda1;
  real<lower=-8, upper=3> log_lambda_I;

  real<lower=-5, upper=2> log_r;
  real<lower=-5, upper=2> log_d;
  real<lower=0.25, upper=12> m;
  real<lower=-5, upper=3> log_kappa;
  real<lower=-6, upper=1.5> log_sigma_z;
}

transformed parameters {
}

model {
  log_lambda1 ~ normal(prior_log_lambda1_mean, prior_log_lambda1_sd);
  log_lambda_I ~ normal(prior_log_lambda_I_mean, prior_log_lambda_I_sd);
  log_r ~ normal(prior_log_r_mean, prior_log_r_sd);
  log_d ~ normal(prior_log_d_mean, prior_log_d_sd);
  m ~ normal(prior_m_mean, prior_m_sd);
  log_kappa ~ normal(prior_log_kappa_mean, prior_log_kappa_sd);
  log_sigma_z ~ normal(prior_log_sigma_z_mean, prior_log_sigma_z_sd);

  {
    real lambda1 = exp(log_lambda1);
    real lambda_I = exp(log_lambda_I);
    real sigma_z = exp(log_sigma_z);

    vector[N_l] base_l;
    vector[N_t] base_t;
    vector[N_ip] base_ip;

    matrix[N_ip, N_z] fIP_z;
    vector[N_ip] fIP2;
    matrix[N_delta, N_z] C_delta_z;
    array[N_z] matrix[N_cat, N_gt] H_cat_g_z;
    vector[N_cat] pi_cat;
    real norm_pi;

    for (jl in 1:N_l) base_l[jl] = viral_shape_one(l_grid[jl], log_r, log_d, m, log_kappa);
    for (jt in 1:N_t) base_t[jt] = viral_shape_one(t_grid[jt], log_r, log_d, m, log_kappa);
    for (ii in 1:N_ip) base_ip[ii] = viral_shape_one(ip_grid[ii], log_r, log_d, m, log_kappa);

    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);
      vector[N_ip] hI = ez * lambda_I * base_ip;
      vector[N_ip] HI = cumulative_trapezoid(hI, dx_ip);
      vector[N_ip] rawIP;
      for (ii in 1:N_ip) rawIP[ii] = hI[ii] * exp(-HI[ii]);
      fIP_z[, jz] = normalize_density(rawIP, dx_ip);
    }

    fIP2 = rep_vector(0, N_ip);
    for (jz in 1:N_z) fIP2 += z_std_w[jz] * fIP_z[, jz];
    fIP2 = normalize_density(fIP2, dx_ip);

    C_delta_z = rep_matrix(0, N_delta, N_z);
    for (jz in 1:N_z) {
      for (qq in 1:N_ip_pair) {
        int dd = pair_delta[qq];
        int ii = pair_i1[qq];
        int jj = pair_i2[qq];
        C_delta_z[dd, jz] += park_weight_ip[ii] * fIP_z[ii, jz] * fIP2[jj] * dx_ip * dx_ip;
      }
    }

    for (jz in 1:N_z) H_cat_g_z[jz] = rep_matrix(0, N_cat, N_gt);
    for (jz in 1:N_z) {
      for (qq in 1:N_coarse) {
        H_cat_g_z[jz][coarse_cat[qq], coarse_g[qq]] += C_delta_z[coarse_delta[qq], jz] * coarse_w[qq];
      }
    }

    pi_cat = rep_vector(0, N_cat);
    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);
      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid(hL, dx_l);
      vector[N_l] rawL;
      vector[N_l] normL;
      vector[N_u] rawU;
      vector[N_u] normU;

      for (jl in 1:N_l) rawL[jl] = hL[jl] * exp(-HL[jl]);
      normL = normalize_density(rawL, dx_l);

      for (jl in 1:N_l) {
        for (ju in 1:N_u) rawU[ju] = base_t[lu_t_idx[jl, ju]];
        normU = normalize_density(rawU, dx_u);

        for (ju in 1:N_u) {
          int kg = gt_idx[jl, ju];
          real lu_mass = normL[jl] * dx_l * normU[ju] * dx_u;
          for (kc in 1:N_cat) {
            pi_cat[kc] += z_std_w[jz] * lu_mass * H_cat_g_z[jz][kc, kg];
          }
        }
      }
    }

    norm_pi = sum(pi_cat);
    for (kc in 1:N_cat) pi_cat[kc] = pi_cat[kc] / fmax(norm_pi, 1e-300);

    for (kc in 1:N_cat) {
      if (counts[kc] > 0) target += counts[kc] * log(fmax(pi_cat[kc], 1e-300));
    }
  }
}

generated quantities {
  vector[N_obs] log_lik;
  vector[N_cat] log_lik_grouped;
  vector[N_gt] generation_prob_grid;
  real mean_gen;
  real sd_gen;

  {
    real lambda1 = exp(log_lambda1);
    real lambda_I = exp(log_lambda_I);
    real sigma_z = exp(log_sigma_z);

    vector[N_l] base_l;
    vector[N_t] base_t;
    vector[N_ip] base_ip;

    matrix[N_ip, N_z] fIP_z;
    vector[N_ip] fIP2;
    matrix[N_delta, N_z] C_delta_z;
    array[N_z] matrix[N_cat, N_gt] H_cat_g_z;
    vector[N_cat] pi_cat;
    real norm_pi;
    vector[N_gt] fG_mix_mass;

    for (jl in 1:N_l) base_l[jl] = viral_shape_one(l_grid[jl], log_r, log_d, m, log_kappa);
    for (jt in 1:N_t) base_t[jt] = viral_shape_one(t_grid[jt], log_r, log_d, m, log_kappa);
    for (ii in 1:N_ip) base_ip[ii] = viral_shape_one(ip_grid[ii], log_r, log_d, m, log_kappa);

    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);
      vector[N_ip] hI = ez * lambda_I * base_ip;
      vector[N_ip] HI = cumulative_trapezoid(hI, dx_ip);
      vector[N_ip] rawIP;
      for (ii in 1:N_ip) rawIP[ii] = hI[ii] * exp(-HI[ii]);
      fIP_z[, jz] = normalize_density(rawIP, dx_ip);
    }

    fIP2 = rep_vector(0, N_ip);
    for (jz in 1:N_z) fIP2 += z_std_w[jz] * fIP_z[, jz];
    fIP2 = normalize_density(fIP2, dx_ip);

    C_delta_z = rep_matrix(0, N_delta, N_z);
    for (jz in 1:N_z) {
      for (qq in 1:N_ip_pair) {
        int dd = pair_delta[qq];
        int ii = pair_i1[qq];
        int jj = pair_i2[qq];
        C_delta_z[dd, jz] += park_weight_ip[ii] * fIP_z[ii, jz] * fIP2[jj] * dx_ip * dx_ip;
      }
    }

    for (jz in 1:N_z) H_cat_g_z[jz] = rep_matrix(0, N_cat, N_gt);
    for (jz in 1:N_z) {
      for (qq in 1:N_coarse) {
        H_cat_g_z[jz][coarse_cat[qq], coarse_g[qq]] += C_delta_z[coarse_delta[qq], jz] * coarse_w[qq];
      }
    }

    pi_cat = rep_vector(0, N_cat);
    fG_mix_mass = rep_vector(0, N_gt);
    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);
      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid(hL, dx_l);
      vector[N_l] rawL;
      vector[N_l] normL;
      vector[N_u] rawU;
      vector[N_u] normU;

      for (jl in 1:N_l) rawL[jl] = hL[jl] * exp(-HL[jl]);
      normL = normalize_density(rawL, dx_l);

      for (jl in 1:N_l) {
        for (ju in 1:N_u) rawU[ju] = base_t[lu_t_idx[jl, ju]];
        normU = normalize_density(rawU, dx_u);

        for (ju in 1:N_u) {
          int kg = gt_idx[jl, ju];
          real lu_mass = normL[jl] * dx_l * normU[ju] * dx_u;
          fG_mix_mass[kg] += z_std_w[jz] * lu_mass;
          for (kc in 1:N_cat) {
            pi_cat[kc] += z_std_w[jz] * lu_mass * H_cat_g_z[jz][kc, kg];
          }
        }
      }
    }

    norm_pi = sum(pi_cat);
    for (kc in 1:N_cat) pi_cat[kc] = pi_cat[kc] / fmax(norm_pi, 1e-300);

    for (n in 1:N_obs) log_lik[n] = log(fmax(pi_cat[obs_cat_pt[n]], 1e-300));
    for (kc in 1:N_cat) log_lik_grouped[kc] = counts[kc] * log(fmax(pi_cat[kc], 1e-300));

    generation_prob_grid = normalize_mass(fG_mix_mass);
    mean_gen = dot_product(gt_grid, generation_prob_grid);
    sd_gen = sqrt(dot_product(square(gt_grid - rep_vector(mean_gen, N_gt)), generation_prob_grid));
  }
}
