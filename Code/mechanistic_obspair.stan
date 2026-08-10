functions {
  real viral_shape_one(real t, real log_r, real log_d, real m, real log_kappa) {
    real rr = exp(log_r);
    real dd = exp(log_d);
    real kappa = exp(log_kappa);

    real a = log(rr / dd) + kappa * (t - m);
    real log1pexp_a;

    if (a > 40) {
      log1pexp_a = a;
    } else {
      log1pexp_a = log1p_exp(a);
    }

    return exp(fmin(rr * t - ((rr + dd) / kappa) * log1pexp_a, 700));
  }

  vector normalize_density(vector raw, real dx) {
    real z = sum(raw) * dx;
    int N = num_elements(raw);
    vector[N] out;

    if (z > 0) {
      out = raw / z;
    } else {
      out = rep_vector(1.0 / (N * dx), N);
    }

    return out;
  }

  vector normalize_mass(vector raw) {
    real z = sum(raw);
    int N = num_elements(raw);
    vector[N] out;

    if (z > 0) {
      out = raw / z;
    } else {
      out = rep_vector(1.0 / N, N);
    }

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
  int<lower=1> N_si;
  array[N_si] int<lower=1> obs_index;
  array[N_si] int<lower=0> counts;

  int<lower=1> N_norm;
  vector[N_norm] si_grid;
  real<lower=0> dx_si;

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
  real<lower=0> dx_t;

  int<lower=1> N_delta;
  vector[N_delta] delta_grid;

  array[N_l] int<lower=1, upper=N_t> l_t_idx;
  array[N_l, N_u] int<lower=1, upper=N_t> lu_t_idx;
  array[N_l, N_u] int<lower=1, upper=N_gt> gt_idx;
  array[N_norm, N_gt] int<lower=1, upper=N_delta> si_g_delta_idx;

  int<lower=1> N_ip_sparse;
  array[N_ip_sparse] int<lower=1, upper=N_delta> ip_sparse_delta;
  array[N_ip_sparse] int<lower=1, upper=N_ip> ip_sparse_i1;
  array[N_ip_sparse] int<lower=1, upper=N_ip> ip_sparse_i2;

  vector[N_ip] park_weight_ip;

  int<lower=1> N_z;
  vector[N_z] z_std_grid;
  vector[N_z] z_std_w;

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
    real dx_lu = dx_l * dx_u;

    vector[N_l] base_l;
    vector[N_t] base_t;
    vector[N_ip] base_ip;

    matrix[N_ip, N_z] fIP_z;
    vector[N_ip] fIP2;

    matrix[N_delta, N_z] C_delta_z;
    matrix[N_gt, N_z] fG_mass_z;
    vector[N_norm] log_unnorm_si;
    real log_norm;

    for (jl in 1:N_l) {
      base_l[jl] = viral_shape_one(l_grid[jl], log_r, log_d, m, log_kappa);
    }

    for (jt in 1:N_t) {
      base_t[jt] = viral_shape_one(t_grid[jt], log_r, log_d, m, log_kappa);
    }

    for (ii in 1:N_ip) {
      base_ip[ii] = viral_shape_one(ip_grid[ii], log_r, log_d, m, log_kappa);
    }

    fG_mass_z = rep_matrix(0, N_gt, N_z);

    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);

      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid(hL, dx_l);
      vector[N_l] rawL = hL .* exp(-HL);
      vector[N_l] normL = normalize_density(rawL, dx_l);

      vector[N_ip] hI = ez * lambda_I * base_ip;
      vector[N_ip] HI = cumulative_trapezoid(hI, dx_ip);
      vector[N_ip] rawIP = hI .* exp(-HI);

      vector[N_u] rawU;
      vector[N_u] normU;

      for (jl in 1:N_l) {
        for (ju in 1:N_u) {
          int idx = lu_t_idx[jl, ju];
          rawU[ju] = base_t[idx];
        }

        normU = normalize_density(rawU, dx_u);

        for (ju in 1:N_u) {
          int kg = gt_idx[jl, ju];
          fG_mass_z[kg, jz] += normL[jl] * normU[ju] * dx_lu;
        }
      }

      fG_mass_z[, jz] = normalize_mass(fG_mass_z[, jz]);
      fIP_z[, jz] = normalize_density(rawIP, dx_ip);
    }

    fIP2 = rep_vector(0, N_ip);
    for (jz in 1:N_z) {
      fIP2 += z_std_w[jz] * fIP_z[, jz];
    }
    fIP2 = normalize_density(fIP2, dx_ip);

    C_delta_z = rep_matrix(0, N_delta, N_z);
    for (jz in 1:N_z) {
      for (qq in 1:N_ip_sparse) {
        int dd = ip_sparse_delta[qq];
        int ii = ip_sparse_i1[qq];
        int jj = ip_sparse_i2[qq];

        C_delta_z[dd, jz] += park_weight_ip[ii]
                              * fIP_z[ii, jz]
                              * fIP2[jj]
                              * dx_ip;
      }
    }

    for (ks in 1:N_norm) {
      real acc_s = 0;

      for (kg in 1:N_gt) {
        int dd = si_g_delta_idx[ks, kg];

        for (jz in 1:N_z) {
          acc_s += z_std_w[jz]
                   * fG_mass_z[kg, jz]
                   * C_delta_z[dd, jz];
        }
      }

      log_unnorm_si[ks] = log(fmax(acc_s, 1e-300));
    }

    log_norm = log_sum_exp(log_unnorm_si) + log(dx_si);

    for (n in 1:N_si) {
      target += counts[n] * (log_unnorm_si[obs_index[n]] - log_norm);
    }
  }
}

generated quantities {
  real mean_gen;
  real sd_gen;

  {
    real lambda1 = exp(log_lambda1);
    real sigma_z = exp(log_sigma_z);
    real dx_lu = dx_l * dx_u;

    vector[N_l] base_l;
    vector[N_t] base_t;
    matrix[N_gt, N_z] fG_mass_z;
    vector[N_gt] fG_mix_mass;

    for (jl in 1:N_l) {
      base_l[jl] = viral_shape_one(l_grid[jl], log_r, log_d, m, log_kappa);
    }

    for (jt in 1:N_t) {
      base_t[jt] = viral_shape_one(t_grid[jt], log_r, log_d, m, log_kappa);
    }

    fG_mass_z = rep_matrix(0, N_gt, N_z);

    for (jz in 1:N_z) {
      real z = sigma_z * z_std_grid[jz];
      real ez = exp(z);

      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid(hL, dx_l);
      vector[N_l] rawL = hL .* exp(-HL);
      vector[N_l] normL = normalize_density(rawL, dx_l);

      vector[N_u] rawU;
      vector[N_u] normU;

      for (jl in 1:N_l) {
        for (ju in 1:N_u) {
          int idx = lu_t_idx[jl, ju];
          rawU[ju] = base_t[idx];
        }

        normU = normalize_density(rawU, dx_u);

        for (ju in 1:N_u) {
          int kg = gt_idx[jl, ju];
          fG_mass_z[kg, jz] += normL[jl] * normU[ju] * dx_lu;
        }
      }

      fG_mass_z[, jz] = normalize_mass(fG_mass_z[, jz]);
    }

    fG_mix_mass = rep_vector(0, N_gt);
    for (jz in 1:N_z) {
      fG_mix_mass += z_std_w[jz] * fG_mass_z[, jz];
    }
    fG_mix_mass = normalize_mass(fG_mix_mass);

    mean_gen = dot_product(gt_grid, fG_mix_mass);
    sd_gen = sqrt(dot_product(square(gt_grid - rep_vector(mean_gen, N_gt)), fG_mix_mass));
  }
}
