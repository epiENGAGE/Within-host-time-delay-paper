functions {
  vector viral_shape_vec(vector t, real log_r, real log_d, real m, real log_kappa) {
    int N = num_elements(t);
    real rr = exp(log_r);
    real dd = exp(log_d);
    real kappa = exp(log_kappa);
    vector[N] a = rep_vector(log(rr / dd) - kappa * m, N) + kappa * t;
    vector[N] log_v = rr * t - ((rr + dd) / kappa) * log1p_exp(a);
    return exp(log_v);
  }

  vector normalize_density_fast(vector raw, real dx) {
    return raw / fmax(sum(raw) * dx, 1e-300);
  }

  vector normalize_mass_fast(vector raw) {
    return raw / fmax(sum(raw), 1e-300);
  }

  vector cumulative_trapezoid_fast(vector y, real dx) {
    int N = num_elements(y);
    vector[N] out;
    out[1] = 0;
    if (N > 1) {
      out[2:N] = cumulative_sum(0.5 * dx * (y[1:(N - 1)] + y[2:N]));
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

  int<lower=1> N_delta;
  vector[N_delta] delta_grid;

  // Compact lookup arrays replacing the large scalar-scatter lookup tables.
  array[N_l] int<lower=1, upper=N_t> lu_start_idx;
  array[N_l] int<lower=1, upper=N_gt> gt_start_idx;

  // For each delta, IP2 = IP1 + delta. These define the two aligned
  // contiguous IP slices used in a single dot product.
  array[N_delta] int<lower=1, upper=N_ip> ip_i1_start;
  array[N_delta] int<lower=1, upper=N_ip> ip_i2_start;
  array[N_delta] int<lower=0, upper=N_ip> ip_overlap_len;

  // For each SI cell, the required C_delta values for increasing GT are a
  // contiguous block in reverse order.
  array[N_norm] int<lower=1, upper=N_delta> si_delta_low_idx;

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

    // Vectorised evaluation of the within-host curve.
    vector[N_l] base_l = viral_shape_vec(l_grid, log_r, log_d, m, log_kappa);
    vector[N_t] base_t = viral_shape_vec(t_grid, log_r, log_d, m, log_kappa);
    vector[N_ip] base_ip = viral_shape_vec(ip_grid, log_r, log_d, m, log_kappa);

    matrix[N_ip, N_z] fIP_z;
    matrix[N_gt, N_z] fG_mass_z;
    vector[N_ip] fIP2;
    matrix[N_delta, N_z] C_delta_z;
    vector[N_norm] unnorm_si;

    fG_mass_z = rep_matrix(0.0, N_gt, N_z);

    // Build LP/IP/GT conditional on each heterogeneity quadrature point.
    for (jz in 1:N_z) {
      real ez = exp(sigma_z * z_std_grid[jz]);

      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid_fast(hL, dx_l);
      vector[N_l] normL = normalize_density_fast(hL .* exp(-HL), dx_l);

      vector[N_ip] hI = ez * lambda_I * base_ip;
      vector[N_ip] HI = cumulative_trapezoid_fast(hI, dx_ip);
      fIP_z[, jz] = normalize_density_fast(hI .* exp(-HI), dx_ip);

      {
        vector[N_gt] gcol = rep_vector(0.0, N_gt);

        // Because l/u/t/GT share the same grid spacing, each conditional-U
        // trajectory is a contiguous slice and each GT contribution is a
        // contiguous slice: no scalar gather/scatter loop is needed.
        for (jl in 1:N_l) {
          int t0 = lu_start_idx[jl];
          int g0 = gt_start_idx[jl];
          vector[N_u] rawU = segment(base_t, t0, N_u);
          vector[N_u] normU = normalize_density_fast(rawU, dx_u);

          gcol[g0:(g0 + N_u - 1)] =
            gcol[g0:(g0 + N_u - 1)] + normL[jl] * normU * dx_lu;
        }

        fG_mass_z[, jz] = normalize_mass_fast(gcol);
      }
    }

    // Infectee IP mixture.
    fIP2 = fIP_z * z_std_w;
    fIP2 = normalize_density_fast(fIP2, dx_ip);

    // IP2-IP1 cross-correlation. The old implementation performed hundreds
    // of thousands of scalar AD scatter-adds. Here every delta is one
    // dot-product over two contiguous IP slices.
    C_delta_z = rep_matrix(0.0, N_delta, N_z);
    for (jz in 1:N_z) {
      vector[N_ip] weighted_ip1 = park_weight_ip .* fIP_z[, jz];

      for (dd in 1:N_delta) {
        int L = ip_overlap_len[dd];
        if (L > 0) {
          C_delta_z[dd, jz] =
            dot_product(
              segment(weighted_ip1, ip_i1_start[dd], L),
              segment(fIP2, ip_i2_start[dd], L)
            ) * dx_ip;
        }
      }
    }

    // SI = GT + delta. For fixed SI, the required delta indices are
    // contiguous and reversed as GT increases, so each z contribution is
    // one dot product rather than an inner scalar loop over GT.
    for (ks in 1:N_norm) {
      real acc_s = 0.0;
      int dlow = si_delta_low_idx[ks];

      for (jz in 1:N_z) {
        acc_s += z_std_w[jz] *
          dot_product(
            fG_mass_z[, jz],
            reverse(segment(C_delta_z[, jz], dlow, N_gt))
          );
      }

      unnorm_si[ks] = fmax(acc_s, 1e-300);
    }

    // Same normalized-grid likelihood as before, but avoid log_sum_exp of
    // already-positive probabilities.
    {
      real norm_si = fmax(sum(unnorm_si) * dx_si, 1e-300);
      for (n in 1:N_si) {
        target += counts[n] * (log(unnorm_si[obs_index[n]]) - log(norm_si));
      }
    }
  }
}
