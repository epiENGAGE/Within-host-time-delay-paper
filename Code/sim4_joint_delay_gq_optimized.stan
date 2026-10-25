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

  array[N_l] int<lower=1, upper=N_t> lu_start_idx;
  array[N_l] int<lower=1, upper=N_gt> gt_start_idx;

  array[N_delta] int<lower=1, upper=N_ip> ip_i1_start;
  array[N_delta] int<lower=1, upper=N_ip> ip_i2_start;
  array[N_delta] int<lower=0, upper=N_ip> ip_overlap_len;

  array[N_norm] int<lower=1, upper=N_delta> si_delta_low_idx;

  vector[N_ip] park_weight_ip;

  int<lower=1> N_z;
  vector[N_z] z_std_grid;
  vector[N_z] z_std_w;

  // Kept in the GQ data schema so the same R data list can be reused.
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
  // Intentionally empty. This model is only used with standalone GQ.
}

generated quantities {
  real mean_lp;
  real sd_lp;
  real mean_ip;
  real sd_ip;
  real mean_gen;
  real sd_gen;
  real mean_si;
  real sd_si;

  {
    real lambda1 = exp(log_lambda1);
    real lambda_I = exp(log_lambda_I);
    real sigma_z = exp(log_sigma_z);
    real dx_lu = dx_l * dx_u;

    vector[N_l] base_l = viral_shape_vec(l_grid, log_r, log_d, m, log_kappa);
    vector[N_t] base_t = viral_shape_vec(t_grid, log_r, log_d, m, log_kappa);
    vector[N_ip] base_ip = viral_shape_vec(ip_grid, log_r, log_d, m, log_kappa);

    matrix[N_l, N_z] fL_z;
    matrix[N_ip, N_z] fIP_z;
    matrix[N_gt, N_z] fG_mass_z;
    vector[N_ip] fIP2;
    matrix[N_delta, N_z] C_delta_z;

    vector[N_l] fL_mass;
    vector[N_ip] fIP_mass;
    vector[N_gt] fG_mass;
    vector[N_norm] p_si_mass;

    fG_mass_z = rep_matrix(0.0, N_gt, N_z);

    for (jz in 1:N_z) {
      real ez = exp(sigma_z * z_std_grid[jz]);

      vector[N_l] hL = ez * lambda1 * base_l;
      vector[N_l] HL = cumulative_trapezoid_fast(hL, dx_l);
      vector[N_l] normL = normalize_density_fast(hL .* exp(-HL), dx_l);

      vector[N_ip] hI = ez * lambda_I * base_ip;
      vector[N_ip] HI = cumulative_trapezoid_fast(hI, dx_ip);

      fL_z[, jz] = normL;
      fIP_z[, jz] = normalize_density_fast(hI .* exp(-HI), dx_ip);

      {
        vector[N_gt] gcol = rep_vector(0.0, N_gt);
        for (jl in 1:N_l) {
          int t0 = lu_start_idx[jl];
          int g0 = gt_start_idx[jl];
          vector[N_u] normU =
            normalize_density_fast(segment(base_t, t0, N_u), dx_u);

          gcol[g0:(g0 + N_u - 1)] =
            gcol[g0:(g0 + N_u - 1)] + normL[jl] * normU * dx_lu;
        }
        fG_mass_z[, jz] = normalize_mass_fast(gcol);
      }
    }

    fIP2 = normalize_density_fast(fIP_z * z_std_w, dx_ip);

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

    fL_mass = normalize_mass_fast((fL_z * z_std_w) * dx_l);
    fIP_mass = normalize_mass_fast(fIP2 * dx_ip);
    fG_mass = normalize_mass_fast(fG_mass_z * z_std_w);

    {
      vector[N_norm] unnorm_si;
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
      p_si_mass = normalize_mass_fast(unnorm_si);
    }

    mean_lp = dot_product(l_grid, fL_mass);
    sd_lp = sqrt(dot_product(square(l_grid - rep_vector(mean_lp, N_l)), fL_mass));

    mean_ip = dot_product(ip_grid, fIP_mass);
    sd_ip = sqrt(dot_product(square(ip_grid - rep_vector(mean_ip, N_ip)), fIP_mass));

    mean_gen = dot_product(gt_grid, fG_mass);
    sd_gen = sqrt(dot_product(square(gt_grid - rep_vector(mean_gen, N_gt)), fG_mass));

    mean_si = dot_product(si_grid, p_si_mass);
    sd_si = sqrt(dot_product(square(si_grid - rep_vector(mean_si, N_norm)), p_si_mass));
  }
}
