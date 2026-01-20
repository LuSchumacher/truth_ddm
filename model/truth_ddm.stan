functions {
  real partial_sum_fullddm(
    array[] real rt_slice, int start, int end, matrix a, matrix v, matrix bias,
    array[] real t0, array[] int subject_id, array[] int resp, array[] int condition
  ) {
    real ans = 0;
    for (i in start:end) {
      if (resp[i] == 1) {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[condition[i], subject_id[i]],
          t0[i], bias[condition[i], subject_id[i]],
          v[condition[i], subject_id[i]], 0, 0, 0
        );
      } else {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[condition[i], subject_id[i]],
          t0[i], 1-bias[condition[i], subject_id[i]],
          -v[condition[i], subject_id[i]], 0, 0, 0
        );
      }
    }
    return ans;
  }
}

data {
  int<lower=1>                    T;
  int<lower=1>                    N;
  array[T] real<lower=0>          rt;
  array[T] int<lower=1>           subject_id;
  array[T] int<lower=0, upper=1>  resp;
  array[T] int<lower=-1, upper=1> truth;
  array[T] int<lower=-1, upper=1> repetition;
  array[T] int<lower=1, upper=4>  condition;
  matrix<lower=0>[4, N]           minRT;
}

parameters {
  // Group-level intercepts and regression betas
  real      v_intercept;
  vector[3] v_betas;
  real      sigma_v;
  real      a_intercept;
  vector[3] a_betas;
  real      sigma_a;
  real      bias_intercept;
  vector[3] bias_betas;
  real      sigma_bias;
  real      ndt_intercept;
  vector[3] ndt_betas;
  real      sigma_ndt;
  real      ndt_var_intercept;
  vector[3] ndt_var_betas;
  real      sigma_ndt_var;

  // Subject-level random effects
  vector[N] z_v;
  vector[N] z_a;
  vector[N] z_bias;
  vector[N] z_ndt;
  vector[N] z_ndt_var;

  vector<lower=0, upper=1>[T] s; // for t0 per trial
  matrix<lower=0, upper=1>[4, N] trel; // for minRT
}

transformed parameters {
  real<lower=0> s_v       = log1p_exp(sigma_v);
  real<lower=0> s_a       = log1p_exp(sigma_a);
  real<lower=0> s_bias    = log1p_exp(sigma_bias);
  real<lower=0> s_ndt     = log1p_exp(sigma_ndt);
  real<lower=0> s_ndt_var = log1p_exp(sigma_ndt_var);
  
  // Condition specific parameters
  vector[4] mu_v;
  vector[4] mu_a;
  vector[4] mu_bias;
  vector[4] mu_ndt;
  vector[4] mu_ndt_var;
  // Subject specific parameters
  matrix[4, N] v;
  matrix[4, N] a;
  matrix[4, N] bias;
  matrix[4, N] ndt;
  matrix[4, N] ndt_var;
  
  for (cell in 1:4) {
    // condition (1-4) coding
    vector[3] x_cell;
    x_cell[1] = (cell == 1 || cell == 3) ? -1 : 1; // repetition
    x_cell[2] = (cell >= 3) ? -1 : 1;              // truth
    x_cell[3] = x_cell[1] * x_cell[2];             // interaction
  
    mu_v[cell]       = v_intercept       + dot_product(v_betas, x_cell);
    mu_a[cell]       = a_intercept       + dot_product(a_betas, x_cell);
    mu_bias[cell]    = bias_intercept    + dot_product(bias_betas, x_cell);
    mu_ndt[cell]     = ndt_intercept     + dot_product(ndt_betas, x_cell);
    mu_ndt_var[cell] = ndt_var_intercept + dot_product(ndt_var_betas, x_cell);
    
    // Subject specific parameters
    // for (subj in 1:N) {
    v[cell]       = (mu_v[cell] + s_v * z_v)';
    a[cell]       = (log1p_exp(mu_a[cell] + s_a * z_a))';
    bias[cell]    = (inv_logit(mu_bias[cell] + s_bias * z_bias))';
    ndt[cell]     = minRT[cell].* trel[cell];
    ndt_var[cell] = (log1p_exp(mu_ndt_var[cell] + s_ndt_var * z_ndt_var))';
    // }
  
  }
}

model {
  v_intercept          ~ normal(0, 1);
  v_betas              ~ normal(0, 0.5);
  sigma_v              ~ normal(0, 1);
  a_intercept          ~ normal(2, 1);
  a_betas              ~ normal(0, 0.5);
  sigma_a              ~ normal(0, 1);
  bias_intercept       ~ normal(0, 1);
  bias_betas           ~ normal(0, 0.5);
  sigma_bias           ~ normal(0, 1);
  ndt_intercept        ~ normal(1, 1);
  ndt_betas            ~ normal(0, 0.5);
  sigma_ndt            ~ normal(0, 1);
  ndt_var_intercept    ~ normal(1, 1);
  ndt_var_betas        ~ normal(0, 0.5);
  sigma_ndt_var        ~ normal(0, 1);
  to_vector(z_v)       ~ std_normal();
  to_vector(z_a)       ~ std_normal();
  to_vector(z_bias)    ~ std_normal();
  to_vector(z_ndt)     ~ std_normal();
  to_vector(z_ndt_var) ~ std_normal();
  
  for (i in 1:4) {
    ndt[i] ~ normal(mu_ndt[i], s_ndt);
  }
  
  s ~ uniform(0, 1);
  array[T] real t0;
  for (i in 1:T)  {
    t0[i] = ndt[condition[i], subject_id[i]] + s[i] * ndt_var[condition[i], subject_id[i]];
  }

  target += reduce_sum(
    partial_sum_fullddm, rt, 1,
    a, v, bias, t0, subject_id, resp, condition
  );
}

generated quantities {
  vector[4] transf_mu_v;
  vector[4] transf_mu_a;
  vector[4] transf_mu_ndt;
  vector[4] transf_mu_bias;
  vector[4] transf_mu_ndt_var;

  matrix[4, N] indiv_transf_v;
  matrix[4, N] indiv_transf_a;
  matrix[4, N] indiv_transf_bias;
  matrix[4, N] indiv_transf_ndt;
  matrix[4, N] indiv_transf_ndt_var;

  for (i in 1:4) {
    transf_mu_v[i]       = mu_v[i];
    transf_mu_a[i]       = log1p_exp(mu_a[i]);
    transf_mu_bias[i]    = inv_logit(mu_bias[i]);
    transf_mu_ndt[i]     = mu_ndt[i];
    transf_mu_ndt_var[i] = log1p_exp(mu_ndt_var[i]);
  }
  
  indiv_transf_v = v;
  indiv_transf_a = a;
  indiv_transf_bias = bias;
  indiv_transf_ndt = ndt;
  indiv_transf_ndt_var = ndt_var;
}
