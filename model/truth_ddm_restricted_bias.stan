functions {
  real partial_sum_fullddm(
    array[] real rt_slice, int start, int end, matrix a, matrix v, matrix bias,
    array[] real t0, array[] int subject_id, array[] int resp,
    array[] int condition, array[] int truth
  ) {
    real ans = 0;
    for (i in start:end) {
      int c = condition[i];
      int t = truth[i];
      if (resp[i] == 1) {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[t, subject_id[i]],
          t0[i], bias[c, subject_id[i]],
          v[t, subject_id[i]], 0, 0, 0
        );
      } else {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[t, subject_id[i]],
          t0[i], 1 - bias[c, subject_id[i]],
          -v[t, subject_id[i]], 0, 0, 0
        );
      }
    }
    return ans;
  }
}

data {
  int<lower=1>                   T;
  int<lower=1>                   N;
  array[T] real<lower=0>         rt;
  array[T] int<lower=1>          subject_id;
  array[T] int<lower=0, upper=1> resp;
  array[T] int<lower=1, upper=2> truth;
  array[T] int<lower=1, upper=4> condition;
  matrix<lower=0>[2, N]          minRT;
}

parameters {
  // Group-level intercepts and regression betas
  real      v_intercept;
  real      v_beta;
  real      sigma_v;
  real      a_intercept;
  real      a_beta;
  real      sigma_a;
  real      bias_intercept;
  vector[3] bias_betas;
  real      sigma_bias;
  real      ndt_intercept;
  real      ndt_beta;
  real      sigma_ndt;
  real      ndt_var_intercept;
  real      ndt_var_beta;
  real      sigma_ndt_var;

  // Subject-level random effects
  vector[N] z_v;
  vector[N] z_a;
  vector[N] z_bias;
  vector[N] z_ndt;
  vector[N] z_ndt_var;

  vector<lower=0, upper=1>[T] s; // for t0 per trial
  matrix<lower=0, upper=1>[2, N] trel; // for minRT
}

transformed parameters {
  real<lower=0> s_v       = log1p_exp(sigma_v);
  real<lower=0> s_a       = log1p_exp(sigma_a);
  real<lower=0> s_bias    = log1p_exp(sigma_bias);
  real<lower=0> s_ndt     = log1p_exp(sigma_ndt);
  real<lower=0> s_ndt_var = log1p_exp(sigma_ndt_var);
  
  // Condition specific parameters
  vector[2] mu_v;
  vector[2] mu_a;
  vector[4] mu_bias;
  vector[2] mu_ndt;
  vector[2] mu_ndt_var;
  // Subject specific parameters
  matrix[2, N] v;
  matrix[2, N] a;
  matrix[4, N] bias;
  matrix[2, N] ndt;
  matrix[2, N] ndt_var;
  
  vector[2] effect_coding;
  effect_coding[1] = -1;
  effect_coding[2] =  1;
  
  for (t in 1:2) {
    mu_v[t]       = v_intercept + v_beta * effect_coding[t];
    mu_a[t]       = a_intercept + a_beta * effect_coding[t];
    mu_ndt[t]     = ndt_intercept + ndt_beta * effect_coding[t];
    mu_ndt_var[t] = ndt_var_intercept + ndt_var_beta * effect_coding[t];
    
    v[t] = mu_v[t] + s_v * z_v;
    a[t] = log1p_exp(mu_a[t] + s_a * z_a);
    ndt[t] = mu_ndt[t] + s_ndt * z_ndt;
    ndt_var[t] = log1p_exp(mu_ndt_var[t] + s_ndt_var * z_ndt_var);
  }
  
  for (cell in 1:4) {
    // condition (1-4) coding
    vector[3] x_cell;
    x_cell[1] = (cell == 1 || cell == 3) ? -1 : 1; // repetition
    x_cell[2] = (cell >= 3) ? -1 : 1;              // truth
    x_cell[3] = x_cell[1] * x_cell[2];             // interaction
    mu_bias[cell] = bias_intercept + dot_product(bias_betas, x_cell);
    bias[cell] = (inv_logit(mu_bias[cell] + s_bias * z_bias))';
  }
}

model {
  v_intercept          ~ normal(0, 1);
  v_beta               ~ normal(0, 0.5);
  sigma_v              ~ normal(0, 1);
  a_intercept          ~ normal(2, 1);
  a_beta               ~ normal(0, 0.5);
  sigma_a              ~ normal(0, 1);
  bias_intercept       ~ normal(0, 1);
  bias_betas           ~ normal(0, 0.5);
  sigma_bias           ~ normal(0, 1);
  ndt_intercept        ~ normal(1, 1);
  ndt_beta             ~ normal(0, 0.5);
  sigma_ndt            ~ normal(0, 1);
  ndt_var_intercept    ~ normal(1, 1);
  ndt_var_beta         ~ normal(0, 0.5);
  sigma_ndt_var        ~ normal(0, 1);
  to_vector(z_v)       ~ std_normal();
  to_vector(z_a)       ~ std_normal();
  to_vector(z_bias)    ~ std_normal();
  to_vector(z_ndt)     ~ std_normal();
  to_vector(z_ndt_var) ~ std_normal();
  
  for (i in 1:2) {
    ndt[i] ~ normal(mu_ndt[i], s_ndt);
  }
  
  s ~ uniform(0, 1);
  array[T] real t0;
  for (i in 1:T)  {
    t0[i] = ndt[truth[i], subject_id[i]] + s[i] * ndt_var[truth[i], subject_id[i]];
  }

  target += reduce_sum(
    partial_sum_fullddm, rt, 1,
    a, v, bias, t0, subject_id, resp, condition, truth
  );
}

generated quantities {
  vector[2] transf_mu_v;
	vector[2] transf_mu_a;
	vector[2] transf_mu_ndt;
	vector[4] transf_mu_bias;
	vector[2] transf_mu_ndt_var;

  for (i in 1:2) {
    transf_mu_v[i]       = mu_v[i];
    transf_mu_a[i]       = log1p_exp(mu_a[i]);
    transf_mu_ndt[i]     = mu_ndt[i];
    transf_mu_ndt_var[i] = log1p_exp(mu_ndt_var[i]);
  }
  for (i in 1:4) {
    transf_mu_bias[i]    = inv_logit(mu_bias[i]);
  }
}
