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
  
  row_vector softplus(row_vector x) {
      int N = num_elements(x);
      row_vector[N] y;
      for (i in 1:N) {
        y[i] = (x[i] > 20) ? x[i] + log1p(exp(-x[i])) : log1p(exp(x[i]));
      }
      return y;
    }

  vector softplus(vector x) {
    int N = num_elements(x);
    vector[N] y;
    for (i in 1:N) {
      y[i] = (x[i] > 20) ? x[i] + log1p(exp(-x[i])) : log1p(exp(x[i]));
    }
    return y;
  }

  real softplus_r(real x) {
     return (x > 20) ? x + log1p(exp(-x)) : log1p(exp(x));
  }
}

data {
  int<lower=1>                   T;
  int<lower=1>                   N;
  array[T] real<lower=0>         rt;
  array[T] int<lower=1>          subject_id;
  array[T] int<lower=0, upper=1> resp;
  array[T] int<lower=1, upper=4> condition;
  matrix<lower=0>[4, N]          minRT;
}

parameters {
  // group-level parameters
  vector[4] mu_v;
  vector[4] sigma_v;
  vector[4] mu_a;
  vector[4] sigma_a;
  vector[4] mu_ndt;
  vector[4] sigma_ndt;
  vector[4] mu_bias;
  vector[4] sigma_bias;
  real mu_ndt_s;
  real sigma_ndt_s;

  matrix[4, N] z_v;
  matrix[4, N] z_a;
  matrix[4, N] z_bias;
  vector[N] z_ndt_s;
  
  vector<lower=0, upper=1>[T] s; // for t0 per trial
  matrix<lower=0, upper=1>[4, N] trel; // for minRT
}

transformed parameters {
  // subject-level parameters
  matrix[4, N] v;
  matrix[4, N] a;
  matrix[4, N] ndt;
  matrix[4, N] bias;
  vector[N] ndt_s;
  
  vector<lower=0>[4] s_v;
  vector<lower=0>[4] s_a;
  vector<lower=0>[4] s_ndt;
  vector<lower=0>[4] s_bias;
  real<lower=0> s_ndt_s;
  
  s_v     = softplus(sigma_v);
  s_a     = softplus(sigma_a);
  s_ndt   = softplus(sigma_ndt);
  s_bias  = softplus(sigma_bias);
  s_ndt_s = softplus_r(sigma_ndt_s);

	for (i in 1:4) {
    v[i]    = mu_v[i] + s_v[i] * z_v[i];
    a[i]    = softplus(mu_a[i] + s_a[i] * z_a[i,]);
    ndt[i]  = minRT[i].* trel[i];
    bias[i] = inv_logit(mu_bias[i] + s_bias[i] * z_bias[i,]);
	}
	
	ndt_s = softplus(mu_ndt_s + s_ndt_s * z_ndt_s);
	
}

model {
  array[T] real t0;
  // priors
  mu_v        ~ normal(2, 2);
  sigma_v     ~ normal(0, 1.5);
  mu_a        ~ normal(5, 3);
  sigma_a     ~ normal(0, 1.5);
  mu_bias     ~ normal(0, 0.5);
  sigma_bias  ~ normal(0, 1.5);
  mu_ndt      ~ normal(1, 2)T[0, ];
  sigma_ndt   ~ normal(0, 1.5);
  mu_ndt_s    ~ normal(1, 1)T[0, ];
  sigma_ndt_s ~ normal(0, 1.5);
  s           ~ uniform(0, 1);

  for (i in 1:4) {
  	z_v[i]    ~ std_normal();
  	z_a[i]    ~ std_normal();
  	z_bias[i] ~ std_normal();
  	ndt[i]    ~ normal(mu_ndt[i], s_ndt[i]);
  }
  
  z_ndt_s ~ std_normal();

  for (i in 1:T)  t0[i] = ndt[condition[i], subject_id[i]] + s[i] * ndt_s[subject_id[i]];
	
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
	real transf_mu_ndt_s;

  for (i in 1:4) {
    transf_mu_v[i]    = mu_v[i];
    transf_mu_ndt[i]  = mu_ndt[i];
    transf_mu_a[i]    = softplus_r(mu_a[i]);
    transf_mu_bias[i] = inv_logit(mu_bias[i]);
  }
  transf_mu_ndt_s = softplus_r(mu_ndt_s);
}
