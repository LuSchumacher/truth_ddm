functions {
  real partial_sum_fullddm(
    array[] real rt_slice, int start, int end, matrix a, matrix v, matrix bias,
    matrix ndt, array[] int subject_id, array[] int resp, array[] int condition
  ) {
    real ans = 0;
    for (i in start:end) {
      if (resp[i] == 1) {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[condition[i], subject_id[i]],
          ndt[condition[i], subject_id[i]], bias[condition[i], subject_id[i]],
          v[condition[i], subject_id[i]], 0, 0, 0
        );
      } else {
        ans += wiener_lpdf(
          rt_slice[i+1-start] | a[condition[i], subject_id[i]],
          ndt[condition[i], subject_id[i]], 1-bias[condition[i], subject_id[i]],
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
  
  matrix[4, N] z_v;
  matrix[4, N] z_a;
  matrix[4, N] z_bias;

  matrix<lower=0, upper=1>[4, N] trel;
}

transformed parameters {
  // subject-level parameters
  matrix[4, N] v;
  matrix[4, N] a;
  matrix[4, N] ndt;
  matrix[4, N] bias;
  
  vector<lower=0>[4] s_v;
  vector<lower=0>[4] s_a;
  vector<lower=0>[4] s_ndt;
  vector<lower=0>[4] s_bias;
  
  s_v    = softmax(sigma_v);
  s_a    = softmax(sigma_a);
  s_ndt  = softmax(sigma_ndt);
  s_bias = softmax(sigma_bias);
  
	vector[4] transf_mu_v;
	vector[4] transf_mu_a;
	vector[4] transf_mu_ndt;
	vector[4] transf_mu_bias;

	for (i in 1:4) {
    v[i]    = mu_v[i] + s_v[i] * z_v[i];
    a[i]    = softplus(mu_a[i] + s_a[i] * z_a[i,]);
    ndt[i]  = minRT[i].* trel[i];
    bias[i] = inv_logit(mu_bias[i] + s_bias[i] * z_bias[i,]);
    
    // for the output
  	transf_mu_v[i]    = mu_v[i];
  	transf_mu_ndt[i]  = mu_ndt[i];
  	transf_mu_a[i]    = softplus_r(mu_a[i]);
  	transf_mu_bias[i] = inv_logit(mu_bias[i]);
	}
}

model {
  // priors
  mu_v       ~ normal(2, 2);
  sigma_v    ~ normal(0, 1.5);
  mu_a       ~ normal(5, 3);
  sigma_a    ~ normal(0, 1.5);
  
  mu_bias    ~ normal(0, 0.5);
  sigma_bias ~ normal(0, 1.5);
  mu_ndt     ~ normal(1, 2)T[0, ];
  sigma_ndt  ~ normal(0, 1.5);
  
  for (i in 1:4) {
  	z_v[i]    ~ std_normal();
  	z_a[i]    ~ std_normal();
  	z_bias[i] ~ std_normal();
  	ndt[i] ~ normal(mu_ndt[i], sigma_ndt[i]);
  }
	
  target += reduce_sum(
    partial_sum_fullddm, rt, 1,
    a, v, bias, ndt, subject_id, resp, condition
  );
}
