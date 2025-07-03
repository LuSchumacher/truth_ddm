library(tidyverse)
library(magrittr)
library(cmdstanr)
library(bayesplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

init_fun <- function(chains = 4){
  L = list()
  for (i in 1:chains) {
    L[[i]] = list(
      mu_v        = 1.0 + runif(4, -0.5, 0.5),
      sigma_v     = 0.1 + runif(4, -0.05, 0.05),
      mu_a        = 2.5 + runif(4, -0.2, 0.2),
      sigma_a     = 0.1 + runif(4, -0.05, 0.05),
      mu_bias     = 0.0 + runif(4, -0.05, 0.05),
      sigma_bias  = 0.5 + runif(4, -0.01, 0.01),
      mu_ndt      = 0.1 + runif(4, -0.01, 0.01),
      sigma_ndt   = 0.01 + runif(4, -0.01, 0.01),
      z_v         = matrix(runif(4*N, -1, 1), nrow=4),
      z_a         = matrix(runif(4*N, -1, 1), nrow=4),
      z_ndt       = matrix(runif(4*N, -0.2, 0.2), nrow=4),
      z_bias      = matrix(runif(4*N, -1, 1), nrow=4),
      trel        = matrix(runif(4*N, 0.01, 0.99), nrow=4),
      s           = runif(`T`, 0, 0.05)
    )
  }
  return(L)
}

PARAM_NAMES <- c(
  "transf_mu_v[1]", "transf_mu_v[2]", "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]", "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]", "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]", "transf_mu_bias[3]", "transf_mu_bias[4]"
)

m_full <- cmdstan_model(
  '../../model/full_model.stan',
  cpp_options = list(stan_threads = T)
)

################################################################################
# SESSION 1
################################################################################
df_session_1 <- read.csv('../../data/data_session_1.csv')

N <- length(unique(df_session_1$id))
`T` <- nrow(df_session_1)
stan_data = list(
  `T`           = nrow(df_session_1),
  N             = N,
  subject_id    = df_session_1$id,
  resp          = df_session_1$resp,
  factual_truth = df_session_1$factual_truth,
  condition     = df_session_1$condition,
  rt            = df_session_1$rt,
  minRT         = tapply(df_session_1$rt, list(df_session_1$condition, df_session_1$id), min)
)

fit_session_1 <- m_full$sample(
  data = stan_data,
  init = init_fun(),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 50,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

print(fit_session_1$summary(variables = PARAM_NAMES), n=100)
mcmc_trace(fit_session_1$draws(inc_warmup = TRUE), n_warmup = 3000, pars=PARAM_NAMES)

fit_session_1$save_object("../../fits/fit_session_1_new.rds")

################################################################################
# SESSION 2
################################################################################
df_session_2 <- read.csv('../../data/data_session_2.csv')

N <- length(unique(df_session_2$id))
`T` <- nrow(df_session_2)
stan_data = list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_session_2$id,
  resp          = df_session_2$resp,
  factual_truth = df_session_2$factual_truth,
  condition     = df_session_2$condition,
  rt            = df_session_2$rt,
  minRT         = tapply(df_session_2$rt, list(df_session_2$condition, df_session_2$id), min)
)

fit_session_2 <- m_full$sample(
  data = stan_data,
  init = init_fun(),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 50,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

fit_session_2$summary(variables = PARAM_NAMES)
fit_session_2$save_object("../../fits/fit_session_2_new.rds")

################################################################################
# EXPERIMENT 2
################################################################################
df_exp_2 <- read.csv('../../data/data_exp_2.csv')

N <- length(unique(df_exp_2$id))
`T` <- nrow(df_exp_2)
stan_data = list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_exp_2$id,
  resp          = df_exp_2$resp,
  factual_truth = df_exp_2$factual_truth,
  condition     = df_exp_2$condition,
  rt            = df_exp_2$rt,
  minRT         = tapply(df_exp_2$rt, list(df_exp_2$condition, df_exp_2$id), min)
)

fit_exp_2 <- m_full$sample(
  data = stan_data,
  init = init_fun(),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 50,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

fit_exp_2$summary(variables = PARAM_NAMES)
fit_exp_2$save_object("../../fits/fit_exp_2_new.rds")

