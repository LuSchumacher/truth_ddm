library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
library(cmdstanr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

PARAM_NAMES <- c("v", "a", "bias", "ndt", "ndt_s")

df <- read.csv('../data/exp2_data.csv')

df %<>%
  mutate(
    resp = resp + 1,
    correct_resp = correct_resp + 1,
    stim_type = stim_type + 1,
    condition = case_when(
      correct_resp == 2 & stim_type == 1 ~ 1,
      correct_resp == 2 & stim_type == 2 ~ 2,
      correct_resp == 1 & stim_type == 1 ~ 3,
      correct_resp == 1 & stim_type == 2 ~ 4
    )
  )

N <- length(unique(df$id))

stan_data = list(
  `T`        = nrow(df),
  N          = N,
  subject_id = df$id,
  resp       = df$resp,
  condition  = df$condition,
  rt         = df$rt
)

init_fun = function(chains=4){
  L = list()
  for (i in 1:chains) {
    L[[i]] = list(
      mu_v       = 1.0 + runif(4, -0.5, 0.5),
      sigma_v    = 0.4 + runif(4, -0.05, 0.05),
      mu_a       = 2.0 + runif(4, -0.2, 0.2),
      sigma_a    = 0.4 + runif(4, -0.05, 0.05),
      mu_bias    = 0.0 + runif(4, -0.05, 0.05),
      sigma_bias = 0.1 + runif(4, -0.05, 0.05),
      mu_ndt     = -1 + runif(4, -0.05, 0.05),
      sigma_ndt  = 0.1 + runif(4, -0.01, 0.01),
      v    = matrix(1.0 + runif(4 * N, -0.5, 0.5), 4, N),
      a    = matrix(2.0 + runif(4 * N, -0.1, 0.1), 4, N),
      bias = 0.5 + runif(4 * N, -0.05, 0.05),
      ndt  = matrix(0.2 + runif(4 * N, -0.01, 0.01), 4, N),
      ndt_s = 0.05 + runif(1, -0.01, 0.01)
    )
  }
  return(L)
}

m_full <- cmdstan_model(
  '../model/model_full_2_without_var.stan',
  cpp_options = list(stan_threads = T)
)

m_full_fit <- m_full$sample(
  data = stan_data,
  init = init_fun(),
  max_treedepth = 10,
  adapt_delta = 0.95,
  refresh = 10,
  iter_sampling = 2000,
  iter_warmup = 2000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE,
  save_cmdstan_config = TRUE
)

m_full_fit$save_object("../fits/exp2_full_model_without_var_2.rds")
