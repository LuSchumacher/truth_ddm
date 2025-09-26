library(tidyverse)
library(magrittr)
library(cmdstanr)
library(bayesplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

init_fun <- function(chains = 4, n) {
  L <- vector("list", chains)
  for (i in 1:chains) {
    L[[i]] <- list(
      # --- v (drift rate) ---
      v_intercept = rnorm(1, 0, 0.1),
      v_betas     = rnorm(3, 0, 0.1),
      sigma_v        = runif(1, 0.2, 0.5),
      # --- a (boundary separation) ---
      a_intercept = runif(1, 1.5, 2.5),
      a_betas     = rnorm(3, 0, 0.1),
      sigma_a        = runif(1, 0.2, 0.5),
      # --- bias ---
      bias_intercept = rnorm(1, 0, 0.2),
      bias_betas     = rnorm(3, 0, 0.2),
      sigma_bias        = runif(1, 0.2, 0.5),
      # --- ndt ---
      ndt_intercept = rnorm(1, -2, 0.05),
      ndt_betas     = rnorm(3, 0, 0.05),
      sigma_ndt        = abs(rnorm(1, 0, 0.05)),
      # --- ndt_var ---
      ndt_var_intercept = rnorm(1, 0, 0.05),
      ndt_var_betas     = rnorm(3, 0, 0.05),
      sigma_ndt_var        = abs(rnorm(1, 0, 0.05)),
      # --- z (subject-level random effects, std_normal) ---
      z_v        = rnorm(n, 0, 0.5),
      z_a        = rnorm(n, 0, 0.5),
      z_bias     = rnorm(n, 0, 0.5),
      z_ndt      = rnorm(n, 0, 0.5),
      z_ndt_var  = rnorm(n, 0, 0.5),
      # --- trial-level t0 variability (if modeled) ---
      s    = runif(T, 0, 0.03),
      trel = matrix(runif(4*n, 0.01, 0.8), nrow=4)
    )
  }
  return(L)
}

PARAM_NAMES <- c(
  "v_intercept", "v_betas[1]", "v_betas[2]", "v_betas[3]",
  "a_intercept", "a_betas[1]", "a_betas[2]", "a_betas[3]",
  "bias_intercept", "bias_betas[1]", "bias_betas[2]", "bias_betas[3]",
  "ndt_intercept", "ndt_betas[1]", "ndt_betas[2]", "ndt_betas[3]",
  "ndt_var_intercept", "ndt_var_betas[1]", "ndt_var_betas[2]", "ndt_var_betas[3]"
)

truth_ddm <- cmdstan_model(
  '../model/truth_ddm.stan',
  cpp_options = list(stan_threads = T)
)

################################################################################
# EXPERIMENT 1: SESSION 1
################################################################################
df_session_1 <- read.csv('../data/data_session_1.csv')
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == 0, -1, 1)
df_session_1$stim_type <- ifelse(df_session_1$stim_type == 0, -1, 1)
interaction_term <- df_session_1$factual_truth * df_session_1$stim_type

N <- length(unique(df_session_1$id))
`T` <- nrow(df_session_1)

# Prepare Stan data list
stan_data <- list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_session_1$id,
  resp          = df_session_1$resp,
  truth         = df_session_1$factual_truth,
  repetition    = df_session_1$stim_type,
  interaction   = interaction_term,
  condition     = df_session_1$condition,
  rt            = df_session_1$rt,
  minRT         = tapply(df_session_1$rt, list(df_session_1$condition, df_session_1$id), min)
)

fit_session_1 <- truth_ddm$sample(
  data = stan_data,
  init = init_fun(n=N),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 25,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

print(fit_session_1$summary(variables = PARAM_NAMES), n=100)
mcmc_trace(fit_session_1$draws(), pars=PARAM_NAMES, facet_args = list(ncol = 4))
fit_session_1$save_object("../fits/fit_session_1.rds")

################################################################################
# EXPERIMENT 1: SESSION 2
################################################################################
df_session_2 <- read.csv('../data/data_session_2.csv')
df_session_2$factual_truth <- ifelse(df_session_2$factual_truth == 0, -1, 1)
df_session_2$stim_type <- ifelse(df_session_2$stim_type == 0, -1, 1)
interaction_term <- df_session_2$factual_truth * df_session_2$stim_type

N <- length(unique(df_session_2$id))
`T` <- nrow(df_session_2)

# Prepare Stan data list
stan_data <- list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_session_2$id,
  resp          = df_session_2$resp,
  truth         = df_session_2$factual_truth,
  repetition    = df_session_2$stim_type,
  interaction   = interaction_term,
  condition     = df_session_2$condition,
  rt            = df_session_2$rt,
  minRT         = tapply(df_session_2$rt, list(df_session_2$condition, df_session_2$id), min)
)

fit_session_2 <- truth_ddm$sample(
  data = stan_data,
  init = init_fun(n=N),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 25,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

print(fit_session_2$summary(variables = PARAM_NAMES), n=100)
mcmc_trace(fit_session_2$draws(), pars=PARAM_NAMES, facet_args = list(ncol = 4))
fit_session_2$save_object("../fits/fit_session_2.rds")

################################################################################
# EXPERIMENT 2
################################################################################
df_exp_2 <- read.csv('../data/data_exp_2.csv')
df_exp_2$factual_truth <- ifelse(df_exp_2$factual_truth == 0, -1, 1)
df_exp_2$stim_type <- ifelse(df_exp_2$stim_type == 0, -1, 1)
interaction_term <- df_exp_2$factual_truth * df_exp_2$stim_type

N <- length(unique(df_exp_2$id))
`T` <- nrow(df_exp_2)

# Prepare Stan data list
stan_data <- list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_exp_2$id,
  resp          = df_exp_2$resp,
  truth         = df_exp_2$factual_truth,
  repetition    = df_exp_2$stim_type,
  interaction   = interaction_term,
  condition     = df_exp_2$condition,
  rt            = df_exp_2$rt,
  minRT         = tapply(df_exp_2$rt, list(df_exp_2$condition, df_exp_2$id), min)
)

fit_exp_2 <- truth_ddm$sample(
  data = stan_data,
  init = init_fun(n=N),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 25,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

print(fit_exp_2$summary(variables = PARAM_NAMES), n=100)
mcmc_trace(fit_exp_2$draws(), pars=PARAM_NAMES, facet_args = list(ncol = 4))
fit_exp_2$save_object("../fits/fit_exp_2.rds")

################################################################################
# EXPERIMENT 2: COMBINED RT (EXPLORATIVE)
################################################################################
df_exp_2 <- read.csv('../data/data_exp_2.csv')
df_exp_2$factual_truth <- ifelse(df_exp_2$factual_truth == 0, -1, 1)
df_exp_2$stim_type <- ifelse(df_exp_2$stim_type == 0, -1, 1)
interaction_term <- df_exp_2$factual_truth * df_exp_2$stim_type

N <- length(unique(df_exp_2$id))
`T` <- nrow(df_exp_2)

# Prepare Stan data list
stan_data <- list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_exp_2$id,
  resp          = df_exp_2$resp,
  truth         = df_exp_2$factual_truth,
  repetition    = df_exp_2$stim_type,
  interaction   = interaction_term,
  condition     = df_exp_2$condition,
  rt            = df_exp_2$rt_total,
  minRT         = tapply(df_exp_2$rt, list(df_exp_2$condition, df_exp_2$id), min)
)

fit_exp_2_rt_total <- truth_ddm$sample(
  data = stan_data,
  init = init_fun(n=N),
  max_treedepth = 12,
  adapt_delta = 0.95,
  refresh = 25,
  iter_sampling = 3000,
  iter_warmup = 3000,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  save_warmup = TRUE
)

print(fit_exp_2_rt_total$summary(variables = PARAM_NAMES), n=100)
mcmc_trace(fit_exp_2_rt_total$draws(), pars=PARAM_NAMES, facet_args = list(ncol = 4))
fit_exp_2_rt_total$save_object("../fits/fit_exp_2_rt_total.rds")
