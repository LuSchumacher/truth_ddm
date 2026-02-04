library(tidyverse)
library(magrittr)
library(cmdstanr)
library(bayesplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

init_fun_bias <- function(chains = 4, n) {
  L <- vector("list", chains)
  for (i in 1:chains) {
    L[[i]] <- list(
      v_intercept.      = rnorm(1, 0, 0.1),
      v_beta            = rnorm(1, 0, 0.1),
      sigma_v           = runif(1, 0.2, 0.5),
      a_intercept       = runif(1, 1.5, 2.5),
      a_beta            = rnorm(1, 0, 0.1),
      sigma_a           = runif(1, 0.2, 0.5),
      bias_intercept    = rnorm(1, 0, 0.2),
      bias_betas        = rnorm(3, 0, 0.2),
      sigma_bias        = runif(1, 0.2, 0.5),
      ndt_intercept     = rnorm(1, -2, 0.05),
      ndt_beta          = rnorm(1, 0, 0.05),
      sigma_ndt         = abs(rnorm(1, 0, 0.05)),
      ndt_var_intercept = rnorm(1, 0, 0.05),
      ndt_var_beta      = rnorm(1, 0, 0.05),
      sigma_ndt_var     = abs(rnorm(1, 0, 0.05)),
      z_v               = rnorm(n, 0, 0.5),
      z_a               = rnorm(n, 0, 0.5),
      z_bias            = rnorm(n, 0, 0.5),
      z_ndt             = rnorm(n, 0, 0.5),
      z_ndt_var         = rnorm(n, 0, 0.5),
      s                 = runif(T, 0, 0.03),
      trel              = matrix(runif(2*n, 0.01, 0.8), nrow=2)
    )
  }
  return(L)
}

init_fun_drift <- function(chains = 4, n) {
  L <- vector("list", chains)
  for (i in 1:chains) {
    L[[i]] <- list(
      v_intercept       = rnorm(1, 0, 0.1),
      v_betas           = rnorm(3, 0, 0.1),
      sigma_v           = runif(1, 0.2, 0.5),
      a_intercept       = runif(1, 1.5, 2.5),
      a_beta            = rnorm(1, 0, 0.1),
      sigma_a           = runif(1, 0.2, 0.5),
      bias_intercept    = rnorm(1, 0, 0.2),
      bias_beta         = rnorm(1, 0, 0.2),
      sigma_bias        = runif(1, 0.2, 0.5),
      ndt_intercept     = rnorm(1, -2, 0.05),
      ndt_beta          = rnorm(1, 0, 0.05),
      sigma_ndt         = abs(rnorm(1, 0, 0.05)),
      ndt_var_intercept = rnorm(1, 0, 0.05),
      ndt_var_beta      = rnorm(1, 0, 0.05),
      sigma_ndt_var     = abs(rnorm(1, 0, 0.05)),
      z_v               = rnorm(n, 0, 0.5),
      z_a               = rnorm(n, 0, 0.5),
      z_bias            = rnorm(n, 0, 0.5),
      z_ndt             = rnorm(n, 0, 0.5),      
      z_ndt_var         = rnorm(n, 0, 0.5),
      s                 = runif(T, 0, 0.03),
      trel              = matrix(runif(2*n, 0.01, 0.8), nrow=2)
    )
  }
  return(L)
}

PARAM_NAMES_BIAS <- c(
  "v_intercept", "v_beta",
  "a_intercept", "a_beta",
  "bias_intercept", "bias_betas[1]", "bias_betas[2]", "bias_betas[3]",
  "ndt_intercept", "ndt_beta",
  "ndt_var_intercept", "ndt_var_beta"
)

PARAM_NAMES_DRIFT <- c(
  "v_intercept", "v_betas[1]", "v_betas[2]", "v_betas[3]",
  "a_intercept", "a_beta",
  "bias_intercept", "bias_beta",
  "ndt_intercept", "ndt_beta",
  "ndt_var_intercept", "ndt_var_beta"
)

truth_ddm_bias <- cmdstan_model(
  '../model/truth_ddm_restricted_bias.stan',
  cpp_options = list(stan_threads = T),
  force_recompile = TRUE
)

truth_ddm_drift <- cmdstan_model(
  '../model/truth_ddm_restricted_drift.stan',
  cpp_options = list(stan_threads = T),
  force_recompile = TRUE
)

# ---------------------------------------------------------------------------- #
# FITTING EXPERIMENT 1: SESSION 1
# ---------------------------------------------------------------------------- #
df_session_1 <- read.csv('../data/data_session_1.csv')
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == 1, -1, 1)
df_session_1$stim_type <- ifelse(df_session_1$stim_type == 0, -1, 1)
interaction_term <- df_session_1$factual_truth * df_session_1$stim_type
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == -1, 1, 2)

N <- length(unique(df_session_1$id))
`T` <- nrow(df_session_1)

stan_data <- list(
  `T`           = `T`,
  N             = N,
  subject_id    = df_session_1$id,
  resp          = df_session_1$resp,
  truth         = df_session_1$factual_truth,
  condition     = df_session_1$condition,
  rt            = df_session_1$rt,
  minRT         = tapply(
    df_session_1$rt,
    list(df_session_1$factual_truth, df_session_1$id),
    min
  )
)

fit_bias <- truth_ddm_bias$sample(
  data = stan_data,
  init = init_fun_bias(n=N),
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

print(fit_bias$summary(variables = PARAM_NAMES_BIAS), n=100)
mcmc_trace(fit_bias$draws(), pars=PARAM_NAMES_BIAS, facet_args = list(ncol = 4))
fit_bias$save_object("../fits/fit_restricted_bias.rds")


fit_drift <- truth_ddm_drift$sample(
  data = stan_data,
  init = init_fun_drift(n=N),
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

print(fit_drift$summary(variables = PARAM_NAMES_DRIFT), n=100)
mcmc_trace(fit_drift$draws(), pars=PARAM_NAMES_DRIFT, facet_args = list(ncol = 4))
fit_drift$save_object("../fits/fit_restricted_drift.rds")

# ---------------------------------------------------------------------------- #
# POSTERIOR RESIMULATON
# ---------------------------------------------------------------------------- #
fit_bias <- readRDS("../fits/fit_restricted_bias.rds")
fit_drift <- readRDS("../fits/fit_restricted_drift.rds")
df_session_1 <- read.csv('../data/data_session_1.csv')
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == 1, -1, 1)
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == -1, 1, 2)

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
SMALLER_FONT <- 4
COLOR_PALETTE <- c('#27374D', '#7A8F7A', '#B70404')

NUM_RESIMS <- 100
NUM_POST_SAMPLES <- 12000
idx <- sample(1:NUM_POST_SAMPLES, NUM_RESIMS, replace = FALSE)

source("0_ddm_simulator.R")

draws_fit_bias <- as_draws_matrix(fit_bias)
draws_fit_drift <- as_draws_matrix(fit_drift)

resimulate_bias <- function(df, params, num_resims = 100) {
  n_trials <- nrow(df)
  n_id <- length(unique(df$id))
  n_draws <- dim(params$v)[3]
  resim_id_vec <- seq_len(num_resims)
  idx_draws <- sample(seq_len(n_draws), num_resims, replace = FALSE)
  pred_data <- vector("list", num_resims)
  for (r_idx in seq_along(idx_draws)) {
    r <- idx_draws[r_idx]
    tmp_list <- vector("list", n_trials)
    for (t in seq_len(n_trials)) {
      id <- df$id[t]
      cond <- df$condition[t]
      truth <- df$factual_truth[t]
      v    <- params$v[truth, id, r]
      a    <- params$a[truth, id, r]
      bias <- params$bias[cond, id, r]
      ndt_base <- params$ndt[truth, id, r]
      ndt_var  <- params$ndt_var[truth, id, r]
      s_val <- params$s[t, r]
      ndt_current <- ndt_base + s_val * ndt_var
      x <- sample_ddm(v = v, a = a, ndt = ndt_current, bias = bias)
      tmp_list[[t]] <- tibble(
        resp = x[1],
        rt = x[2],
        id = id,
        condition = cond,
        resim_id = r_idx,
        stim_type = df$stim_type[t],
        factual_truth = df$factual_truth[t],
        correct = df$correct[t]
      )
    }
    pred_data[[r_idx]] <- bind_rows(tmp_list)
    cat("Resimulation", r_idx, "of", num_resims, "done\n")
  }
  bind_rows(pred_data)
}

resimulate_drift <- function(df, params, num_resims = 100) {
  n_trials <- nrow(df)
  n_id <- length(unique(df$id))
  n_draws <- dim(params$v)[3]
  resim_id_vec <- seq_len(num_resims)
  idx_draws <- sample(seq_len(n_draws), num_resims, replace = FALSE)
  pred_data <- vector("list", num_resims)
  for (r_idx in seq_along(idx_draws)) {
    r <- idx_draws[r_idx]
    tmp_list <- vector("list", n_trials)
    for (t in seq_len(n_trials)) {
      id <- df$id[t]
      cond <- df$condition[t]
      truth <- df$factual_truth[t]
      v    <- params$v[cond, id, r]
      a    <- params$a[truth, id, r]
      bias <- params$bias[truth, id, r]
      ndt_base <- params$ndt[truth, id, r]
      ndt_var  <- params$ndt_var[truth, id, r]
      s_val <- params$s[t, r]
      ndt_current <- ndt_base + s_val * ndt_var
      x <- sample_ddm(v = v, a = a, ndt = ndt_current, bias = bias)
      tmp_list[[t]] <- tibble(
        resp = x[1],
        rt = x[2],
        id = id,
        condition = cond,
        resim_id = r_idx,
        stim_type = df$stim_type[t],
        factual_truth = df$factual_truth[t],
        correct = df$correct[t]
      )
    }
    pred_data[[r_idx]] <- bind_rows(tmp_list)
    cat("Resimulation", r_idx, "of", num_resims, "done\n")
  }
  bind_rows(pred_data)
}

extract_param_3d_subset <- function(fit, name, n_id, idx, n_cond = 4) {
  dm <- fit$draws(variables = name, format = "draws_matrix") |> as.matrix()
  dm <- dm[idx, , drop = FALSE]
  cols <- colnames(dm)
  m <- stringr::str_match(cols, paste0("^", name, "\\[(\\d+),(\\d+)\\]$"))
  i <- as.integer(m[,2]); j <- as.integer(m[,3])
  n_draw <- length(idx)
  arr <- array(NA_real_, dim = c(n_cond, n_id, n_draw))
  for (k in seq_along(cols)) {
    arr[i[k], j[k], ] <- dm[, k]
  }
  arr
}

extract_param_trialwise_subset <- function(fit, name, idx) {
  dm <- fit$draws(variables = name, format = "draws_matrix") |> as.matrix()
  dm <- dm[idx, , drop = FALSE]
  cols <- colnames(dm)
  m <- stringr::str_match(cols, paste0("^", name, "\\[(\\d+)\\]$"))
  t_idx <- as.integer(m[,2])
  n_draw <- length(idx)
  n_trial <- max(t_idx)
  out <- matrix(NA_real_, nrow = n_trial, ncol = n_draw)
  for (k in seq_along(cols)) {
    out[t_idx[k], ] <- dm[, k]
  }
  out
}

extract_all_params <- function(fit, df, idx, n_cond = 4) {
  n_id <- length(unique(df$id))
  list(
    v       = extract_param_3d_subset(fit, "v",       n_id, idx, n_cond),
    a       = extract_param_3d_subset(fit, "a",       n_id, idx, n_cond),
    bias    = extract_param_3d_subset(fit, "bias",    n_id, idx, n_cond),
    ndt     = extract_param_3d_subset(fit, "ndt",     n_id, idx, n_cond),
    ndt_var = extract_param_3d_subset(fit, "ndt_var", n_id, idx, n_cond),
    s       = extract_param_trialwise_subset(fit, "s", idx)
  )
}

num_id_session_1 <- length(unique(df_session_1$id))
params_fit_bias <- extract_all_params(fit_bias, df_session_1, idx)
pred_data_fit_bias <- resimulate_bias(df_session_1, params_fit_bias)

params_fit_drift <- extract_all_params(fit_drift, df_session_1, idx)
pred_data_fit_drift <- resimulate_drift(df_session_1, params_fit_drift)


# ---------------------------------------------------------------------------- #
# RT QUANTILES
# ---------------------------------------------------------------------------- #
quantiles <- seq(0.1, 0.9, 0.1)

emp_rt_summaries <- df_session_1 %>%
  group_by(stim_type) %>%
  summarise(
    data = list({
      qs <- quantile(rt, probs = quantiles, na.rm = TRUE)
      tibble(
        quantile = quantiles,
        rt = qs
      )
    }),
    .groups = "drop"
  ) %>%
  unnest(data) %>%
  mutate(
    data_type = "Observed",
    quantile = factor(quantile)
  )

summarise_pred_quantiles <- function(df, label) {
  df %>%
    group_by(resim_id, stim_type) %>%
    summarise(
      data = list({
        qs <- quantile(rt, probs = quantiles)
        tibble(
          quantile = quantiles,
          rt = qs
        )
      }),
      .groups = "drop"
    ) %>%
    unnest(data) %>%
    group_by(stim_type, quantile) %>%
    summarise(
      mean = mean(rt),
      q_lower = quantile(rt, 0.025),
      q_upper = quantile(rt, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      data_type = label,
      quantile = factor(quantile)
    )
}

pred_rt_bias  <- summarise_pred_quantiles(
  pred_data_fit_bias,
  "Only bias modulation"
)

pred_rt_drift <- summarise_pred_quantiles(
  pred_data_fit_drift,
  "Only drift modulation"
)

plot_pred <- dplyr::bind_rows(
  dplyr::mutate(pred_rt_bias, source = "Only bias modulation"),
  dplyr::mutate(pred_rt_drift, source = "Only drift modulation")
)
plot_obs <- dplyr::mutate(emp_rt_summaries, source = "Observed")

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #
pd <- position_dodge(width = 0.8)

ggplot(
  plot_pred,
  aes(
    x = quantile,
    y = mean,
    ymin = q_lower,
    ymax = q_upper,
    color = source,
    group = source
  )
) +
  geom_pointrange(
    position = pd,
    linewidth = 1.2, size = 0.5
  ) +
  geom_point(
    data = plot_obs,
    inherit.aes = FALSE,
    aes(
      x = quantile,
      y = rt,
      color = source,
      group = source
    ),
    position = pd,
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "Only bias modulation" = COLOR_PALETTE[1],
      "Only drift modulation" = COLOR_PALETTE[2],
      "Observed" = COLOR_PALETTE[3]
    )
  ) +
  facet_wrap(
    ~ stim_type,
    nrow = 1,
    labeller = labeller(
      stim_type = c(
        "0" = "New\nstatements",
        "1" = "Repeated\nstatements"
      )
    )
  ) +
  labs(
    x = "Quantile",
    y = "Response time (s)",
    color = ""
  ) +
  ggthemes::theme_tufte(base_size = FONT_SIZE_2 - SMALLER_FONT) +
  theme(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.line = element_line(linewidth = 0.5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3 - SMALLER_FONT, vjust = 0.5),
    axis.text.y = element_text(size = FONT_SIZE_3 - SMALLER_FONT),
    strip.text.x = element_text(size = FONT_SIZE_2 - SMALLER_FONT, angle = 0),
    strip.text.y = element_text(size = FONT_SIZE_2 - SMALLER_FONT, angle = 0),
    panel.grid.major = element_line(color = scales::alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = scales::alpha("gray70", 0.15)),
    panel.background = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(0.2, "cm")
  )

ggsave(
  '../plots/07_rt_quantiles_restricted_models.jpeg',
  device = 'jpeg', dpi = 300,
  width = 10, height = 6
)

# ---------------------------------------------------------------------------- #
# CUMMULATIVE ACCURACY X RT QUANTILES
# ---------------------------------------------------------------------------- #
quantiles <- seq(0.1, 0.9, 0.1)

caf_obs <- df_session_1 %>%
  group_by(stim_type) %>%
  summarise(
    data = list({
      qs <- quantile(rt, probs = quantiles, na.rm = TRUE)
      tibble(
        quantile = quantiles,
        rt = qs,
        cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x], na.rm = TRUE))
      )
    }),
    .groups = "drop"
  ) %>%
  unnest(data) %>%
  mutate(
    quantile = as_factor(quantile),
    source = "Observed"
  )

caf_pred_bias <- pred_data_fit_bias %>%
  group_by(resim_id, stim_type) %>%
  summarise(
    data = list({
      qs <- quantile(rt, probs = quantiles)
      tibble(
        quantile = quantiles,
        rt = qs,
        cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x]))
      )
    }),
    .groups = "drop"
  ) %>%
  unnest(data) %>%
  group_by(stim_type, quantile) %>%
  summarise(
    rt = mean(rt),
    mean = mean(cum_acc),
    q_lower = quantile(cum_acc, 0.025),
    q_upper = quantile(cum_acc, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    quantile = as_factor(quantile),
    source = "Only bias modulation"
  )

caf_pred_drift <- pred_data_fit_drift %>%
  group_by(resim_id, stim_type) %>%
  summarise(
    data = list({
      qs <- quantile(rt, probs = quantiles)
      tibble(
        quantile = quantiles,
        rt = qs,
        cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x]))
      )
    }),
    .groups = "drop"
  ) %>%
  unnest(data) %>%
  group_by(stim_type, quantile) %>%
  summarise(
    rt = mean(rt),
    mean = mean(cum_acc),
    q_lower = quantile(cum_acc, 0.025),
    q_upper = quantile(cum_acc, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    quantile = as_factor(quantile),
    source = "Only drift modulation"
  )

ggplot(
  caf_obs,
  aes(
    x = rt,
    y = cum_acc,
    group = stim_type
  )
) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.6,
    aes(color = "Observed")
  ) +
  geom_point(
    shape = 4,
    size = 3,
    stroke = 1.3,
    aes(color = "Observed")
  ) +
  geom_line(
    data = caf_pred_bias,
    aes(
      x = rt,
      y = mean,
      group = stim_type,
      color = "Only bias modulation"
    ),
    inherit.aes = FALSE,
    linewidth = 1.0,
    linetype = "dashed"
  ) +
  geom_point(
    data = caf_pred_bias,
    aes(
      x = rt,
      y = mean,
      color = "Only bias modulation"
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = 1.3
  ) +
  geom_line(
    data = caf_pred_drift,
    aes(
      x = rt,
      y = mean,
      group = stim_type,
      color = "Only drift modulation"
    ),
    inherit.aes = FALSE,
    linewidth = 1.0,
    linetype = "dashed"
  ) +
  geom_point(
    data = caf_pred_drift,
    aes(
      x = rt,
      y = mean,
      color = "Only drift modulation"
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3,
    stroke = 1.3
  ) +
  facet_grid(
    stim_type ~ .,
    labeller = labeller(
      stim_type = c(
        "0" = "New\nstatements",
        "1" = "Repeated\nstatements"
      )
    )
  ) +
  scale_color_manual(
    name = "",
    breaks = c(
      "Observed",
      "Only bias modulation",
      "Only drift modulation"
    ),
    values = c(
      "Only bias modulation" = COLOR_PALETTE[1],
      "Only drift modulation" = COLOR_PALETTE[2],
      "Observed" = COLOR_PALETTE[3]
    )
  ) +
  labs(
    x = "Response time (s)",
    y = "Cumulative accuracy"
  ) +
  ggthemes::theme_tufte(base_size = FONT_SIZE_2 - SMALLER_FONT) +
  theme(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.line = element_line(linewidth = 0.5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3 - SMALLER_FONT),
    axis.text.y = element_text(size = FONT_SIZE_3 - SMALLER_FONT),
    strip.text.x = element_text(size = FONT_SIZE_2 - SMALLER_FONT),
    strip.text.y = element_text(size = FONT_SIZE_2 - SMALLER_FONT),
    panel.grid.major = element_line(color = scales::alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = scales::alpha("gray70", 0.15)),
    panel.background = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )
