library(tidyverse)
library(magrittr)
library(LaplacesDemon)
library(cmdstanr)
library(patchwork)
library(posterior)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("0_ddm_simulator.R")

softplus <- function(x) {
  y <- ifelse(x > 20, x + log1p(exp(-x)), log1p(exp(x)))
  return(y)
}

df <- read_csv('../data/data_session_1.csv') 
fit_session_1 <- readRDS("../fits/fit_session_1.rds")

NUM_SIM <- 50
NUM_SUBS <- length(unique(df$id))

NUM_CONDS <- 4

GROUP_PARAM_NAMES <- c(
  "transf_mu_v[1]", "transf_mu_v[2]",
  "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]",
  "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]",
  "transf_mu_bias[3]", "transf_mu_bias[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]",
  "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_ndt_var[1]", "transf_mu_ndt_var[2]",
  "transf_mu_ndt_var[3]", "transf_mu_ndt_var[4]"
)

PARAM_LABELS <- c(
  "v"       = "italic(v)",
  "a"       = "italic(a)",
  "bias"    = "beta",
  "ndt"     = "tau",
  "ndt_var" = "tau[var]"
)

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18

set.seed(1991)

# ---------------------------------------------------------------------------- #
# PARAMETER SAMPLING
# ---------------------------------------------------------------------------- #
group_level_draws <- fit_session_1$draws(
  format = "draws_df",
  variables = GROUP_PARAM_NAMES
)

group_df <- as_draws_df(group_level_draws)
print(colnames(group_df))
num_chains <- 4
num_iter_per_chain <- 3000

group_df <- group_df %>%
  mutate(draw_across = (.chain - 1) * num_iter_per_chain + .iteration)

group_long <- group_df %>%
  select(-.iteration, -.draw) %>%
  pivot_longer(
    cols = all_of(GROUP_PARAM_NAMES),
    names_to = c("parameter", "condition"),
    names_pattern = "(.*)\\[(\\d+)\\]",
    values_to = "value"
  ) %>%
  mutate(condition = as.integer(condition),
         parameter = str_remove(parameter, "^transf_")) %>%
  select(draw = draw_across, parameter, condition, value)

draws <- fit_session_1$draws(format = "draws_df")

v_intercept <- draws %>% select(v_intercept) %>% pull()
v_betas <- draws %>% select(starts_with("v_betas")) %>% as.matrix()
sigma_v <- draws %>% select(sigma_v) %>% pull()
z_v <- draws %>% select(starts_with("z_v[")) %>% as.matrix()

a_intercept <- draws %>% select(a_intercept) %>% pull()
a_betas <- draws %>% select(starts_with("a_betas")) %>% as.matrix()
sigma_a <- draws %>% select(sigma_a) %>% pull()
z_a <- draws %>% select(starts_with("z_a[")) %>% as.matrix()

bias_intercept <- draws %>% select(bias_intercept) %>% pull()
bias_betas <- draws %>% select(starts_with("bias_betas")) %>% as.matrix()
sigma_bias <- draws %>% select(sigma_bias) %>% pull()
z_bias <- draws %>% select(starts_with("z_bias[")) %>% as.matrix()

ndt_intercept <- draws %>% select(ndt_intercept) %>% pull()
ndt_betas <- draws %>% select(starts_with("ndt_betas")) %>% as.matrix()
sigma_ndt <- draws %>% select(sigma_ndt) %>% pull()
z_ndt <- draws %>% select(starts_with("z_ndt[")) %>% as.matrix()

ndt_var_intercept <- draws %>% select(ndt_var_intercept) %>% pull()
ndt_var_betas <- draws %>% select(starts_with("ndt_var_betas")) %>% as.matrix()
sigma_ndt_var <- draws %>% select(sigma_ndt_var) %>% pull()
z_ndt_var <- draws %>% select(starts_with("z_ndt_var[")) %>% as.matrix()

num_draws <- length(v_intercept)

x_cells <- matrix(c(
  -1,  1, -1,  1,  # repetition
  1,  1, -1, -1,  # truth
  -1,  1,  1, -1   # interaction
), nrow=3, byrow=TRUE)

v_samples <- array(NA, dim = c(num_draws, NUM_CONDS, NUM_SUBS))
a_samples <- array(NA, dim = c(num_draws, NUM_CONDS, NUM_SUBS))
bias_samples <- array(NA, dim = c(num_draws, NUM_CONDS, NUM_SUBS))
ndt_samples <- array(NA, dim = c(num_draws, NUM_CONDS, NUM_SUBS))
ndt_var_samples <- array(NA, dim = c(num_draws, NUM_CONDS, NUM_SUBS))

for (draw_i in seq_len(num_draws)) {
  # Calculate condition means mu for each param
  mu_v <- v_intercept[draw_i] + as.numeric(v_betas[draw_i, ]) %*% x_cells
  mu_a <- a_intercept[draw_i] + as.numeric(a_betas[draw_i, ]) %*% x_cells
  mu_bias <- bias_intercept[draw_i] + as.numeric(bias_betas[draw_i, ]) %*% x_cells
  mu_ndt <- ndt_intercept[draw_i] + as.numeric(ndt_betas[draw_i, ]) %*% x_cells
  mu_ndt_var <- ndt_var_intercept[draw_i] + as.numeric(ndt_var_betas[draw_i, ]) %*% x_cells
  
  for (cond_i in seq_len(NUM_CONDS)) {
    for (subj_i in seq_len(NUM_SUBS)) {
      # Extract subject-level random effects for this draw and subject
      z_v_s <- z_v[draw_i, subj_i]
      z_a_s <- z_a[draw_i, subj_i]
      z_bias_s <- z_bias[draw_i, subj_i]
      z_ndt_s <- z_ndt[draw_i, subj_i]
      z_ndt_var_s <- z_ndt_var[draw_i, subj_i]
      
      # Apply transformations same as Stan model
      s_v <- softplus(sigma_v[draw_i])
      s_a <- softplus(sigma_a[draw_i])
      s_bias <- softplus(sigma_bias[draw_i])
      s_ndt <- softplus(sigma_ndt[draw_i])
      s_ndt_var <- softplus(sigma_ndt_var[draw_i])
      
      v_samples[draw_i, cond_i, subj_i] <- mu_v[cond_i] + s_v * z_v_s
      a_samples[draw_i, cond_i, subj_i] <- softplus(mu_a[cond_i] + s_a * z_a_s)
      bias_samples[draw_i, cond_i, subj_i] <- plogis(mu_bias[cond_i] + s_bias * z_bias_s)
      ndt_samples[draw_i, cond_i, subj_i] <- mu_ndt[cond_i] + s_ndt * z_ndt_s
      ndt_var_samples[draw_i, cond_i, subj_i] <- softplus(mu_ndt_var[cond_i] + s_ndt_var * z_ndt_var_s)
    }
  }
}

num_chains <- 4
num_iter_per_chain <- 3000

# Add draw_across as absolute draw index
draws <- draws %>%
  mutate(draw_across = (.chain - 1) * num_iter_per_chain + .iteration)

# Helper function to convert 3D array [draw, condition, subject] to long df
array_to_long_df <- function(arr, param_name) {
  dims <- dim(arr)
  df <- as.data.frame.table(arr, responseName = "value") %>%
    rename(draw = 1, condition = 2, subject = 3) %>%
    mutate(
      draw = as.integer(draw),
      condition = as.integer(condition),
      subject = as.integer(subject),
      parameter = param_name
    )
  return(df)
}

# Convert all parameter arrays to long format
v_long       <- array_to_long_df(v_samples, "v")
a_long       <- array_to_long_df(a_samples, "a")
bias_long    <- array_to_long_df(bias_samples, "bias")
ndt_long     <- array_to_long_df(ndt_samples, "ndt")
ndt_var_long <- array_to_long_df(ndt_var_samples, "ndt_var")

# Combine all into one dataframe
individuals_long <- bind_rows(v_long, a_long, bias_long, ndt_long, ndt_var_long)

# Join chain and iteration info from draws
draw_info <- draws %>% select(.chain, .iteration, draw_across)

individuals_long <- individuals_long %>%
  left_join(draw_info, by = c("draw" = "draw_across")) %>%
  select(
    draw,
    parameter,
    condition,
    subject,
    value
  )

sampled_draws <- sample(unique(draws$draw_across), NUM_SIM)
group_params <- group_long %>%
  filter(draw %in% sampled_draws) %>% 
  mutate(draw = dense_rank(draw)) %>% 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) %>%
  arrange(draw, condition) %>% 
  rename(
    sim_id = draw
  )
individual_params <- individuals_long %>%
  filter(draw %in% sampled_draws) %>% 
  mutate(draw = dense_rank(draw)) %>% 
  rename(
    sim_id = draw,
    id = subject
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) %>%
  arrange(sim_id, id, condition)

write_csv(
  group_params,
  "../generated_data/param_recovery_new/true_group_params.csv"
)
write_csv(
  individual_params,
  "../generated_data/param_recovery_new/true_individual_params.csv"
)

# ---------------------------------------------------------------------------- #
# DATA SIMULATION
# ---------------------------------------------------------------------------- #
df_sim <- df %>%
  mutate(trial = row_number()) %>%
  tidyr::crossing(sim_id = 1:NUM_SIM)

df_sim <- df_sim %>%
  left_join(
    individual_params,
    by = c("sim_id", "id", "condition")
  )

sim_data <- df_sim %>%
  rowwise() %>%
  mutate(
    t0  = ndt + runif(1) * ndt_var,
    sim = list(sample_ddm(v = v, a = a, ndt = t0, bias = bias)),
    resp = sim[1],
    rt   = sim[2]
  ) %>%
  ungroup() %>% 
  select(
    sim_id,
    id,
    condition,
    stim_type,
    factual_truth,
    resp,
    rt
  )

write_csv(sim_data, "../generated_data/param_recovery_new/sim_data_recovery.csv")

# ---------------------------------------------------------------------------- #
# MODEL FITTING
# ---------------------------------------------------------------------------- #
sim_data <- read_csv("../generated_data/param_recovery_new/sim_data_recovery.csv")

init_fun <- function(chains = 4, n) {
  L <- vector("list", chains)
  for (i in 1:chains) {
    L[[i]] <- list(
      v_intercept = rnorm(1, 0, 0.1),
      v_betas     = rnorm(3, 0, 0.1),
      sigma_v        = runif(1, 0.2, 0.5),
      a_intercept = runif(1, 1.5, 2.5),
      a_betas     = rnorm(3, 0, 0.1),
      sigma_a        = runif(1, 0.2, 0.5),
      bias_intercept = rnorm(1, 0, 0.2),
      bias_betas     = rnorm(3, 0, 0.2),
      sigma_bias        = runif(1, 0.2, 0.5),
      ndt_intercept = rnorm(1, -2, 0.05),
      ndt_betas     = rnorm(3, 0, 0.05),
      sigma_ndt        = abs(rnorm(1, 0, 0.05)),
      ndt_var_intercept = rnorm(1, 0, 0.05),
      ndt_var_betas     = rnorm(3, 0, 0.05),
      sigma_ndt_var        = abs(rnorm(1, 0, 0.05)),
      z_v        = rnorm(n, 0, 0.5),
      z_a        = rnorm(n, 0, 0.5),
      z_bias     = rnorm(n, 0, 0.5),
      z_ndt      = rnorm(n, 0, 0.5),
      z_ndt_var  = rnorm(n, 0, 0.5),
      s    = runif(T, 0, 0.03),
      trel = matrix(runif(4*n, 0.01, 0.8), nrow=4)
    )
  }
  return(L)
}

truth_ddm <- cmdstan_model(
  '../model/truth_ddm.stan',
  cpp_options = list(stan_threads = T),
  force_recompile = TRUE
)

for (i in 1:NUM_SIM) {
  fitting_data <- sim_data %>%
    filter(sim_id == i)
  fitting_data$factual_truth <- ifelse(fitting_data$factual_truth == 0, -1, 1)
  fitting_data$stim_type <- ifelse(fitting_data$stim_type == 0, -1, 1)
  interaction_term <- fitting_data$factual_truth * fitting_data$stim_type

  N <- length(unique(fitting_data$id))
  `T` <- nrow(fitting_data)

  stan_data <- list(
    `T`           = `T`,
    N             = N,
    subject_id    = fitting_data$id,
    resp          = fitting_data$resp,
    truth         = fitting_data$factual_truth,
    repetition    = fitting_data$stim_type,
    interaction   = interaction_term,
    condition     = fitting_data$condition,
    rt            = fitting_data$rt,
    minRT         = tapply(
      fitting_data$rt, list(fitting_data$condition, fitting_data$id), min
    )
  )

  fit_sim_data <- truth_ddm$sample(
    data = stan_data,
    init = init_fun(n=N),
    max_treedepth = 8,
    adapt_delta = 0.85,
    refresh = 100,
    iter_sampling = 1000,
    iter_warmup = 1000,
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 2
  )

  fit_sim_data$save_object(paste0("../fits/param_recovery_new/fit_sim_data_", i, ".rds"))

  group_param_estimates <- fit_sim_data$summary(variables = GROUP_PARAM_NAMES)
  group_param_estimates <- group_param_estimates %>%
    select(name = 1, median = 3) %>%
    mutate(
      parameter = str_extract(name, "(?<=transf_mu_)[a-z_]+"),
      condition = as.integer(str_extract(name, "\\d+")),
      sim_id = i
    ) %>%
    select(sim_id, condition, parameter, median) %>%
    arrange(sim_id, parameter, condition)
  write_csv(
    group_param_estimates,
    paste0("../generated_data/param_recovery_new/group_param_estimates_recovery_", i, ".csv")
  )
  
  draws <- fit_sim_data$draws()
  draws_df <- as_draws_df(draws) %>%
    select(starts_with("indiv_transf_"))
  
  individual_param_estimates <- draws_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "param",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_match(param, "indiv_transf_(.*)\\[")[,2],
      condition = as.integer(str_match(param, "\\[(\\d+),")[,2]),
      subject_id = as.integer(str_match(param, ",(\\d+)\\]")[,2])
    ) %>%
    group_by(parameter, condition, subject_id) %>%
    summarise(
      median = median(value),
      .groups = "drop"
    ) %>%
    mutate(sim_id = i) %>%
    select(
      sim_id,
      parameter,
      condition,
      subject_id,
      median
    )
  write_csv(
    individual_param_estimates,
    paste0("../generated_data/param_recovery_new/individual_param_estimates_recovery_", i, ".csv")
  )
  print(paste("Fitting Nr.", i, "finished"))
}

# ---------------------------------------------------------------------------- #
# EVALUATION: GROUP-LEVEL PARAMETERS
# ---------------------------------------------------------------------------- #
group_params <- read_csv(
  "../generated_data/param_recovery_new/true_group_params.csv"
)

recovery_params <- group_params %>%
  select(
    sim_id, condition, mu_v,
    mu_a, mu_bias, mu_ndt, mu_ndt_var
  ) %>%
  pivot_longer(
    cols = -c(sim_id, condition),
    names_to = "parameter",
    values_to = "value_true"
  ) %>%
  mutate(parameter = case_when(
    parameter == "mu_v" ~ "v",
    parameter == "mu_a" ~ "a",
    parameter == "mu_ndt" ~ "ndt",
    parameter == "mu_ndt_var" ~ "ndt_var",
    parameter == "mu_bias" ~ "bias",
  )) %>%
  arrange(sim_id, parameter, condition)

path <- "../generated_data/param_recovery_new/"
files <- list.files(path, pattern = "group_param_estimates_recovery")
pred_params <- read_csv(paste0(path, files)) %>% 
  arrange(sim_id, parameter, condition)

recovery_params <- recovery_params %>%
  mutate(
    value_pred = pred_params$median,
    param_label = recode(parameter, !!!PARAM_LABELS),
    param_label = factor(param_label),
    parameter = factor(parameter, c("v","a","bias","ndt","ndt_var"))
  )

corr_scores <- recovery_params %>%
  group_by(parameter) %>%
  summarise(
    r = cor(value_true, value_pred)
  )

recovery_params <- recovery_params %>%
  left_join(corr_scores, by = "parameter") %>% 
  mutate(
    param_label = PARAM_LABELS[parameter],
    corr_label = paste0("r = ", sprintf("%.2f", r))
  )

combined_plot <- ggplot(recovery_params, aes(x = value_true, y = value_pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
  facet_wrap(
    ~ parameter,
    scales = "free",
    ncol = 5,
    labeller = as_labeller(PARAM_LABELS, default = label_parsed)
  ) +
  geom_text(
    data = recovery_params %>% group_by(parameter) %>% summarize(
      x = min(value_true),
      y = max(value_pred),
      corr_label = unique(corr_label)
    ),
    mapping = aes(x = x, y = y, label = corr_label),
    hjust = 0, vjust = 1,
    size = 6,
    fontface = "italic",
    family = "Palatino",
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(x = "True", y = "Estimated") +
  scale_x_continuous(n.breaks = 4) +
  ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(linewidth = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_2),
    strip.text.y = element_text(size = FONT_SIZE_2, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    plot.title = element_text(size = FONT_SIZE_1,
                              hjust = 0.5,
                              face = 'bold'),
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank(),
    legend.spacing.y = unit(0.25, 'cm'),
    panel.spacing = unit(1., "lines")
  )

ggsave(
  "../plots/06_parameter_recovery_new_plot.jpeg",
  plot = combined_plot,
  width = 14, height = 3.5, dpi = 300,
  device = "jpeg"
)

# ---------------------------------------------------------------------------- #
# EVALUATION: INDIVIDUAL PARAMETERS
# ---------------------------------------------------------------------------- #
individual_params <- read_csv(
  "../generated_data/param_recovery_new/true_individual_params.csv"
)

recovery_params <- individual_params %>%
  pivot_longer(
    cols = -c(sim_id, id, condition),
    names_to = "parameter",
    values_to = "value_true"
  ) %>%
  arrange(sim_id, id, parameter, condition)

path <- "../generated_data/param_recovery_new/"
files <- list.files(path, pattern = "individual_param_estimates_recovery")
pred_params <- read_csv(paste0(path, files)) %>% 
  arrange(sim_id, subject_id, parameter, condition)

recovery_params <- recovery_params %>%
  mutate(
    value_pred = pred_params$median,
    param_label = recode(parameter, !!!PARAM_LABELS),
    param_label = factor(param_label),
    parameter = factor(parameter, c("v","a","bias","ndt","ndt_var"))
  )

corr_scores <- recovery_params %>%
  group_by(parameter) %>%
  summarise(
    r = cor(value_true, value_pred)
  )

recovery_params <- recovery_params %>%
  left_join(corr_scores, by = "parameter") %>% 
  mutate(
    param_label = PARAM_LABELS[parameter],
    corr_label = paste0("r = ", sprintf("%.2f", r))
  )

combined_plot <- ggplot(recovery_params, aes(x = value_true, y = value_pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
  facet_wrap(
    ~ parameter,
    scales = "free",
    ncol = 5,
    labeller = as_labeller(PARAM_LABELS, default = label_parsed)
  ) +
  geom_text(
    data = recovery_params %>% group_by(parameter) %>% summarize(
      x = min(value_true),
      y = max(value_pred),
      corr_label = unique(corr_label)
    ),
    mapping = aes(x = x, y = y, label = corr_label),
    hjust = 0, vjust = 1,
    size = 6,
    fontface = "italic",
    family = "Palatino",
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(x = "True", y = "Estimated") +
  scale_x_continuous(n.breaks = 4) +
  ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(linewidth = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_2),
    strip.text.y = element_text(size = FONT_SIZE_2, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    plot.title = element_text(size = FONT_SIZE_1,
                              hjust = 0.5,
                              face = 'bold'),
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank(),
    legend.spacing.y = unit(0.25, 'cm'),
    panel.spacing = unit(1., "lines")
  )

ggsave(
  "../plots/06_parameter_recovery_new_individual_plot.jpeg",
  plot = combined_plot,
  width = 14, height = 3.5, dpi = 300,
  device = "jpeg"
)
