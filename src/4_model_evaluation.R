library(tidyverse)
library(magrittr)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(rempsyc)
library(flextable)
library(purrr)

options(scipen = 999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

summarize_param <- function(draws, param) {
  x <- draws[, param]
  med <- median(x)
  ci <- quantile(x, c(0.025, 0.975))
  sprintf("%.2f [%.2f, %.2f]", med, ci[1], ci[2])
}

calc_BF <- function(posterior_samples, prior_sd = 0.5) {
  post_density <- density(posterior_samples, bw = "nrd0")
  posterior_at_zero <- approx(post_density$x, post_density$y, xout = 0)$y
  if(is.na(posterior_at_zero) || posterior_at_zero < 1e-10) {
    posterior_at_zero <- 1e-10
  }
  prior_at_zero <- dnorm(0, mean = 0, sd = prior_sd)
  BF10 <- prior_at_zero / posterior_at_zero
  return(BF10)
}

compute_BFs <- function(fit, param_names, prior_sd = 0.5) {
  sapply(param_names, function(param) {
    samples <- as.numeric(fit$draws(param))
    calc_BF(samples, prior_sd)
  })
}

extract_posteriors <- function(fit, param_names, session_label) {
  df_list <- lapply(param_names, function(p) {
    data.frame(
      beta = p,
      value = as.numeric(fit$draws(p)),
      session = session_label
    )
  })
  bind_rows(df_list)
}

summarize_param <- function(fit, param_name, cond_indices) {
  draws <- as.matrix(fit)[, param_name]
  sapply(cond_indices, function(idx) {
    median_val <- round(median(draws[, idx]), 2)
    ci <- round(quantile(draws[, idx], probs = c(0.025, 0.975)), 2)
    sprintf("%0.2f [%0.2f, %0.2f]", median_val, ci[1], ci[2])
  })
}

fit_session_1 <- readRDS("../fits/fit_session_1.rds")
fit_session_2 <- readRDS("../fits/fit_session_2.rds")
fit_exp_2 <- readRDS("../fits/fit_exp_2.rds")

draws_session_1 <- as_draws_matrix(fit_session_1)
draws_session_2 <- as_draws_matrix(fit_session_2)
draws_exp_2 <- as_draws_matrix(fit_exp_2)

PARAM_NAMES <- c(
  "v_betas[1]", "v_betas[2]", "v_betas[3]",
  "a_betas[1]", "a_betas[2]", "a_betas[3]",
  "bias_betas[1]", "bias_betas[2]", "bias_betas[3]",
  "ndt_betas[1]", "ndt_betas[2]", "ndt_betas[3]",
  "ndt_var_betas[1]", "ndt_var_betas[2]", "ndt_var_betas[3]"
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
COLOR_PALETTE <- c('#27374D', '#B70404')

################################################################################
# POSTERIOR SUMMARIES
################################################################################
summary_session_1 <- sapply(PARAM_NAMES, summarize_param, draws = draws_session_1)
summary_session_2 <- sapply(PARAM_NAMES, summarize_param, draws = draws_session_2)
summary_exp_2     <- sapply(PARAM_NAMES, summarize_param, draws = draws_exp_2)

posterior_table <- tibble(
  parameter = PARAM_NAMES,
  Session_1 = summary_session_1,
  Session_2 = summary_session_2,
  Exp_2 = summary_exp_2
) %>% 
  mutate(
    effect = str_extract(parameter, "(?<=\\[)\\d+(?=\\])"),
    parameter = str_remove(parameter, "_[^_]+$"),
    effect = case_when(
      effect == "1" ~ "Repetition",
      effect == "2" ~ "Factual truth",
      effect == "3" ~ "Interaction",
      TRUE ~ NA_character_
    )
  ) %>% 
  relocate(effect, .after = parameter)

write_csv(posterior_table, "../tables/effects_post_summaries_table.csv")

param_names <- c(
  "transf_mu_v[1]", "transf_mu_v[2]", "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]", "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]", "transf_mu_bias[3]", "transf_mu_bias[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]", "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_ndt_var[1]", "transf_mu_ndt_var[2]", "transf_mu_ndt_var[3]", "transf_mu_ndt_var[4]"
)

summary_session_1 <- sapply(param_names, summarize_param, draws = draws_session_1)
summary_session_2 <- sapply(param_names, summarize_param, draws = draws_session_2)
summary_exp_2     <- sapply(param_names, summarize_param, draws = draws_exp_2)

posterior_table <- tibble(
  Parameter = param_names,
  Session_1 = summary_session_1,
  Session_2 = summary_session_2,
  Exp_2 = summary_exp_2
) %>% 
  mutate(
    Condition_index = str_extract(Parameter, "(?<=\\[)\\d+(?=\\])"),
    Condition_index = as.integer(Condition_index),
    Parameter = str_remove(Parameter, "transf_mu_"),
    Parameter = str_remove(Parameter, "\\[\\d+\\]")
  )

condition_map <- tibble(
  Condition_index = 1:4,
  Factual_Truth = c("True", "True", "False", "False"),
  Repetition    = c("New", "Repeated", "New", "Repeated")
)

posterior_table <- posterior_table %>%
  left_join(condition_map, by = "Condition_index") %>%
  select(Parameter, Factual_Truth, Repetition, Session_1, Session_2, Exp_2)

write_csv(posterior_table, "../tables/post_summaries_table.csv")

################################################################################
# BAYES FACTORS
################################################################################
BF_session_1 <- compute_BFs(fit_session_1, PARAM_NAMES)
BF_session_2 <- compute_BFs(fit_session_2, PARAM_NAMES)
BF_exp_2     <- compute_BFs(fit_exp_2, PARAM_NAMES)

BF_table <- tibble(
  Parameter = PARAM_NAMES,
  BF_session_1 = BF_session_1,
  BF_session_2 = BF_session_2,
  BF_exp_2 = BF_exp_2
) %>% 
  mutate(
    param = sub("_betas\\[.*\\]", "", Parameter),
    beta_index = gsub(".*\\[|\\]", "", Parameter),
    effect = case_when(
      beta_index == "1" ~ "Repetition",
      beta_index == "2" ~ "Truth",
      beta_index == "3" ~ "Interaction"
    ),
    BF_session_1 = round(BF_session_1, 2),
    BF_session_2 = round(BF_session_2, 2),
    BF_exp_2     = round(BF_exp_2, 2)
  ) %>%
  select(param, effect, BF_session_1, BF_session_2, BF_exp_2) %>%
  arrange(factor(param, levels = c("v","a","bias","ndt","ndt_var")),
          factor(effect, levels = c("Repetition","Truth","Interaction")))

write_csv(BF_table, "../tables/ddm_bayes_factors.csv")

################################################################################
# PLOT EFFECTS: EXPERIMENT 1
################################################################################
post_s1 <- extract_posteriors(fit_session_1, PARAM_NAMES, "Session 1")
post_s2 <- extract_posteriors(fit_session_2, PARAM_NAMES, "Session 2")
post <- bind_rows(post_s1, post_s2)

post <- post %>%
  mutate(
    param = case_when(
      grepl("^v_betas", beta) ~ "v",
      grepl("^a_betas", beta) ~ "a",
      grepl("^bias_betas", beta) ~ "bias",
      grepl("^ndt_betas", beta) ~ "ndt",
      grepl("^ndt_var_betas", beta) ~ "ndt_var"
    ),
    effect = case_when(
      grepl("\\[1\\]", beta) ~ "Repetition status",
      grepl("\\[2\\]", beta) ~ "Factual truth",
      grepl("\\[3\\]", beta) ~ "Interaction"
    ),
    param = factor(param, levels = c("v", "a", "bias", "ndt", "ndt_var")),
    effect = factor(effect, levels = c("Repetition status", "Factual truth", "Interaction")),
    value_diff = 2 * value
  )

effects_plot_exp1 <- ggplot(post, aes(x = value_diff, fill = session, color = session)) +
  geom_density(alpha = 0.65, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(
    param ~ effect,
    scales = "free_x",
    labeller = labeller(param = as_labeller(PARAM_LABELS, default = label_parsed))
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_color_manual(values = COLOR_PALETTE, guide = "none") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(x = "Effect", y = "Density",
       fill = "") +
  ggthemes::theme_tufte(base_size = FONT_SIZE_2) + 
  theme(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.line = element_line(linewidth = 0.5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3, vjust = 0.5),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_2),
    strip.text.y = element_text(size = FONT_SIZE_2, hjust = 0, angle = 0),
    panel.grid.major = element_line(color = scales::alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = scales::alpha("gray70", 0.15)),
    panel.background = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(0.2, "cm")
  )

ggsave(
  '../plots/02_effects_exp_1.jpeg',
  effects_plot_exp1,
  device = 'jpeg', dpi = 300,
  width = 12, height = 9
)

################################################################################
# PLOT EFFECTS: EXPERIMENT 2
################################################################################
post_exp2 <- extract_posteriors(fit_exp_2, PARAM_NAMES, "Exp 2")

post_exp2 <- post_exp2 %>%
  mutate(
    param = case_when(
      grepl("^v_betas", beta) ~ "v",
      grepl("^a_betas", beta) ~ "a",
      grepl("^bias_betas", beta) ~ "bias",
      grepl("^ndt_betas", beta) ~ "ndt",
      grepl("^ndt_var_betas", beta) ~ "ndt_var"
    ),
    effect = case_when(
      grepl("\\[1\\]", beta) ~ "Repetition status",
      grepl("\\[2\\]", beta) ~ "Factual truth",
      grepl("\\[3\\]", beta) ~ "Interaction"
    ),
    param = factor(param, levels = c("v", "a", "bias", "ndt", "ndt_var")),
    effect = factor(effect, levels = c("Repetition status", "Factual truth", "Interaction")),
    value_diff = 2 * value
  )

effects_plot_exp2 <- ggplot(post_exp2, aes(x = value_diff)) +
  geom_density(fill = COLOR_PALETTE[1], color = COLOR_PALETTE[1], alpha = 0.65, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(
    param ~ effect,
    scales = "free_x",
    labeller = labeller(param = as_labeller(PARAM_LABELS, default = label_parsed))
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(x = "Effect", y = "Density", fill = "") +
  ggthemes::theme_tufte(base_size = FONT_SIZE_2) +
  theme(
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.line = element_line(linewidth = 0.5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3, vjust = 0.5),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_2),
    strip.text.y = element_text(size = FONT_SIZE_2, hjust = 0, angle = 0),
    panel.grid.major = element_line(color = scales::alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = scales::alpha("gray70", 0.15)),
    panel.background = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "none"
  )

ggsave(
  '../plots/02_effects_exp_2.jpeg',
  effects_plot_exp2,
  device = 'jpeg', dpi = 300,
  width = 12, height = 9
)
