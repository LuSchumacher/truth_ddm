library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
library(cmdstanr)
library(rempsyc)
library(flextable)
library(LaplacesDemon)
library(truncnorm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fit_session_1 <- readRDS("../../fits/fit_session_1_ndt_var.rds")
fit_session_2 <- readRDS("../../fits/fit_session_2_ndt_var.rds")
fit_exp_2 <- readRDS("../../fits/fit_exp_2_ndt_var.rds")

PARAM_NAMES <- c(
  "transf_mu_v[1]", "transf_mu_v[2]", "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]", "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]", "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]", "transf_mu_bias[3]", "transf_mu_bias[4]",
  "transf_mu_ndt_s"
)

NUM_POST_SAMPLES <- 12000

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
COLOR_PALETTE <- c('#27374D', '#B70404', '#6F1E29', '#B78F94')

custom_labeller <- labeller(
  parameter = label_parsed,
  effect = label_value
)

softplus <- function(x) {
  log1p(exp(x))
}

density_at_zero <- function(x) {
  dens_obj <- density(x, from = min(0, min(x)) - 1, to = max(0, max(x)) + 1)
  approx(dens_obj$x, dens_obj$y, xout = 0)$y
}

NUM_PRIOR_SAMPLES <- 1e7

mu_v1_prior <- rnorm(NUM_PRIOR_SAMPLES, 2, 2)
mu_v2_prior <- rnorm(NUM_PRIOR_SAMPLES, 2, 2)
mu_v3_prior <- rnorm(NUM_PRIOR_SAMPLES, 2, 2)
mu_v4_prior <- rnorm(NUM_PRIOR_SAMPLES, 2, 2)

contrast_v_prior <- ((mu_v1_prior + mu_v2_prior) / 2) - ((mu_v3_prior + mu_v4_prior) / 2)
p_v_prior <- density_at_zero(contrast_v_prior)

mu_a1_prior <- rnorm(NUM_PRIOR_SAMPLES, 5, 3)
mu_a2_prior <- rnorm(NUM_PRIOR_SAMPLES, 5, 3)
mu_a3_prior <- rnorm(NUM_PRIOR_SAMPLES, 5, 3)
mu_a4_prior <- rnorm(NUM_PRIOR_SAMPLES, 5, 3)

a1_prior <- softplus(mu_a1_prior)
a2_prior <- softplus(mu_a2_prior)
a3_prior <- softplus(mu_a3_prior)
a4_prior <- softplus(mu_a4_prior)

contrast_a_prior <- ((a1_prior + a2_prior) / 2) - ((a3_prior + a4_prior) / 2)
p_a_prior <- density_at_zero(contrast_a_prior)

mu_bias1_prior <- rnorm(NUM_PRIOR_SAMPLES, 0, 0.5)
mu_bias2_prior <- rnorm(NUM_PRIOR_SAMPLES, 0, 0.5)
mu_bias3_prior <- rnorm(NUM_PRIOR_SAMPLES, 0, 0.5)
mu_bias4_prior <- rnorm(NUM_PRIOR_SAMPLES, 0, 0.5)

bias1_prior <- invlogit(mu_bias1_prior)
bias2_prior <- invlogit(mu_bias2_prior)
bias3_prior <- invlogit(mu_bias3_prior)
bias4_prior <- invlogit(mu_bias4_prior)

contrast_bias_prior <- ((bias1_prior + bias2_prior) / 2) - ((bias3_prior + bias4_prior) / 2)
p_bias_prior <- density_at_zero(contrast_bias_prior)

mu_ndt1_prior <- rtruncnorm(NUM_PRIOR_SAMPLES, a=0, mean=1, sd=2)
mu_ndt2_prior <- rtruncnorm(NUM_PRIOR_SAMPLES, a=0, mean=1, sd=2)
mu_ndt3_prior <- rtruncnorm(NUM_PRIOR_SAMPLES, a=0, mean=1, sd=2)
mu_ndt4_prior <- rtruncnorm(NUM_PRIOR_SAMPLES, a=0, mean=1, sd=2)

contrast_ndt_prior <- ((mu_ndt1_prior + mu_ndt2_prior) / 2) - ((mu_ndt3_prior + mu_ndt4_prior) / 2)
p_ndt_prior <- density_at_zero(contrast_ndt_prior)

#------------------------------------------------------------------------------#
# Session 1
#------------------------------------------------------------------------------#
s1_post_samples <- fit_session_1$draws(
  variables = PARAM_NAMES,
  inc_warmup = FALSE,
  format = "draws_matrix"
)
s1_post_samples <- as_tibble(s1_post_samples) %>% 
  mutate(draw = 1:NUM_POST_SAMPLES) %>% 
  pivot_longer(
    cols = starts_with("transf"),
    names_to = c("parameter", "condition"),
    names_pattern = "(mu_\\w+)\\[(\\d+)\\]"
  )
s1_post_samples <- s1_post_samples %>% 
  mutate(
    factual_truth = case_when(
      condition == 1 | condition == 2 ~ "True",
      condition == 3 | condition == 4 ~ "False"
    ),
    stim_type = case_when(
      condition == 1 | condition == 3 ~ "New statement",
      condition == 2 | condition == 4 ~ "Repeated statement"
    ),
    value = as.numeric(value)
  )

s1_post_samples$parameter[is.na(s1_post_samples$parameter)] <- "mu_ndt_sd"

s1_post_samples <- s1_post_samples %>%
  mutate(parameter = factor(
    parameter,
    levels = c("mu_v", "mu_a", "mu_ndt", "mu_bias", "mu_ndt_sd"),
    labels = c(
      expression(mu[v]), expression(mu[a]),
      expression(mu[ndt]), expression(mu[bias]),
      expression(mu[ndt[sd]])
    )
    )
  )

post_summaries_session_1 <- s1_post_samples %>%
  group_by(parameter, condition) %>%
  summarise(
    median = median(value),
    q5 = quantile(value, probs = 0.025),
    q95 = quantile(value, probs = 0.975),
    .groups = "drop"
  )  %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(post_summaries_session_1, "../../data/post_summaries_session_1_ndt_var.csv")

estimates_plot_s1 <- s1_post_samples %>% 
  ggplot(aes(x = factual_truth, y = value, colour = stim_type, fill = stim_type)) +
  geom_violin() +
  facet_wrap(~parameter, scales = "free", labeller = custom_labeller) +
  theme_classic() +
  scale_color_manual(guide = "none", values = COLOR_PALETTE) +
  scale_fill_manual(name = "Repetition status", values = COLOR_PALETTE) + 
  labs(
    x = 'Factual truth',
    y = 'Value',
    title = 'Group−level Posterior Distributions: Experiment 1, Session 1'
  ) +
  ggthemes::theme_tufte() + 
  theme(axis.line = element_line(linewidth = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text = element_text(size = FONT_SIZE_1),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid.major = element_line(color = alpha("gray70", 0.3)),
        panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b =05, l = 0))
  )

ggsave(
  '../../plots/01_param_estimates_s1_ndt_var.jpeg',
  estimates_plot_s1,
  device = 'jpeg', dpi = 300,
  width = 12, height = 6
)

s1_post_samples <- s1_post_samples %>% 
  filter(parameter != "mu[ndt[sd]]")

## Contrasts
contrast_df_s1 <- tibble()
for (name in unique(s1_post_samples$parameter)){
  # main effect repetition
  effect_repetition <- ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 2] +
                           s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 4]) / 2) -
    ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 1] +
        s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 3]) / 2)
  # main effect factual truth
  effect_truth <- ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 1] +
                      s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 2]) / 2) -
    ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 3] +
        s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 4]) / 2)
  # interaction
  effect_interaction <- ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 1] +
                            s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 4]) / 2) -
    ((s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 2] +
        s1_post_samples$value[s1_post_samples$parameter == name & s1_post_samples$condition == 3]) / 2)
  tmp_df <- tibble(
    "parameter" = rep(name, NUM_POST_SAMPLES*3),
    "effect" = rep(c("Repetition status", "Factual truth", "Interaction"), each = NUM_POST_SAMPLES),
    "value" = c(effect_repetition, effect_truth, effect_interaction)
  )
  contrast_df_s1 <- rbind(contrast_df_s1, tmp_df)
}

contrast_session_01 <- contrast_df_s1 %>%
  group_by(parameter, effect) %>%
  summarise(
    median = median(value),
    ci_low = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    p_post = density_at_zero(value)
  ) %>% 
  ungroup()

# Bayes factors
contrast_session_01$bf <- NA
# threshold
contrast_session_01$bf[1] <- p_a_prior / contrast_session_01$p_post[1]
contrast_session_01$bf[2] <- p_a_prior / contrast_session_01$p_post[2]
contrast_session_01$bf[3] <- p_a_prior / contrast_session_01$p_post[3]
# bias
contrast_session_01$bf[4] <- p_bias_prior / contrast_session_01$p_post[4]
contrast_session_01$bf[5] <- p_bias_prior / contrast_session_01$p_post[5]
contrast_session_01$bf[6] <- p_bias_prior / contrast_session_01$p_post[6]
# ndt
contrast_session_01$bf[7] <- p_ndt_prior / contrast_session_01$p_post[7]
contrast_session_01$bf[8] <- p_ndt_prior / contrast_session_01$p_post[8]
contrast_session_01$bf[9] <- p_ndt_prior / contrast_session_01$p_post[9]
# drift
contrast_session_01$bf[10] <- p_v_prior / contrast_session_01$p_post[10]
contrast_session_01$bf[11] <- p_v_prior / contrast_session_01$p_post[11]
contrast_session_01$bf[12] <- p_v_prior / contrast_session_01$p_post[12]

contrast_session_01 <- contrast_session_01 %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(contrast_session_01, "../../data/contrast_session_01_ndt_var.csv")

#------------------------------------------------------------------------------#
# Session 2
#------------------------------------------------------------------------------#
s2_post_samples <- fit_session_2$draws(
  variables = PARAM_NAMES,
  inc_warmup = FALSE,
  format = "draws_matrix"
)

s2_post_samples <- as_tibble(s2_post_samples) %>% 
  mutate(draw = 1:NUM_POST_SAMPLES) %>% 
  pivot_longer(
    cols = starts_with("transf"),
    names_to = c("parameter", "condition"),
    names_pattern = "(mu_\\w+)\\[(\\d+)\\]"
  )

s2_post_samples <- s2_post_samples %>% 
  mutate(
    factual_truth = case_when(
      condition == 1 | condition == 2 ~ "True",
      condition == 3 | condition == 4 ~ "False"
    ),
    stim_type = case_when(
      condition == 1 | condition == 3 ~ "New statement",
      condition == 2 | condition == 4 ~ "Repeated statement"
    )
  )

s2_post_samples$parameter[is.na(s2_post_samples$parameter)] <- "mu_ndt_sd"

s2_post_samples <- s2_post_samples %>%
  mutate(parameter = factor(
    parameter,
    levels = c("mu_v", "mu_a", "mu_ndt", "mu_bias", "mu_ndt_sd"),
    labels = c(
      expression(mu[v]), expression(mu[a]),
      expression(mu[ndt]), expression(mu[bias]),
      expression(mu[ndt[sd]])
    )
  )
  )

post_summaries_session_2 <- s2_post_samples %>%
  group_by(parameter, condition) %>%
  summarise(
    median = median(value),
    q5 = quantile(value, probs = 0.025),
    q95 = quantile(value, probs = 0.975),
    .groups = "drop"
  )  %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(post_summaries_session_2, "../../data/post_summaries_session_2_ndt_var.csv")

estimates_plot_s2 <- s2_post_samples %>% 
  ggplot(aes(x = factual_truth, y = value, colour = stim_type, fill = stim_type)) +
  geom_violin() +
  facet_wrap(~parameter, scales = "free", labeller = custom_labeller) +
  theme_classic() +
  scale_color_manual(guide = "none", values = COLOR_PALETTE) +
  scale_fill_manual(name = "Repetition status", values = COLOR_PALETTE) + 
  labs(
    x = 'Factual truth',
    y = 'Value',
    title = 'Group−level Posterior Distributions: Experiment 1, Session 2'
  ) +
  ggthemes::theme_tufte() + 
  theme(axis.line = element_line(linewidth = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text = element_text(size = FONT_SIZE_1),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid.major = element_line(color = alpha("gray70", 0.3)),
        panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b =05, l = 0))
  )

ggsave(
  '../../plots/01_param_estimates_s2_ndt_var.jpeg',
  estimates_plot_s2,
  device = 'jpeg', dpi = 300,
  width = 12, height = 6
)

s2_post_samples <- s2_post_samples %>% 
  filter(parameter != "mu[ndt[sd]]")

## Contrasts
contrast_df_s2 <- tibble()
for (name in unique(s2_post_samples$parameter)){
  # main effect repetition
  effect_repetition <- ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 2] +
                           s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 4]) / 2) -
    ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 1] +
        s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 3]) / 2)
  # main effect factual truth
  effect_truth <- ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 1] +
                      s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 2]) / 2) -
    ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 3] +
        s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 4]) / 2)
  # interaction
  effect_interaction <- ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 1] +
                            s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 4]) / 2) -
    ((s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 2] +
        s2_post_samples$value[s2_post_samples$parameter == name & s2_post_samples$condition == 3]) / 2)
  tmp_df <- tibble(
    "parameter" = rep(name, NUM_POST_SAMPLES*3),
    "effect" = rep(c("Repetition status", "Factual truth", "Interaction"), each = NUM_POST_SAMPLES),
    "value" = c(effect_repetition, effect_truth, effect_interaction)
  )
  contrast_df_s2 <- rbind(contrast_df_s2, tmp_df)
}

contrast_session_2 <- contrast_df_s2 %>%
  group_by(parameter, effect) %>%
  summarise(
    median = median(value),
    ci_low = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    p_post = density_at_zero(value)
  ) %>% 
  ungroup()

# Bayes factors
contrast_session_2$bf <- NA
# threshold
contrast_session_2$bf[1] <- p_a_prior / contrast_session_2$p_post[1]
contrast_session_2$bf[2] <- p_a_prior / contrast_session_2$p_post[2]
contrast_session_2$bf[3] <- p_a_prior / contrast_session_2$p_post[3]
# bias
contrast_session_2$bf[4] <- p_bias_prior / contrast_session_2$p_post[4]
contrast_session_2$bf[5] <- p_bias_prior / contrast_session_2$p_post[5]
contrast_session_2$bf[6] <- p_bias_prior / contrast_session_2$p_post[6]
# ndt
contrast_session_2$bf[7] <- p_ndt_prior / contrast_session_2$p_post[7]
contrast_session_2$bf[8] <- p_ndt_prior / contrast_session_2$p_post[8]
contrast_session_2$bf[9] <- p_ndt_prior / contrast_session_2$p_post[9]
# drift
contrast_session_2$bf[10] <- p_v_prior / contrast_session_2$p_post[10]
contrast_session_2$bf[11] <- p_v_prior / contrast_session_2$p_post[11]
contrast_session_2$bf[12] <- p_v_prior / contrast_session_2$p_post[12]

contrast_session_2 <- contrast_session_2 %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(contrast_session_2, "../../data/contrast_session_02_ndt_var.csv")

contrast_df_s1$session <- 1
contrast_df_s2$session <- 2
contrast_df <- rbind(contrast_df_s1, contrast_df_s2)

contrast_df <- contrast_df %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("mu[v]", "mu[a]", "mu[ndt]", "mu[bias]"),
      labels = c(expression(mu[v]), expression(mu[a]), expression(mu[ndt]), expression(mu[bias]))
    ),
    session = factor(session)
  )

contrast_plot_exp_1 <- contrast_df %>% 
  ggplot(aes(x = value, fill = session)) +
  geom_density(alpha = 0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(parameter ~ effect, labeller = custom_labeller) +
  scale_x_continuous(limits = c(-0.7, 0.5)) +
  scale_fill_manual(values = COLOR_PALETTE) +
  labs(
    x = "Value",
    y = "Density",
    title = 'Group-level Parameter Contrasts: Experiment 1'
  ) +
  ggthemes::theme_tufte() + 
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3,
                                   angle = 45,
                                   vjust = 0.5),
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
        panel.spacing = unit(1., "lines")) + 
  guides(fill = guide_legend(title = "Session"))

ggsave(
  '../../plots/02_contrast_plot_exp_1_ndt_var.jpeg',
  contrast_plot_exp_1,
  device = 'jpeg', dpi = 300,
  width = 12, height = 7
)

#------------------------------------------------------------------------------#
# Experiment 2
#------------------------------------------------------------------------------#
exp2_post_samples <- fit_exp_2$draws(
  variables = PARAM_NAMES,
  inc_warmup = FALSE,
  format = "draws_matrix"
)
exp2_post_samples <- as_data_frame(exp2_post_samples) %>% 
  mutate(draw = 1:NUM_POST_SAMPLES) %>% 
  pivot_longer(
    cols = starts_with("transf"),
    names_to = c("parameter", "condition"),
    names_pattern = "(mu_\\w+)\\[(\\d+)\\]"
  )
exp2_post_samples <- exp2_post_samples %>% 
  mutate(
    factual_truth = case_when(
      condition == 1 | condition == 2 ~ "True",
      condition == 3 | condition == 4 ~ "False"
    ),
    stim_type = case_when(
      condition == 1 | condition == 3 ~ "New statement",
      condition == 2 | condition == 4 ~ "Repeated statement"
    ),
    value = as.numeric(value)
  )

exp2_post_samples$parameter[is.na(exp2_post_samples$parameter)] <- "mu_ndt_sd"

exp2_post_samples <- exp2_post_samples %>%
  mutate(parameter = factor(
    parameter,
    levels = c("mu_v", "mu_a", "mu_ndt", "mu_bias", "mu_ndt_sd"),
    labels = c(
      expression(mu[v]), expression(mu[a]),
      expression(mu[ndt]), expression(mu[bias]),
      expression(mu[ndt[sd]])
    )
  )
  )

post_summaries_exp_2 <- exp2_post_samples %>%
  group_by(parameter, condition) %>%
  summarise(
    median = median(value),
    q5 = quantile(value, probs = 0.025),
    q95 = quantile(value, probs = 0.975),
    .groups = "drop"
  )  %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(post_summaries_exp_2, "../../data/post_summaries_exp_2_ndt_var.csv")

estimates_plot_exp2 <- exp2_post_samples %>% 
  ggplot(aes(x = factual_truth, y = value, colour = stim_type, fill = stim_type)) +
  geom_violin() +
  facet_wrap(~parameter, scales = "free", labeller = custom_labeller) +
  theme_classic() +
  scale_color_manual(guide = "none", values = COLOR_PALETTE) +
  scale_fill_manual(name = "Repetition status", values = COLOR_PALETTE) + 
  labs(
    x = 'Factual truth',
    y = 'Value',
    title = 'Group−level Posterior Distributions: Experiment 2'
  ) +
  ggthemes::theme_tufte() + 
  theme(axis.line = element_line(linewidth = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text = element_text(size = FONT_SIZE_1),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid.major = element_line(color = alpha("gray70", 0.3)),
        panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b =05, l = 0))
  )

ggsave(
  '../../plots/01_param_estimates_exp2_ndt_var.jpeg',
  estimates_plot_exp2,
  device = 'jpeg', dpi = 300,
  width = 12, height = 6
)

exp2_post_samples <- exp2_post_samples %>% 
  filter(parameter != "mu[ndt[sd]]")

## Contrasts
contrast_df_exp2 <- tibble()
for (name in unique(exp2_post_samples$parameter)){
  # main effect repetition
  effect_repetition <- ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 2] +
                           exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 4]) / 2) -
    ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 1] +
        exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 3]) / 2)
  # main effect factual truth
  effect_truth <- ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 1] +
                      exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 2]) / 2) -
    ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 3] +
        exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 4]) / 2)
  # interaction
  effect_interaction <- ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 1] +
                            exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 4]) / 2) -
    ((exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 2] +
        exp2_post_samples$value[exp2_post_samples$parameter == name & exp2_post_samples$condition == 3]) / 2)
  tmp_df <- tibble(
    "parameter" = rep(name, NUM_POST_SAMPLES*3),
    "effect" = rep(c("Repetition status", "Factual truth", "Interaction"), each = NUM_POST_SAMPLES),
    "value" = c(effect_repetition, effect_truth, effect_interaction)
  )
  contrast_df_exp2 <- rbind(contrast_df_exp2, tmp_df)
}

contrast_exp2 <- contrast_df_exp2 %>%
  group_by(parameter, effect) %>%
  summarise(
    median = median(value),
    ci_low = quantile(value, 0.025),
    ci_high = quantile(value, 0.975),
    p_post = density_at_zero(value)
  ) %>% 
  ungroup()

# Bayes factors
contrast_exp2$bf <- NA
# threshold
contrast_exp2$bf[1] <- p_a_prior / contrast_exp2$p_post[1]
contrast_exp2$bf[2] <- p_a_prior / contrast_exp2$p_post[2]
contrast_exp2$bf[3] <- p_a_prior / contrast_exp2$p_post[3]
# bias
contrast_exp2$bf[4] <- p_bias_prior / contrast_exp2$p_post[4]
contrast_exp2$bf[5] <- p_bias_prior / contrast_exp2$p_post[5]
contrast_exp2$bf[6] <- p_bias_prior / contrast_exp2$p_post[6]
# ndt
contrast_exp2$bf[7] <- p_ndt_prior / contrast_exp2$p_post[7]
contrast_exp2$bf[8] <- p_ndt_prior / contrast_exp2$p_post[8]
contrast_exp2$bf[9] <- p_ndt_prior / contrast_exp2$p_post[9]
# drift
contrast_exp2$bf[10] <- p_v_prior / contrast_exp2$p_post[10]
contrast_exp2$bf[11] <- p_v_prior / contrast_exp2$p_post[11]
contrast_exp2$bf[12] <- p_v_prior / contrast_exp2$p_post[12]

contrast_exp2 <- contrast_exp2 %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(contrast_exp2, "../../data/contrast_exp2_ndt_var.csv")

contrast_df_exp2 <- contrast_df_exp2 %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("mu[v]", "mu[a]", "mu[ndt]", "mu[bias]"),
      labels = c(expression(mu[v]), expression(mu[a]), expression(mu[ndt]), expression(mu[bias]))
    )
  )

contrast_plot_exp_2 <- contrast_df_exp2 %>% 
  ggplot(aes(x = value)) +
  geom_density(alpha = 0.75, fill = COLOR_PALETTE[1]) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(parameter ~ effect, labeller = custom_labeller) +
  scale_x_continuous(limits = c(-0.7, 0.5)) + 
  labs(
    x = "Value",
    y = "Density",
    title = 'Group-level Parameter Contrasts: Experiment 2'
  ) +
  ggthemes::theme_tufte() + 
  theme(
    axis.line = element_line(size = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3,
                               angle = 45,
                               vjust = 0.5),
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
  '../../plots/02_contrast_plot_exp_2_ndt_var.jpeg',
  contrast_plot_exp_2,
  device = 'jpeg', dpi = 300,
  width = 12, height = 7
)
