library(brms)
library(tidyverse)
library(latex2exp)
library(lme4)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

PARAM_NAMES <- c(
  TeX('Intercept'), TeX(r'($\beta_{Session}$)'),
  TeX(r'($\beta_{Repetition \, status}$)'), TeX(r'($\beta_{Factual \, truth}$)'),
  TeX(r'($\beta_{Session:Repetition \, status}$)'),
  TeX(r'($\beta_{Repetition \, status:Factual \, truth}$)'), TeX(r'($\beta_{Session:Factual \, truth}$)'),
  TeX(r'($\beta_{Session:Repetition \, status:Factual \, truth}$)')
)

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
COLOR_PALETTE <- c('#27374D', '#B70404')

#------------------------------------------------------------------------------#
# Retention interval: Experiment 1
#------------------------------------------------------------------------------#
df_session_1 <- read_csv('../data/data_session_1.csv')
df_session_2 <- read_csv('../data/data_session_2.csv')

df_exp1 <- rbind(df_session_1, df_session_2) %>% 
  mutate(
    session = as_factor(session),
    stim_type = as_factor(stim_type),
    factual_truth = as_factor(factual_truth)
  )

contrasts(df_exp1$session) <- rbind(-1.0, 1.0)
contrasts(df_exp1$factual_truth) <- rbind(-1.0, 1.0)
contrasts(df_exp1$stim_type) <- rbind(-1.0, 1.0)

formula <- resp ~ stim_type * session * factual_truth + (1 + stim_type * session * factual_truth | id)

priors <- prior(normal(0, 1.0), class = b)
model_exp1 <- brm(
  formula = formula,
  data = df_exp1,
  family = bernoulli("logit"),
  prior = priors,
  iter = 4000,
  cores = 4,
  chains = 4,
  file = "../fits/hlm_fit_exp1",
  sample_prior = "yes"
)

model_exp1

# hypothesis test
hyp_session <- "session1 = 0"
test_session <- hypothesis(model_exp1, hyp_session)$hypothesis
bf_session <- 1 / test_session["Evid.Ratio"]

hyp_rep <- "stim_type1 = 0"
test_rep <- hypothesis(model_exp1, hyp_rep)$hypothesis
bf_rep <- 1 / test_rep["Evid.Ratio"]

hyp_factual_truth <- "factual_truth1 = 0"
test_factual_truth <- hypothesis(model_exp1, hyp_factual_truth)$hypothesis
bf_factual_truth <- 1 / test_factual_truth["Evid.Ratio"]

hyp_type_session <- "stim_type1:session1 = 0"
test_type_session <- hypothesis(model_exp1, hyp_type_session)$hypothesis
bf_type_session <- 1 / test_type_session["Evid.Ratio"]

hyp_type_truth <- "stim_type1:factual_truth1 = 0"
test_type_truth <- hypothesis(model_exp1, hyp_type_truth)$hypothesis
bf_type_truth <- 1 / test_type_truth["Evid.Ratio"]

hyp_session_truth <- "session1:factual_truth1 = 0"
test_session_truth <- hypothesis(model_exp1, hyp_session_truth)$hypothesis
bf_session_truth <- 1 / test_session_truth["Evid.Ratio"]

hyp_three_way <- "stim_type1:session1:factual_truth1 = 0"
test_three_way <- hypothesis(model_exp1, hyp_three_way)$hypothesis
bf_three_way <- 1 / test_three_way["Evid.Ratio"]

# ----- Extract posterior samples -----
fixed_effects <- brms::fixef(model_exp1, summary = FALSE) %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "post_sample") %>%
  mutate(parameter = factor(parameter))

# ----- Compute posterior density at 0 -----
fixed_effects$y <- NA
params <- unique(fixed_effects$parameter)
for (param in params) {
  post_samples <- fixed_effects$post_sample[fixed_effects$parameter == param]
  d <- density(post_samples)
  # rule = 2 allows linear extrapolation if 0 is outside range
  fixed_effects$y[fixed_effects$parameter == param] <- approx(d$x, d$y, xout = 0, rule = 2)$y
}

# ----- Compute Bayes factors -----
# (assuming you've already computed bf_session, bf_rep, etc.)
fixed_effects$bf <- NA
fixed_effects$bf[fixed_effects$parameter == "session1"] <- paste("BF =", round(bf_session, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "factual_truth1"] <- paste("BF =", round(bf_factual_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:session1"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "stim_type1:factual_truth1"] <- paste("BF =", round(bf_type_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "session1:factual_truth1"] <- paste("BF =", round(bf_session_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:session1:factual_truth1"] <- paste("BF =", round(bf_three_way, 2))

# ----- Set text positions for labels -----
fixed_effects$text_loc_x <- NA
fixed_effects$text_loc_x[fixed_effects$parameter == "session1"] <- 0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1"] <- 0.0
fixed_effects$text_loc_x[fixed_effects$parameter == "factual_truth1"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:session1"] <- 0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:factual_truth1"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "session1:factual_truth1"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:session1:factual_truth1"] <- -0.3

# ----- Factor levels & labels -----
fixed_effects$parameter <- factor(
  fixed_effects$parameter,
  levels = c(
    'Intercept', 'session1', 'stim_type1', 'factual_truth1',
    'stim_type1:session1', 'stim_type1:factual_truth1',
    'session1:factual_truth1', 'stim_type1:session1:factual_truth1'
  ),
  labels = PARAM_NAMES
)

# ----- Prior density at 0 -----
prior_sd <- 0.5
prior_density_0 <- dnorm(0, mean = 0, sd = prior_sd)

# Threshold for drawing posterior dot
dot_threshold <- 0.01

bf_plot <- fixed_effects %>%
  filter(parameter != PARAM_NAMES[1]) %>%
  ggplot(aes(x = post_sample)) +
  # posterior histogram
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "#969696",
                 color = "black", alpha = 0.7) +
  # posterior density
  geom_density(linewidth = 1, color = COLOR_PALETTE[1], n = 5000) +
  # prior density (dashed)
  stat_function(fun = dnorm, args = list(mean = 0, sd = prior_sd),
                xlim = c(-0.5, 0.5), linewidth = 1, linetype = "dashed") +
  # posterior dot at 0 (only if density > threshold)
  geom_point(data = fixed_effects %>% filter(y > dot_threshold),
             aes(x = 0, y = y),
             size = 3, shape = 21, fill = COLOR_PALETTE[2]) +
  # prior dot at 0
  geom_point(aes(x = 0, y = prior_density_0), size = 3, shape = 21, fill = COLOR_PALETTE[2]) +
  # Bayes factor text (fixed y = 10)
  geom_text(aes(label = bf, x = text_loc_x),
            y = 10,
            family = "serif",
            size = 5,
            alpha = 0.8) +
  # axis scaling
  scale_x_continuous(expand = c(0.05, 0.05), limits = c(-0.5, 0.5)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 16)) +
  facet_wrap(~parameter, labeller = label_parsed, ncol = 3) +
  ggthemes::theme_tufte() +
  labs(x = 'Posterior estimate', y = "Density") +
  theme(
    axis.line = element_line(linewidth = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_2),
    axis.text.y = element_text(size = FONT_SIZE_2),
    strip.text.x = element_text(size = FONT_SIZE_2),
    strip.text.y = element_text(size = FONT_SIZE_2, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.spacing = unit(1, "lines")
  )

ggsave(
  "../plots/01_glm_exp1.jpeg",
  bf_plot,
  dpi=300, width = 12, height = 8
)


post_summaries <- fixed_effects %>% 
  group_by(parameter) %>% 
  summarise(
    median = median(post_sample),
    q5 = quantile(post_sample, probs = 0.025),
    q95 = quantile(post_sample, probs = 0.975),
    .groups = "drop"
  )  %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(post_summaries, "../tables/hlm_post_summaries_exp1.csv")

#------------------------------------------------------------------------------#
# Retention interval: Experiment 2
#------------------------------------------------------------------------------#
df_exp2 <- read.csv('../data/data_exp_2.csv') %>% 
  mutate(
    stim_type = as_factor(stim_type),
    factual_truth = as_factor(factual_truth)
  )

contrasts(df_exp2$factual_truth) <- rbind(-1.0, 1.0)
contrasts(df_exp2$stim_type) <- rbind(-1.0, 1.0)

formula <- resp ~ stim_type * factual_truth + (1 + stim_type * factual_truth | id)

priors <- prior(normal(0, 1.0), class = b)
model_exp2 <- brm(
  formula = formula,
  data = df_exp2,
  family = bernoulli("logit"),
  prior = priors,
  autocor = priors,
  iter = 4000,
  cores = 4,
  chains = 4,
  file = "../fits/hlm_fit_exp2",
  sample_prior = "yes"
)

model_exp2

# ----- Hypothesis tests -----
hyp_rep <- "stim_type1 = 0"
bf_rep <- 1 / hypothesis(model_exp2, hyp_rep)$hypothesis["Evid.Ratio"]

hyp_factual_truth <- "factual_truth1 = 0"
bf_factual_truth <- 1 / hypothesis(model_exp2, hyp_factual_truth)$hypothesis["Evid.Ratio"]

hyp_type_truth <- "stim_type1:factual_truth1 = 0"
bf_type_truth <- 1 / hypothesis(model_exp2, hyp_type_truth)$hypothesis["Evid.Ratio"]

# ----- Extract posterior samples -----
fixed_effects <- brms::fixef(model_exp2, summary = FALSE) %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "post_sample") %>%
  mutate(parameter = factor(parameter))

# ----- Posterior density at 0 -----
fixed_effects$y <- NA
params <- unique(fixed_effects$parameter)
for (param in params) {
  post_samples <- fixed_effects$post_sample[fixed_effects$parameter == param]
  d <- density(post_samples)
  fixed_effects$y[fixed_effects$parameter == param] <- approx(d$x, d$y, xout = 0, rule = 2)$y
}

# Threshold for posterior dot
dot_threshold <- 0.01

# ----- Bayes factors -----
fixed_effects$bf <- NA
fixed_effects$bf[fixed_effects$parameter == "stim_type1"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "factual_truth1"] <- paste("BF =", round(bf_factual_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:factual_truth1"] <- paste("BF =", round(bf_type_truth, 2))

# ----- Text positions -----
fixed_effects$text_loc_x <- NA
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "factual_truth1"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:factual_truth1"] <- -0.3

# ----- Factor levels & labels -----
fixed_effects$parameter <- factor(
  fixed_effects$parameter,
  levels = c('Intercept', 'stim_type1', 'factual_truth1', 'stim_type1:factual_truth1'),
  labels = PARAM_NAMES[c(1,3,4,6)]
)

# ----- Prior density at 0 -----
prior_sd <- 1.0
prior_density_0 <- dnorm(0, mean = 0, sd = prior_sd)

# ----- Plot -----
bf_plot <- fixed_effects %>%
  filter(parameter != PARAM_NAMES[1]) %>%
  ggplot(aes(x = post_sample)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "#969696", color = "black", alpha = 0.7) +
  geom_density(linewidth = 1, color = COLOR_PALETTE[1], n = 5000) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = prior_sd),
                xlim = c(-1.5, 1.5), linewidth = 1, linetype = "dashed") +
  geom_point(data = fixed_effects %>% filter(y > dot_threshold),
             aes(x = 0, y = y), size = 3, shape = 21, fill = COLOR_PALETTE[2]) +
  geom_point(aes(x = 0, y = prior_density_0), size = 3, shape = 21, fill = COLOR_PALETTE[2]) +
  geom_text(aes(label = bf, x = text_loc_x),
            y = 10, family = "serif", size = 5, alpha = 0.8) +
  scale_x_continuous(expand = c(0.05, 0.05), limits = c(-0.5, 0.5)) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,16)) +
  facet_wrap(~parameter, labeller = label_parsed, ncol = 3) +
  ggthemes::theme_tufte() +
  labs(x = "Posterior estimate", y = "Density") +
  theme(
    axis.line = element_line(linewidth = 0.5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_3),
    strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    panel.grid.major = element_line(color = alpha("gray70",0.3)),
    panel.grid.minor = element_line(color = alpha("gray70",0.15)),
    panel.background = element_blank(),
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10)),
    panel.spacing = unit(1,"lines")
  )

ggsave("../plots/01_glm_exp2.jpeg", bf_plot, dpi = 300, width = 12, height = 3)

# ----- Posterior summaries -----
post_summaries <- fixed_effects %>%
  group_by(parameter) %>%
  summarise(
    median = median(post_sample),
    q5 = quantile(post_sample, probs = 0.025),
    q95 = quantile(post_sample, probs = 0.975),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x,2)))

write_csv(post_summaries, "../tables/hlm_post_summaries_exp2.csv")

#------------------------------------------------------------------------------#
# Reading time: Experiment 2
#------------------------------------------------------------------------------#
formula <- rt_question ~ stim_type * factual_truth + (1 + stim_type * factual_truth | id)

priors <- prior(normal(0, 1.0), class = b)
model_reading_time_exp2 <- brm(
  formula = formula,
  family = lognormal(),
  data = df_exp2,
  prior = priors,
  iter = 4000,
  cores = 4,
  chains = 4,
  file = "../fits/hlm_fit_reading_time_exp2",
  sample_prior = "yes"
)

model_reading_time_exp2

hyp_rep <- "stim_type1 = 0"
test_rep <- hypothesis(model_reading_time_exp2, hyp_rep)$hypothesis
bf_rep <- round(1 / test_rep["Evid.Ratio"], 2)

hyp_factual_truth <- "factual_truth1 = 0"
test_factual_truth <- hypothesis(model_reading_time_exp2, hyp_factual_truth)$hypothesis
bf_factual_truth <- round(1 / test_factual_truth["Evid.Ratio"], 2)

hyp_type_truth <- "stim_type1:factual_truth1 = 0"
test_type_truth <- hypothesis(model_reading_time_exp2, hyp_type_truth)$hypothesis
bf_type_truth <- round(1 / test_type_truth["Evid.Ratio"], 2)

