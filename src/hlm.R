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
# Experiment 1
#------------------------------------------------------------------------------#
df_session_1 <- read.csv('../data/data_session_1.csv')
df_session_2 <- read.csv('../data/data_session_2.csv')

df_exp1 <- rbind(df_session_1, df_session_2) %>% 
  mutate(
    session = as_factor(session),
    stim_type = as_factor(stim_type),
    factual_truth = as_factor(factual_truth)
  )

contrasts(df_exp1$session) <- contr.sum(nlevels(df_exp1$session)) / 2
colnames(contrasts(df_exp1$session)) <- levels(df_exp1$session)[2]
contrasts(df_exp1$factual_truth) <- contr.sum(nlevels(df_exp1$factual_truth)) / 2

formula <- resp ~ stim_type * session * factual_truth + (1 + stim_type * session * factual_truth | id)

priors <- prior(normal(0, 1.0), class = b)
model_exp1 <- brm(
  formula = formula,
  data = df_exp1,
  family = bernoulli("logit"),
  prior = priors,
  autocor = priors,
  iter = 4000,
  cores = 4,
  chains = 4,
  file = "../fits/hlm_fit_exp1",
  sample_prior = "yes"
)

model_exp1

# hypothesis test
hyp_session <- "session2 = 0"
test_session <- hypothesis(model_exp1, hyp_session)$hypothesis
bf_session <- 1 / test_session["Evid.Ratio"]

hyp_rep <- "stim_type1 = 0"
test_rep <- hypothesis(model_exp1, hyp_rep)$hypothesis
bf_rep <- 1 / test_rep["Evid.Ratio"]

hyp_factual_truth <- "factual_truth1 = 0"
test_factual_truth <- hypothesis(model_exp1, hyp_factual_truth)$hypothesis
bf_factual_truth <- 1 / test_factual_truth["Evid.Ratio"]

hyp_type_session <- "stim_type1:session2 = 0"
test_type_session <- hypothesis(model_exp1, hyp_type_session)$hypothesis
bf_type_session <- 1 / test_type_session["Evid.Ratio"]

hyp_type_truth <- "stim_type1:factual_truth1 = 0"
test_type_truth <- hypothesis(model_exp1, hyp_type_truth)$hypothesis
bf_type_truth <- 1 / test_type_truth["Evid.Ratio"]

hyp_session_truth <- "session2:factual_truth1 = 0"
test_session_truth <- hypothesis(model_exp1, hyp_session_truth)$hypothesis
bf_session_truth <- 1 / test_session_truth["Evid.Ratio"]

hyp_three_way <- "stim_type1:session2:factual_truth1 = 0"
test_three_way <- hypothesis(model_exp1, hyp_three_way)$hypothesis
bf_three_way <- 1 / test_three_way["Evid.Ratio"]

fixed_effects <- brms::fixef(model_exp1, summary = FALSE) %>% 
  as_data_frame()

fixed_effects <- fixed_effects %>% 
  pivot_longer(everything(), names_to = "parameter", values_to = "post_sample") %>% 
  mutate(parameter = factor(parameter))

fixed_effects$y <- NA
params <- unique(fixed_effects$parameter)
for (param in params) {
  post_samples <- fixed_effects$post_sample[fixed_effects$parameter == param]
  d <- density(post_samples)
  posterior_density_0 <- approx(d$x, d$y, xout = 0)$y
  fixed_effects$y[fixed_effects$parameter == param] <- posterior_density_0
}
fixed_effects$y <- ifelse(fixed_effects$y > 6, NA, fixed_effects$y)

fixed_effects$bf <- NA
fixed_effects$bf[fixed_effects$parameter == "session2"] <- paste("BF =", round(bf_session, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "factual_truth1"] <- paste("BF =", round(bf_factual_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:session2"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "stim_type1:factual_truth1"] <- paste("BF =", round(bf_type_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "session2:factual_truth1"] <- paste("BF =", round(bf_session_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:session2:factual_truth1"] <- paste("BF =", round(bf_three_way, 2))

# set text position
fixed_effects$text_loc_x <- NA
fixed_effects$text_loc_x[fixed_effects$parameter == "session2"] <- -0.8
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1"] <- -0.5
fixed_effects$text_loc_x[fixed_effects$parameter == "factual_truth1"] <- -0.5
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:session2"] <- -0.3
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:factual_truth1"] <- -0.8
fixed_effects$text_loc_x[fixed_effects$parameter == "session2:factual_truth1"] <- -0.9
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:session2:factual_truth1"] <- -0.8

fixed_effects$parameter <- factor(
  fixed_effects$parameter,
  levels = c(
    'Intercept', 'session2', 'stim_type1', 'factual_truth1',
    'stim_type1:session2', 'stim_type1:factual_truth1',
    'session2:factual_truth1', 'stim_type1:session2:factual_truth1'
  ),
  labels = PARAM_NAMES
)

# plotting
bf_plot <- fixed_effects %>%
  filter(
    parameter != PARAM_NAMES[1]
  ) %>% 
  ggplot(aes(x=post_sample)) +
  geom_histogram(
    aes(y=after_stat(density)),
    bins = 30, fill = "#969696",
    color = "black", alpha=0.7
  ) +
  geom_density(linewidth = 1, color = COLOR_PALETTE[1],
               n = 5000) +
  stat_function(
    fun = dnorm, args = list(mean = 0, sd = .5),
    xlim = c(-1.5, 1.5), linewidth = 1,
    n = 500, linetype = "dashed"
  ) +
  geom_point(
    aes(x = 0, y = y), size = 3,
    shape = 21, fill = COLOR_PALETTE[2]
  ) +
  geom_point(
    aes(x = 0, y = 0.79788), size = 3,
    shape = 21, fill = COLOR_PALETTE[2]
  ) +
  geom_text(
    aes(label = bf, x = text_loc_x),
    y = 2.0,
    family = "serif",
    size = 5,
    alpha = 0.8
  ) +
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  facet_wrap(
    ~parameter,
    labeller = label_parsed,
    ncol = 3
  ) +
  ggthemes::theme_tufte() +
  labs(
    x = 'Posterior estimate',
    y = "Density",
    title = 'Group-level Parameter Estimates: Experiment 1'
  ) +
  theme(
    axis.line = element_line(size = .5, color = "#969696"),
      axis.ticks = element_line(color = "#969696"),
      axis.text.x = element_text(size = FONT_SIZE_3),
      axis.text.y = element_text(size = FONT_SIZE_3),
      strip.text.x = element_text(size = FONT_SIZE_3),
      strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
      text = element_text(size = FONT_SIZE_2),
      plot.title = element_text(
        size = FONT_SIZE_1, hjust = 0.5, face = 'bold'
      ),
      panel.grid.major = element_line(color = alpha("gray70", 0.3)),
      panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
      panel.background = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      panel.spacing = unit(1, "lines")
  )

bf_plot

ggsave(
  "../plots/bf_plot_exp1.pdf",
  dpi=300, width = 10, height = 8
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

write_csv(post_summaries, "../data/post_summaries_hlm_exp1.csv")

#------------------------------------------------------------------------------#
# Experiment 2
#------------------------------------------------------------------------------#
df_exp2 <- read.csv('../data/data_exp_2.csv') %>% 
  mutate(
    stim_type = as_factor(stim_type),
    factual_truth = as_factor(factual_truth)
  )

contrasts(df_exp2$factual_truth) <- contr.sum(nlevels(df_exp2$factual_truth)) / 2

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

# hypothesis test
hyp_rep <- "stim_type1 = 0"
test_rep <- hypothesis(model_exp2, hyp_rep)$hypothesis
bf_rep <- 1 / test_rep["Evid.Ratio"]

hyp_factual_truth <- "factual_truth1 = 0"
test_factual_truth <- hypothesis(model_exp2, hyp_factual_truth)$hypothesis
bf_factual_truth <- 1 / test_factual_truth["Evid.Ratio"]

hyp_type_truth <- "stim_type1:factual_truth1 = 0"
test_type_truth <- hypothesis(model_exp2, hyp_type_truth)$hypothesis
bf_type_truth <- 1 / test_type_truth["Evid.Ratio"]

fixed_effects <- brms::fixef(model_exp2, summary = FALSE) %>% 
  as_data_frame()

fixed_effects <- fixed_effects %>% 
  pivot_longer(everything(), names_to = "parameter", values_to = "post_sample") %>% 
  mutate(parameter = factor(parameter))

fixed_effects$y <- NA
params <- unique(fixed_effects$parameter)
for (param in params) {
  post_samples <- fixed_effects$post_sample[fixed_effects$parameter == param]
  d <- density(post_samples)
  posterior_density_0 <- approx(d$x, d$y, xout = 0)$y
  fixed_effects$y[fixed_effects$parameter == param] <- posterior_density_0
}
fixed_effects$y <- ifelse(fixed_effects$y > 6, NA, fixed_effects$y)

fixed_effects$bf <- NA
fixed_effects$bf[fixed_effects$parameter == "stim_type1"] <- paste("BF > 1000")
fixed_effects$bf[fixed_effects$parameter == "factual_truth1"] <- paste("BF =", round(bf_factual_truth, 2))
fixed_effects$bf[fixed_effects$parameter == "stim_type1:factual_truth1"] <- paste("BF =", round(bf_type_truth, 2))

# set text position
fixed_effects$text_loc_x <- NA
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1"] <- -0.5
fixed_effects$text_loc_x[fixed_effects$parameter == "factual_truth1"] <- -0.7
fixed_effects$text_loc_x[fixed_effects$parameter == "stim_type1:factual_truth1"] <- -0.8

fixed_effects$parameter <- factor(
  fixed_effects$parameter,
  levels = c(
    'Intercept', 'stim_type1', 'factual_truth1', 'stim_type1:factual_truth1'
  ),
  labels = PARAM_NAMES[c(1, 3, 4, 6)]
)

# plotting
bf_plot <- fixed_effects %>%
  filter(
    parameter != PARAM_NAMES[1]
  ) %>% 
  ggplot(aes(x=post_sample)) +
  geom_histogram(
    aes(y=after_stat(density)),
    bins = 30, fill = "#969696",
    color = "black", alpha=0.7
  ) +
  geom_density(linewidth = 1, color = COLOR_PALETTE[1],
               n = 5000) +
  stat_function(
    fun = dnorm, args = list(mean = 0, sd = .5),
    xlim = c(-1.5, 1.5), linewidth = 1,
    n = 500, linetype = "dashed"
  ) +
  geom_point(
    aes(x = 0, y = y), size = 3,
    shape = 21, fill = COLOR_PALETTE[2]
  ) +
  geom_point(
    aes(x = 0, y = 0.79788), size = 3,
    shape = 21, fill = COLOR_PALETTE[2]
  ) +
  geom_text(
    aes(label = bf, x = text_loc_x),
    y = 2.0,
    family = "serif",
    size = 5,
    alpha = 0.8
  ) +
  scale_x_continuous(expand = c(0.05, 0.05)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  facet_wrap(
    ~parameter,
    labeller = label_parsed,
    ncol = 3
  ) +
  ggthemes::theme_tufte() +
  labs(
    x = 'Posterior estimate',
    y = "Density",
    title = 'Group-level Parameter Estimates: Experiment 2'
  ) +
  theme(
    axis.line = element_line(size = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_3),
    strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    plot.title = element_text(
      size = FONT_SIZE_1, hjust = 0.5, face = 'bold'
    ),
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.spacing = unit(1, "lines")
  )

bf_plot

ggsave(
  "../plots/bf_plot_exp2.pdf",
  bf_plot,
  dpi=300, width = 10, height = 4
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

write_csv(post_summaries, "../data/post_summaries_hlm_exp2.csv")


# Model on reading time
formula <- rt_question ~ stim_type * factual_truth + (1 + stim_type * factual_truth | id)

priors <- prior(normal(0, 1.0), class = b)
model_reading_time_exp2 <- brm(
  formula = formula,
  family = lognormal(),
  data = df_exp2,
  prior = priors,
  iter = 10000,
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

