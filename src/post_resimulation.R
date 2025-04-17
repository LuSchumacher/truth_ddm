library(tidyverse)
library(magrittr)
library(bayesplot)
library(tidybayes)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ddm_simulator.R")

df_s1 <- read_csv('../data/data_session_1.csv')
df_s2 <- read_csv('../data/data_session_2.csv')
df_exp2 <- read_csv('../data/data_exp_2.csv')

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
COLOR_PALETTE <- c('#27374D', '#B70404')

NUM_RESIMS <- 200
idx <- sample(1:8000, NUM_RESIMS, replace = FALSE)

#------------------------------------------------------------------------------#
# Re-simulation: Session 1
#------------------------------------------------------------------------------#
fit_session_1 <- readRDS("../fits/fit_session_1")
df_posterior_s1 <- as_tibble(
  fit_session_1$draws(
    inc_warmup = FALSE,
    format = "draws_matrix"
  )
) %>% 
  mutate(draw = row_number()) %>% 
  select(
    draw,
    matches("^v\\[[1-4],\\s?\\d+\\]$"),
    matches("^a\\[[1-4],\\s?\\d+\\]$"),
    matches("^ndt\\[[1-4],\\s?\\d+\\]$"),
    matches("^bias\\[[1-4],\\s?\\d+\\]$")
  ) %>% 
  pivot_longer(
    cols = -draw,
    names_to = "param_full",
    values_to = "value"
  ) %>% 
  mutate(
    parameter = str_extract(param_full, "^[a-z]+"),  # Capture parameter name (v, a, etc.)
    condition = as.integer(str_extract(param_full, "(?<=\\[)\\d+")),  # Extract condition (1-4)
    participant = as.integer(str_extract(param_full, "(?<=,\\s?)\\d+")),  # Extract participant number
    value = as.numeric(value)
  ) %>%
  select(draw, parameter, condition, participant, value)


pred_data_s1 <- tibble()
for (i in 1:length(idx)) {
  tmp_post <- df_posterior_s1 %>% 
    filter(draw == idx[i])
  tmp_pred_data <- tibble()
  for (j in 1:nrow(df_s1)) {
    v <- tmp_post$value[
      tmp_post$parameter == "v" &
        tmp_post$participant == df_s1$id[j] &
        tmp_post$condition == df_s1$condition[j]
    ]
    a <- tmp_post$value[
      tmp_post$parameter == "a" &
        tmp_post$participant == df_s1$id[j] &
        tmp_post$condition == df_s1$condition[j]
    ]
    ndt <- tmp_post$value[
      tmp_post$parameter == "ndt" &
        tmp_post$participant == df_s1$id[j] &
        tmp_post$condition == df_s1$condition[j]
    ]
    bias <- tmp_post$value[
      tmp_post$parameter == "bias" &
        tmp_post$participant == df_s1$id[j] &
        tmp_post$condition == df_s1$condition[j]
    ]
    x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
    names(x) <- c("resp", "rt")
    tmp_pred_data <- tmp_pred_data %>% 
      bind_rows(x)
  }
  tmp_pred_data$id <- df_s1$id
  tmp_pred_data$condition <- df_s1$condition
  tmp_pred_data$resim_id <- i
  tmp_pred_data$stim_type <- df_s1$stim_type
  tmp_pred_data$factual_truth <- df_s1$factual_truth
  pred_data_s1 <- pred_data_s1 %>% 
    bind_rows(tmp_pred_data)
  print(paste("Resimulation Nr.", i, "finished"))
}

write_csv(pred_data_s1, "../data/pred_data_s1.csv")

emp_rt_summaries_s1 <- df_s1 %>%
  group_by(stim_type) %>%
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  mutate(
    source = "Observed",
    session = 1,
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )

pred_rt_summaries_s1 <- pred_data_s1 %>% 
  group_by(resim_id, stim_type) %>% 
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  group_by(quantile, stim_type) %>% 
  summarise(
    mean = median(rt),
    q_lower = as.numeric(quantile(rt, probs = 0.025)),
    q_upper = as.numeric(quantile(rt, probs = 0.975)),
    .groups = "drop"
  ) %>% 
  mutate(
    source = "Re-simulated",
    session = 1,
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )

#------------------------------------------------------------------------------#
# Re-simulation: Session 2
#------------------------------------------------------------------------------#
fit_session_2 <- readRDS("../fits/fit_session_2")
df_posterior_s2 <- as_tibble(
  fit_session_2$draws(
    inc_warmup = FALSE,
    format = "draws_matrix"
  )
) %>% 
  mutate(draw = row_number()) %>% 
  select(
    draw,
    matches("^v\\[[1-4],\\s?\\d+\\]$"),
    matches("^a\\[[1-4],\\s?\\d+\\]$"),
    matches("^ndt\\[[1-4],\\s?\\d+\\]$"),
    matches("^bias\\[[1-4],\\s?\\d+\\]$")
  ) %>% 
  pivot_longer(
    cols = -draw,
    names_to = "param_full",
    values_to = "value"
  ) %>% 
  mutate(
    parameter = str_extract(param_full, "^[a-z]+"),  # Capture parameter name (v, a, etc.)
    condition = as.integer(str_extract(param_full, "(?<=\\[)\\d+")),  # Extract condition (1-4)
    participant = as.integer(str_extract(param_full, "(?<=,\\s?)\\d+")),  # Extract participant number
    value = as.numeric(value)
  ) %>%
  select(draw, parameter, condition, participant, value)

pred_data_s2 <- tibble()
for (i in 1:length(idx)) {
  tmp_post <- df_posterior_s2 %>% 
    filter(draw == idx[i])
  tmp_pred_data <- tibble()
  for (j in 1:nrow(df_s2)) {
    v <- tmp_post$value[
      tmp_post$parameter == "v" &
        tmp_post$participant == df_s2$id[j] &
        tmp_post$condition == df_s2$condition[j]
    ]
    a <- tmp_post$value[
      tmp_post$parameter == "a" &
        tmp_post$participant == df_s2$id[j] &
        tmp_post$condition == df_s2$condition[j]
    ]
    ndt <- tmp_post$value[
      tmp_post$parameter == "ndt" &
        tmp_post$participant == df_s2$id[j] &
        tmp_post$condition == df_s2$condition[j]
    ]
    bias <- tmp_post$value[
      tmp_post$parameter == "bias" &
        tmp_post$participant == df_s2$id[j] &
        tmp_post$condition == df_s2$condition[j]
    ]
    x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
    names(x) <- c("resp", "rt")
    tmp_pred_data <- tmp_pred_data %>% 
      bind_rows(x)
  }
  tmp_pred_data$id <- df_s2$id
  tmp_pred_data$condition <- df_s2$condition
  tmp_pred_data$resim_id <- i
  tmp_pred_data$stim_type <- df_s2$stim_type
  tmp_pred_data$factual_truth <- df_s2$factual_truth
  pred_data_s2 <- pred_data_s2 %>% 
    bind_rows(tmp_pred_data)
  print(paste("Resimulation Nr.", i, "finished"))
}

write_csv(pred_data_s2, "../data/pred_data_s2.csv")

emp_rt_summaries_s2 <- df_s2 %>%
  group_by(stim_type) %>%
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  mutate(
    source = "Observed",
    session = 2,
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )

pred_rt_summaries_s2 <- pred_data_s2 %>% 
  group_by(resim_id, stim_type) %>% 
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  group_by(quantile, stim_type) %>% 
  summarise(
    mean = median(rt),
    q_lower = as.numeric(quantile(rt, probs = 0.025)),
    q_upper = as.numeric(quantile(rt, probs = 0.975)),
    .groups = "drop"
  ) %>% 
  mutate(
    source = "Re-simulated",
    session = 2,
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )

#------------------------------------------------------------------------------#
# Re-simulation: Experiment 2
#------------------------------------------------------------------------------#
fit_exp_2 <- readRDS("../fits/fit_exp_2")
df_posterior_exp2 <- as_tibble(
  fit_exp_2$draws(
    inc_warmup = FALSE,
    format = "draws_matrix"
  )
) %>% 
  mutate(draw = row_number()) %>% 
  select(
    draw,
    matches("^v\\[[1-4],\\s?\\d+\\]$"),
    matches("^a\\[[1-4],\\s?\\d+\\]$"),
    matches("^ndt\\[[1-4],\\s?\\d+\\]$"),
    matches("^bias\\[[1-4],\\s?\\d+\\]$")
  ) %>% 
  pivot_longer(
    cols = -draw,
    names_to = "param_full",
    values_to = "value"
  ) %>% 
  mutate(
    parameter = str_extract(param_full, "^[a-z]+"),  # Capture parameter name (v, a, etc.)
    condition = as.integer(str_extract(param_full, "(?<=\\[)\\d+")),  # Extract condition (1-4)
    participant = as.integer(str_extract(param_full, "(?<=,\\s?)\\d+")),  # Extract participant number
    value = as.numeric(value)
  ) %>%
  select(draw, parameter, condition, participant, value)

pred_data_exp2 <- tibble()
for (i in 1:length(idx)) {
  tmp_post <- df_posterior_exp2 %>% 
    filter(draw == idx[i])
  tmp_pred_data <- tibble()
  for (j in 1:nrow(df_exp2)) {
    v <- tmp_post$value[
      tmp_post$parameter == "v" &
        tmp_post$participant == df_exp2$id[j] &
        tmp_post$condition == df_exp2$condition[j]
    ]
    a <- tmp_post$value[
      tmp_post$parameter == "a" &
        tmp_post$participant == df_exp2$id[j] &
        tmp_post$condition == df_exp2$condition[j]
    ]
    ndt <- tmp_post$value[
      tmp_post$parameter == "ndt" &
        tmp_post$participant == df_exp2$id[j] &
        tmp_post$condition == df_exp2$condition[j]
    ]
    bias <- tmp_post$value[
      tmp_post$parameter == "bias" &
        tmp_post$participant == df_exp2$id[j] &
        tmp_post$condition == df_exp2$condition[j]
    ]
    x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
    names(x) <- c("resp", "rt")
    tmp_pred_data <- tmp_pred_data %>% 
      bind_rows(x)
  }
  tmp_pred_data$id <- df_exp2$id
  tmp_pred_data$condition <- df_exp2$condition
  tmp_pred_data$resim_id <- i
  tmp_pred_data$stim_type <- df_exp2$stim_type
  tmp_pred_data$factual_truth <- df_exp2$factual_truth
  pred_data_exp2 <- pred_data_exp2 %>% 
    bind_rows(tmp_pred_data)
  print(paste("Resimulation Nr.", i, "finished"))
}
  
write_csv(pred_data_exp2, "../data/pred_data_exp2")
  
emp_rt_summaries_exp2 <- df_exp2 %>%
  group_by(stim_type) %>%
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  mutate(
    source = "Observed",
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )

pred_rt_summaries_exp2 <- pred_data_exp2 %>% 
  group_by(resim_id, stim_type) %>% 
  reframe(
    quantile = seq(0.1, 0.9, by = 0.2),
    rt = quantile(rt, probs = quantile)
  ) %>% 
  group_by(quantile, stim_type) %>% 
  summarise(
    mean = median(rt),
    q_lower = as.numeric(quantile(rt, probs = 0.025)),
    q_upper = as.numeric(quantile(rt, probs = 0.975)),
    .groups = "drop"
  ) %>% 
  mutate(
    source = "Re-simulated",
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements")
  )
  
#------------------------------------------------------------------------------#
# RT quantiles plot: Experiment 1
#------------------------------------------------------------------------------#
emp_rt_summaries <- bind_rows(emp_rt_summaries_s1, emp_rt_summaries_s2)
emp_rt_summaries %<>%
  mutate(session = ifelse(
    session == 1, 'Exp. session 1\n(10 minute interval)', 'Exp. session 2\n(1 week interval)')
  )
pred_rt_summaries <- bind_rows(pred_rt_summaries_s1, pred_rt_summaries_s2)
pred_rt_summaries %<>%
  mutate(session = ifelse(
    session == 1, 'Exp. session 1\n(10 minute interval)', 'Exp. session 2\n(1 week interval)')
  )

pred_rt_summaries %>% 
  ggplot(aes(x = as.factor(quantile), y = mean)) +
  geom_pointrange(
    aes(
      ymin = q_lower,
      ymax = q_upper,
      color = "Re-simulated",
    ),
    size = 0.4, alpha = 0.75
  ) +
  geom_point(
    data = emp_rt_summaries,
    aes(x = as.factor(quantile), y = rt, color = "Observed"),
    size = 1.6, alpha = 0.75
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed" = COLOR_PALETTE[2],
      "Re-simulated" = COLOR_PALETTE[1]
    )
  ) +
  facet_grid(stim_type ~ session) +
  labs(
    x = "Quantile",
    y = "Response time (s)",
    title = "Observed and Re-simulated Response Time Quantiles: Experiment 1"
  ) +
  ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(size = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_3),
    strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    plot.title = element_text(size = FONT_SIZE_1,
                              hjust = 0.5,
                              face = 'bold'),
    legend.position = "right", # Places legend outside
    legend.box = "vertical",
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank()
  )

ggsave(
  '../plots/03_rt_quantiles_plot_exp1.pdf',
  device = 'pdf', dpi = 300,
  width = 12, height = 6
)

#------------------------------------------------------------------------------#
# RT quantiles plot: Experiment 2
#------------------------------------------------------------------------------#
pred_rt_summaries_exp2 %>% 
  ggplot(aes(x = as.factor(quantile), y = mean)) +
  geom_pointrange(
    aes(
      ymin = q_lower,
      ymax = q_upper,
      color = "Re-simulated",
    ),
    size = 0.4, alpha = 0.75
  ) +
  geom_point(
    data = emp_rt_summaries_exp2,
    aes(x = as.factor(quantile), y = rt, color = "Observed"),
    size = 1.6, alpha = 0.75
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed" = COLOR_PALETTE[2],
      "Re-simulated" = COLOR_PALETTE[1]
    )
  ) +
  facet_wrap(~stim_type) +
  labs(
    x = "Quantile",
    y = "Response time (s)",
    title = "Observed and Re-simulated Response Time Quantiles: Experiment 2"
  ) +
  ggthemes::theme_tufte() +
  theme(
    axis.line = element_line(size = .5, color = "#969696"),
    axis.ticks = element_line(color = "#969696"),
    axis.text.x = element_text(size = FONT_SIZE_3),
    axis.text.y = element_text(size = FONT_SIZE_3),
    strip.text.x = element_text(size = FONT_SIZE_3),
    strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
    text = element_text(size = FONT_SIZE_2),
    plot.title = element_text(size = FONT_SIZE_1,
                              hjust = 0.5,
                              face = 'bold'),
    legend.position = "right", # Places legend outside
    legend.box = "vertical",
    panel.grid.major = element_line(color = alpha("gray70", 0.3)),
    panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
    panel.background = element_blank()
  )

ggsave(
  '../plots/03_rt_quantiles_plot_exp2.pdf',
  device = 'pdf', dpi = 300,
  width = 12, height = 6
)

#------------------------------------------------------------------------------#
# Response probability plot: Experiment 1
#------------------------------------------------------------------------------#
pred_data_s1 %<>%
  mutate(session = 1)
pred_data_s2 %<>%
  mutate(session = 2)
pred_data_exp1 <- bind_rows(
  pred_data_s1, pred_data_s2
)
emp_data_exp1 <- bind_rows(df_s1, df_s2)

indiviudal_emp_resp_prob <- emp_data_exp1 %>% 
  group_by(
    id, session, stim_type
  ) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  mutate(
    session = ifelse(session == 1, 'Exp. session 1\n(10 minute interval)', 'Exp. session 2\n(1 week interval)'),
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Observed'
  )

pred_resp_prob <- pred_data_exp1 %>% 
  group_by(
    id, session, stim_type
  ) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  group_by(session, stim_type) %>% 
  summarise(
    lower_ci = quantile(resp_prob, 0.025),
    upper_ci = quantile(resp_prob, 0.975),
    resp_prob = mean(resp_prob),
    .groups = "drop"
  ) %>% 
  mutate(
    session = ifelse(session == 1, 'Exp. session 1\n(10 minute interval)', 'Exp. session 2\n(1 week interval)'),
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Re-simulated'
  )

group_emp_resp_prob <- emp_data_exp1 %>% 
  group_by(id, session, stim_type) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  group_by(session, stim_type) %>% 
  summarise(
    sd_resp_prob = sd(resp_prob),
    resp_prob = mean(resp_prob),
    lower_ci = resp_prob - sd_resp_prob,
    upper_ci = resp_prob + sd_resp_prob,
    .groups = "drop"
  ) %>% 
  mutate(
    session = ifelse(session == 1, 'Exp. session 1\n(10 minute interval)', 'Exp. session 2\n(1 week interval)'),
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Observed'
  )

summary_df <- bind_rows(pred_resp_prob, group_emp_resp_prob)

resp_prob_plot <- summary_df %>% 
  ggplot(aes(x = stim_type,
             y = resp_prob,
             color = type)) +
  # scatter for individuals
  geom_jitter(aes(x = stim_type, y = resp_prob),
              color = COLOR_PALETTE[1],
              data = indiviudal_emp_resp_prob,
              width = 0.1,
              size = 2,
              alpha = 0.4) + 
  geom_pointinterval(aes(x = stim_type, y = resp_prob,
                         ymin = lower_ci,
                         ymax = upper_ci),
                     size = 10,
                     alpha = 1.0,
                     position = position_dodge(width = 0.75)) +
  facet_grid(~session) + 
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0.0, 0.0)) +
  scale_color_manual(values = COLOR_PALETTE) +
  scale_fill_manual(values = COLOR_PALETTE) +
  ggthemes::theme_tufte() + 
  ylab("Probability for a\n\"true\" response")+
  xlab("Repetition status")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text.x = element_text(size = FONT_SIZE_3),
        strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid.major = element_line(color = alpha("gray70", 0.3)),
        panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
        panel.background = element_blank()
  ) + 
  ggtitle('Empirical and Re-simulated Response Probabilities: Experiment 1') + 
  guides(fill = guide_legend(title=""),
         color = guide_legend(title=""))

ggsave(
  '../plots/04_resp_prob_plot_exp1.pdf',
  resp_prob_plot,
  device = 'pdf', dpi = 300,
  width = 12, height = 6
)

#------------------------------------------------------------------------------#
# Response probability plot: Experiment 2
#------------------------------------------------------------------------------#
indiviudal_emp_resp_prob_exp2 <- df_exp2 %>% 
  group_by(
    id, stim_type
  ) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  mutate(
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Observed'
  )

pred_resp_prob_exp2 <- pred_data_exp2 %>% 
  group_by(
    id, stim_type
  ) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  group_by(stim_type) %>% 
  summarise(
    lower_ci = quantile(resp_prob, 0.025),
    upper_ci = quantile(resp_prob, 0.975),
    resp_prob = mean(resp_prob),
    .groups = "drop"
  ) %>% 
  mutate(
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Re-simulated'
  )

group_emp_resp_prob_exp2 <- df_exp2 %>% 
  group_by(id, stim_type) %>% 
  summarise(resp_prob = mean(resp), .groups = "drop") %>% 
  group_by(stim_type) %>% 
  summarise(
    sd_resp_prob = sd(resp_prob),
    resp_prob = mean(resp_prob),
    lower_ci = resp_prob - sd_resp_prob,
    upper_ci = resp_prob + sd_resp_prob,
    .groups = "drop"
  ) %>% 
  mutate(
    stim_type = ifelse(stim_type == 0, "New\n statements", "Repeated\n statements"),
    type = 'Observed'
  )

summary_df_exp2 <- bind_rows(pred_resp_prob_exp2, group_emp_resp_prob_exp2)

resp_prob_plot_exp2 <- summary_df_exp2 %>% 
  ggplot(aes(x = stim_type,
             y = resp_prob,
             color = type)) +
  # scatter for individuals
  geom_jitter(aes(x = stim_type, y = resp_prob),
              color = COLOR_PALETTE[1],
              data = indiviudal_emp_resp_prob_exp2,
              width = 0.1,
              size = 2,
              alpha = 0.4) + 
  geom_pointinterval(aes(x = stim_type, y = resp_prob,
                         ymin = lower_ci,
                         ymax = upper_ci),
                     size = 10,
                     alpha = 1.0,
                     position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0.0, 0.0)) +
  scale_color_manual(values = COLOR_PALETTE) +
  scale_fill_manual(values = COLOR_PALETTE) +
  ggthemes::theme_tufte() + 
  ylab("Probability for a\n\"true\" response")+
  xlab("Repetition status")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text.x = element_text(size = FONT_SIZE_3),
        strip.text.y = element_text(size = FONT_SIZE_3, angle = 0),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid.major = element_line(color = alpha("gray70", 0.3)),
        panel.grid.minor = element_line(color = alpha("gray70", 0.15)),
        panel.background = element_blank()
  ) + 
  ggtitle('Empirical and Re-simulated Response Probabilities: Experiment 2\n') + 
  guides(fill = guide_legend(title=""),
         color = guide_legend(title=""))

ggsave(
  '../plots/04_resp_prob_plot_exp2.pdf',
  resp_prob_plot_exp2,
  device = 'pdf', dpi = 300,
  width = 12, height = 6
)
