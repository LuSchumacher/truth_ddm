library(tidyverse)
library(cmdstanr)
library(posterior)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("0_ddm_simulator.R")

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
SMALLER_FONT <- 4
COLOR_PALETTE <- c('#27374D', '#7A8F7A', '#B70404')

NUM_RESIMS <- 100
NUM_POST_SAMPLES <- 12000

fit_session_1 <- readRDS("../fits/fit_session_1.rds")
draws_session_1 <- as_draws_matrix(fit_session_1)

df_session_1 <- read.csv('../data/data_session_1.csv')
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == 1, -1, 1)
df_session_1$factual_truth <- ifelse(df_session_1$factual_truth == -1, 1, 2)

# ---------------------------------------------------------------------------- #
# HELPER FUNCTION
# ---------------------------------------------------------------------------- #
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
      v    <- params$v_collapsed[truth, id, r]
      a    <- params$a_collapsed[truth, id, r]
      bias <- params$bias[cond, id, r]
      ndt_base <- params$ndt_collapsed[truth, id, r]
      ndt_var  <- params$ndt_var_collapsed[truth, id, r]
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
      a    <- params$a_collapsed[truth, id, r]
      bias <- params$bias_collapsed[truth, id, r]
      ndt_base <- params$ndt_collapsed[truth, id, r]
      ndt_var  <- params$ndt_var_collapsed[truth, id, r]
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

# ---------------------------------------------------------------------------- #
# POSTERIOR PREPARATION
# ---------------------------------------------------------------------------- #
idx <- sample(1:NUM_POST_SAMPLES, NUM_RESIMS, replace = FALSE)

post_session_1 <- extract_all_params(fit_session_1, df_session_1, idx)

post_session_1$v_collapsed <- array(NA, dim = c(2, 75, 100))
post_session_1$v_collapsed[1, , ] <- (post_session_1$v[1, , ] +
                                        post_session_1$v[2, , ]) / 2
post_session_1$v_collapsed[2, , ] <- (post_session_1$v[3, , ] +
                                        post_session_1$v[4, , ]) / 2

post_session_1$a_collapsed <- array(NA, dim = c(2, 75, 100))
post_session_1$a_collapsed[1, , ] <- (post_session_1$a[1, , ] +
                                        post_session_1$a[2, , ]) / 2
post_session_1$a_collapsed[2, , ] <- (post_session_1$a[3, , ] +
                                        post_session_1$a[4, , ]) / 2

post_session_1$bias_collapsed <- array(NA, dim = c(2, 75, 100))
post_session_1$bias_collapsed[1, , ] <- (post_session_1$bias[1, , ] +
                                        post_session_1$bias[2, , ]) / 2
post_session_1$bias_collapsed[2, , ] <- (post_session_1$bias[3, , ] +
                                        post_session_1$bias[4, , ]) / 2

post_session_1$ndt_collapsed <- array(NA, dim = c(2, 75, 100))
post_session_1$ndt_collapsed[1, , ] <- (post_session_1$ndt[1, , ] +
                                        post_session_1$ndt[2, , ]) / 2
post_session_1$ndt_collapsed[2, , ] <- (post_session_1$ndt[3, , ] +
                                        post_session_1$ndt[4, , ]) / 2

post_session_1$ndt_var_collapsed <- array(NA, dim = c(2, 75, 100))
post_session_1$ndt_var_collapsed[1, , ] <- (post_session_1$ndt_var[1, , ] +
                                        post_session_1$ndt_var[2, , ]) / 2
post_session_1$ndt_var_collapsed[2, , ] <- (post_session_1$ndt_var[3, , ] +
                                        post_session_1$ndt_var[4, , ]) / 2

pred_data_fit_bias <- resimulate_bias(df_session_1, post_session_1)
pred_data_fit_drift <- resimulate_drift(df_session_1, post_session_1)

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

plot_pred <- bind_rows(
  mutate(pred_rt_bias, source = "Only bias modulation"),
  mutate(pred_rt_drift, source = "Only drift modulation")
)
plot_obs <- mutate(emp_rt_summaries, source = "Observed")

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
  # geom_line(
  #   linetype = "dotted",
  # ) +
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
  # geom_line(
  #   data = plot_obs,
  #   inherit.aes = FALSE,
  #   aes(
  #     x = quantile,
  #     y = rt,
  #     color = source,
  #     group = source
  #   ),
  #   linetype = "dotted"
  # ) +
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
# quantiles <- seq(0.1, 0.9, 0.1)
# 
# caf_obs <- df_session_1 %>%
#   group_by(stim_type) %>%
#   summarise(
#     data = list({
#       qs <- quantile(rt, probs = quantiles, na.rm = TRUE)
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x], na.rm = TRUE))
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   unnest(data) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Observed"
#   )
# 
# caf_pred_bias <- pred_data_fit_bias %>%
#   group_by(resim_id, stim_type) %>%
#   summarise(
#     data = list({
#       qs <- quantile(rt, probs = quantiles)
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x]))
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   unnest(data) %>%
#   group_by(stim_type, quantile) %>%
#   summarise(
#     rt = mean(rt),
#     mean = mean(cum_acc),
#     q_lower = quantile(cum_acc, 0.025),
#     q_upper = quantile(cum_acc, 0.975),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Only bias modulation"
#   )
# 
# caf_pred_drift <- pred_data_fit_drift %>%
#   group_by(resim_id, stim_type) %>%
#   summarise(
#     data = list({
#       qs <- quantile(rt, probs = quantiles)
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = map_dbl(qs, ~ mean(correct[rt <= .x]))
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   unnest(data) %>%
#   group_by(stim_type, quantile) %>%
#   summarise(
#     rt = mean(rt),
#     mean = mean(cum_acc),
#     q_lower = quantile(cum_acc, 0.025),
#     q_upper = quantile(cum_acc, 0.975),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Only drift modulation"
#   )
# 
# ggplot(
#   caf_obs,
#   aes(
#     x = rt,
#     y = cum_acc,
#     group = stim_type
#   )
# ) +
#   geom_line(
#     linewidth = 0.8,
#     alpha = 0.6,
#     aes(color = "Observed"),
#     linetype = "dotted"
#   ) +
#   geom_point(
#     shape = 4,
#     size = 3,
#     stroke = 1.3,
#     aes(color = "Observed")
#   ) +
#   geom_line(
#     data = caf_pred_bias,
#     aes(
#       x = rt,
#       y = mean,
#       group = stim_type,
#       color = "Only bias modulation"
#     ),
#     inherit.aes = FALSE,
#     linewidth = 1.0,
#     linetype = "dotted"
#   ) +
#   geom_point(
#     data = caf_pred_bias,
#     aes(
#       x = rt,
#       y = mean,
#       color = "Only bias modulation"
#     ),
#     inherit.aes = FALSE,
#     shape = 4,
#     size = 3,
#     stroke = 1.3
#   ) +
#   geom_line(
#     data = caf_pred_drift,
#     aes(
#       x = rt,
#       y = mean,
#       group = stim_type,
#       color = "Only drift modulation"
#     ),
#     inherit.aes = FALSE,
#     linewidth = 1.0,
#     linetype = "dotted"
#   ) +
#   geom_point(
#     data = caf_pred_drift,
#     aes(
#       x = rt,
#       y = mean,
#       color = "Only drift modulation"
#     ),
#     inherit.aes = FALSE,
#     shape = 4,
#     size = 3,
#     stroke = 1.3
#   ) +
#   facet_grid(
#     stim_type ~ .,
#     labeller = labeller(
#       stim_type = c(
#         "0" = "New\nstatements",
#         "1" = "Repeated\nstatements"
#       )
#     )
#   ) +
#   scale_color_manual(
#     name = "",
#     breaks = c(
#       "Observed",
#       "Only bias modulation",
#       "Only drift modulation"
#     ),
#     values = c(
#       "Only bias modulation" = COLOR_PALETTE[1],
#       "Only drift modulation" = COLOR_PALETTE[2],
#       "Observed" = COLOR_PALETTE[3]
#     )
#   ) +
#   labs(
#     x = "Response time (s)",
#     y = "Cumulative accuracy"
#   ) +
#   ggthemes::theme_tufte(base_size = FONT_SIZE_2 - SMALLER_FONT) +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 12)),
#     axis.title.y = element_text(margin = margin(r = 12)),
#     axis.line = element_line(linewidth = 0.5, color = "#969696"),
#     axis.ticks = element_line(color = "#969696"),
#     axis.text.x = element_text(size = FONT_SIZE_3 - SMALLER_FONT, vjust = 0.5),
#     axis.text.y = element_text(size = FONT_SIZE_3 - SMALLER_FONT),
#     strip.text.x = element_text(size = FONT_SIZE_2 - SMALLER_FONT, angle = 0),
#     strip.text.y = element_text(size = FONT_SIZE_2 - SMALLER_FONT, angle = 0),
#     panel.grid.major = element_line(color = scales::alpha("gray70", 0.3)),
#     panel.grid.minor = element_line(color = scales::alpha("gray70", 0.15)),
#     panel.background = element_blank(),
#     panel.spacing = unit(1.2, "lines"),
#     legend.position = "bottom",
#     legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
#     legend.spacing.y = unit(0.2, "cm")
#   )
# # 
# # ggsave(
# #   '../plots/07_cum_acc_rt_quantiles_restricted_models.jpeg',
# #   device = 'jpeg', dpi = 300,
# #   width = 10, height = 6
# # )
# 
# quantiles <- seq(0.1, 0.9, 0.1)
# N_total <- nrow(df_session_1)
# 
# caf_obs <- df_session_1 %>%
#   group_by(stim_type) %>%
#   summarise(
#     data = list({
#       d <- pick(rt, correct)   # 👈 modern replacement for cur_data()
#       qs <- quantile(d$rt, probs = quantiles, na.rm = TRUE)
#       
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = purrr::map_dbl(
#           qs,
#           ~ sum(d$correct == 1 & d$rt <= .x, na.rm = TRUE) / N_total
#         )
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   tidyr::unnest(data) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Observed"
#   )
# 
# 
# caf_pred_bias <- pred_data_fit_bias %>%
#   group_by(resim_id, stim_type) %>%
#   summarise(
#     data = list({
#       d <- pick(rt, correct)
#       qs <- quantile(d$rt, probs = quantiles)
#       
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = purrr::map_dbl(
#           qs,
#           ~ sum(d$correct == 1 & d$rt <= .x) / nrow(d)
#         )
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   tidyr::unnest(data) %>%
#   group_by(stim_type, quantile) %>%
#   summarise(
#     rt = mean(rt),
#     mean = mean(cum_acc),
#     q_lower = quantile(cum_acc, 0.025),
#     q_upper = quantile(cum_acc, 0.975),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Only bias modulation"
#   )
# 
# 
# caf_pred_drift <- pred_data_fit_drift %>%
#   group_by(resim_id, stim_type) %>%
#   summarise(
#     data = list({
#       d <- pick(rt, correct)
#       qs <- quantile(d$rt, probs = quantiles)
#       
#       tibble(
#         quantile = quantiles,
#         rt = qs,
#         cum_acc = purrr::map_dbl(
#           qs,
#           ~ sum(d$correct == 1 & d$rt <= .x) / nrow(d)
#         )
#       )
#     }),
#     .groups = "drop"
#   ) %>%
#   tidyr::unnest(data) %>%
#   group_by(stim_type, quantile) %>%
#   summarise(
#     rt = mean(rt),
#     mean = mean(cum_acc),
#     q_lower = quantile(cum_acc, 0.025),
#     q_upper = quantile(cum_acc, 0.975),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     quantile = as_factor(quantile),
#     source = "Only drift modulation"
#   )
# 
# 
# ggplot(
#   caf_obs,
#   aes(
#     x = rt,
#     y = cum_acc,
#     group = stim_type
#   )
# ) +
# geom_line(
#   linewidth = 0.8,
#   alpha = 0.6,
#   linetype = "dotted",
#   aes(color = "Observed")
# ) +
#   geom_point(
#     shape = 4,
#     size = 3,
#     stroke = 1.3,
#     aes(color = "Observed")
#   ) +
# geom_line(
#   data = caf_pred_bias,
#   aes(
#     x = rt,
#     y = mean,
#     group = stim_type,
#     color = "Only bias modulation"
#   ),
#   inherit.aes = FALSE,
#   linewidth = 1.0,
#   linetype = "dotted"
# ) +
#   geom_point(
#     data = caf_pred_bias,
#     aes(
#       x = rt,
#       y = mean,
#       group = stim_type,
#       color = "Only bias modulation"
#     ),
#     inherit.aes = FALSE,
#     shape = 4,
#     size = 3,
#     stroke = 1.3
#   ) +
# geom_line(
#   data = caf_pred_drift,
#   aes(
#     x = rt,
#     y = mean,
#     group = stim_type,
#     color = "Only drift modulation"
#   ),
#   inherit.aes = FALSE,
#   linewidth = 1.0,
#   linetype = "dotted"
# ) +
# geom_point(
#   data = caf_pred_drift,
#   aes(
#     x = rt,
#     y = mean,
#     group = stim_type,
#     color = "Only drift modulation"
#   ),
#   inherit.aes = FALSE,
#   shape = 4,
#   size = 3,
#   stroke = 1.3
# ) +
# facet_grid(
#   stim_type ~ .,
#   labeller = labeller(
#     stim_type = c(
#       "0" = "New\nstatements",
#       "1" = "Repeated\nstatements"
#     )
#   )
# ) +
# scale_color_manual(
#   name = "",
#   breaks = c(
#     "Observed",
#     "Only bias modulation",
#     "Only drift modulation"
#   ),
#   values = c(
#     "Observed" = COLOR_PALETTE[3],
#     "Only bias modulation" = COLOR_PALETTE[1],
#     "Only drift modulation" = COLOR_PALETTE[2]
#   )
# ) +
# labs(
#   x = "Response time (s)",
#   y = "Cumulative accuracy"
# ) +
# ggthemes::theme_tufte(base_size = FONT_SIZE_2 - SMALLER_FONT) +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 12)),
#     axis.title.y = element_text(margin = margin(r = 12)),
#     axis.line = element_line(linewidth = 0.5, color = "#969696"),
#     axis.ticks = element_line(color = "#969696"),
#     axis.text.x = element_text(
#       size = FONT_SIZE_3 - SMALLER_FONT,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = FONT_SIZE_3 - SMALLER_FONT),
#     strip.text.x = element_text(
#       size = FONT_SIZE_2 - SMALLER_FONT,
#       angle = 0
#     ),
#     strip.text.y = element_text(
#       size = FONT_SIZE_2 - SMALLER_FONT,
#       angle = 0
#     ),
#     panel.grid.major = element_line(
#       color = scales::alpha("gray70", 0.3)
#     ),
#     panel.grid.minor = element_line(
#       color = scales::alpha("gray70", 0.15)
#     ),
#     panel.background = element_blank(),
#     panel.spacing = unit(1.2, "lines"),
#     legend.position = "bottom",
#     legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
#     legend.spacing.y = unit(0.2, "cm")
#   )








vector[4] transf_mu_v;
vector[4] transf_mu_a;
vector[4] transf_mu_ndt;
vector[4] transf_mu_bias;
vector[4] transf_mu_ndt_var;

params <- c(
  "transf_mu_v", "transf_mu_a", "transf_mu_ndt",
  "transf_mu_bias", "transf_mu_ndt_var"
)

idx <- sample(1:NUM_POST_SAMPLES, NUM_RESIMS, replace = FALSE)

s_val <- extract_all_params(fit_session_1, df_session_1, idx)$s

post_draws <- fit_session_1$draws(variables = params,format = "df") %>% 
  select(-c(.iteration, .chain)) %>% 
  filter(.draw %in% idx)

bias_not_collapsed <- array(NA, dim = c(4, NUM_RESIMS))
bias_not_collapsed[1, ] <- post_draws$`transf_mu_bias[1]`
bias_not_collapsed[2, ] <- post_draws$`transf_mu_bias[2]`
bias_not_collapsed[3, ] <- post_draws$`transf_mu_bias[3]`
bias_not_collapsed[4, ] <- post_draws$`transf_mu_bias[4]`

v_not_collapsed <- array(NA, dim = c(4, NUM_RESIMS))
v_not_collapsed[1, ] <- post_draws$`transf_mu_v[1]`
v_not_collapsed[2, ] <- post_draws$`transf_mu_v[2]`
v_not_collapsed[3, ] <- post_draws$`transf_mu_v[3]`
v_not_collapsed[4, ] <- post_draws$`transf_mu_v[4]`

v_collapsed <- array(NA, dim = c(2, NUM_RESIMS))
v_collapsed[1, ] <- (post_draws$`transf_mu_v[1]` + post_draws$`transf_mu_v[2]`) / 2
v_collapsed[2, ] <- (post_draws$`transf_mu_v[3]` + post_draws$`transf_mu_v[4]`) / 2

a_collapsed <- array(NA, dim = c(2, NUM_RESIMS))
a_collapsed[1, ] <- (post_draws$`transf_mu_a[1]` + post_draws$`transf_mu_a[2]`) / 2
a_collapsed[2, ] <- (post_draws$`transf_mu_a[3]` + post_draws$`transf_mu_a[4]`) / 2

ndt_collapsed <- array(NA, dim = c(2, NUM_RESIMS))
ndt_collapsed[1, ] <- (post_draws$`transf_mu_ndt[1]` + post_draws$`transf_mu_ndt[2]`) / 2
ndt_collapsed[2, ] <- (post_draws$`transf_mu_ndt[3]` + post_draws$`transf_mu_ndt[4]`) / 2

bias_collapsed <- array(NA, dim = c(2, NUM_RESIMS))
bias_collapsed[1, ] <- (post_draws$`transf_mu_bias[1]` + post_draws$`transf_mu_bias[2]`) / 2
bias_collapsed[2, ] <- (post_draws$`transf_mu_bias[3]` + post_draws$`transf_mu_bias[4]`) / 2

ndt_var_collapsed <- array(NA, dim = c(2, NUM_RESIMS))
ndt_var_collapsed[1, ] <- (post_draws$`transf_mu_ndt_var[1]` + post_draws$`transf_mu_ndt_var[2]`) / 2
ndt_var_collapsed[2, ] <- (post_draws$`transf_mu_ndt_var[3]` + post_draws$`transf_mu_ndt_var[4]`) / 2

#################
# BIAS
#################
n_trials <- nrow(df_session_1)
bias_df <- tibble()
for (i in 1:NUM_RESIMS) {
  for (t in 1:n_trials) {
    cond <- df_session_1$condition[t]
    truth <- df_session_1$factual_truth[t]
    v <- v_collapsed[truth, i]
    a <- a_collapsed[truth, i]
    bias <- bias_not_collapsed[cond, i]
    ndt <- ndt_collapsed[truth, i]
    ndt_var <- ndt_var_collapsed[truth, i]
    s <- s_val[t, i]
    ndt_current <- ndt + s * ndt_var
    x <- sample_ddm(v = v, a = a, ndt = ndt_current, bias = bias)
    bias_df <- bind_rows(
      bias_df,
      tibble(
        resp = x[1],
        rt = x[2],
        condition = cond,
        resim_id = i,
        stim_type = df_session_1$stim_type[t],
        factual_truth = df_session_1$factual_truth[t],
        correct = df_session_1$correct[t]
      )
    )
  }
  cat("Resimulation", i, "of", NUM_RESIMS, "done\n")
}

#################
# DRIFT
#################
n_trials <- nrow(df_session_1)
drift_df <- tibble()
for (i in 1:NUM_RESIMS) {
  for (t in 1:n_trials) {
    cond <- df_session_1$condition[t]
    truth <- df_session_1$factual_truth[t]
    v <- v_not_collapsed[cond, i]
    a <- a_collapsed[truth, i]
    bias <- bias_collapsed[truth, i]
    ndt <- ndt_collapsed[truth, i]
    ndt_var <- ndt_var_collapsed[truth, i]
    s <- s_val[t, i]
    ndt_current <- ndt + s * ndt_var
    x <- sample_ddm(v = v, a = a, ndt = ndt_current, bias = bias)
    drift_df <- bind_rows(
      drift_df,
      tibble(
        resp = x[1],
        rt = x[2],
        condition = cond,
        resim_id = i,
        stim_type = df_session_1$stim_type[t],
        factual_truth = df_session_1$factual_truth[t],
        correct = df_session_1$correct[t]
      )
    )
  }
  cat("Resimulation", i, "of", NUM_RESIMS, "done\n")
}


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
  bias_df,
  "Only bias modulation"
)

pred_rt_drift <- summarise_pred_quantiles(
  drift_df,
  "Only drift modulation"
)

plot_pred <- bind_rows(
  mutate(pred_rt_bias, source = "Only bias modulation"),
  mutate(pred_rt_drift, source = "Only drift modulation")
)
plot_obs <- mutate(emp_rt_summaries, source = "Observed")



pd <- position_dodge(width = 0.8)

grp_lvl_plot <- ggplot(
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
  # geom_line(
  #   linetype = "dotted",
  # ) +
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
  # geom_line(
  #   data = plot_obs,
  #   inherit.aes = FALSE,
  #   aes(
  #     x = quantile,
  #     y = rt,
  #     color = source,
  #     group = source
  #   ),
  #   linetype = "dotted"
  # ) +
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
  '../plots/07_rt_quantiles_restricted_models_group_level.jpeg',
  grp_lvl_plot,
  device = 'jpeg', dpi = 300,
  width = 10, height = 6
)


