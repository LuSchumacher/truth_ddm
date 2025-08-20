library(tidyverse)
library(magrittr)
library(bayesplot)
library(tidybayes)
library(progress)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("0_ddm_simulator.R")

df_session_1 <- read_csv('../data/data_session_1.csv')
df_session_2 <- read_csv('../data/data_session_2.csv')
df_exp_2 <- read_csv('../data/data_exp_2.csv')

fit_session_1 <- readRDS("../fits/fit_session_1.rds")
fit_session_2 <- readRDS("../fits/fit_session_2.rds")
fit_exp_2 <- readRDS("../fits/fit_exp_2.rds")

draws_session_1 <- as_draws_matrix(fit_session_1)
draws_session_2 <- as_draws_matrix(fit_session_2)
draws_exp_2 <- as_draws_matrix(fit_exp_2)

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
COLOR_PALETTE <- c('#27374D', '#B70404')

NUM_RESIMS <- 200
NUM_POST_SAMPLES <- 12000
idx <- sample(1:NUM_POST_SAMPLES, NUM_RESIMS, replace = FALSE)

################################################################################
# HELPER FUNCTIONS
################################################################################
resimulate_session <- function(df, params, num_resims = 200, n_cond = 4) {
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
      
      v    <- params$v[cond, id, r]
      a    <- params$a[cond, id, r]
      bias <- params$bias[cond, id, r]
      ndt_base <- params$ndt[cond, id, r]
      ndt_var  <- params$ndt_var[cond, id, r]
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
        factual_truth = df$factual_truth[t]
      )
    }
    
    pred_data[[r_idx]] <- bind_rows(tmp_list)
    cat("Resimulation", r_idx, "of", num_resims, "done\n")
  }
  
  bind_rows(pred_data)
}

extract_param_3d_subset <- function(fit, name, n_id, idx, n_cond = 4) {
  dm <- fit$draws(variables = name, format = "draws_matrix") |> as.matrix()
  dm <- dm[idx, , drop = FALSE]  # only keep sampled draws
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
  dm <- dm[idx, , drop = FALSE]  # subset draws
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

################################################################################
# RE-SIMULATION
################################################################################
# num_id_session_1 <- length(unique(df_session_1$id))
# params_session_1 <- extract_all_params(fit_session_1, df_session_1, idx)
# pred_data_session_1 <- resimulate_session(df_session_1, params_session_1)
# write_csv(pred_data_session_1, "../fits/pred_data_session_1.csv")

# num_id_session_2 <- length(unique(df_session_2$id))
# params_session_2 <- extract_all_params(fit_session_2, df_session_2, idx)
# pred_data_session_2 <- resimulate_session(df_session_2, params_session_2)
# write_csv(pred_data_session_2, "../fits/pred_data_session_2.csv")

# num_id_exp_2 <- length(unique(df_exp_2$id))
# params_exp_2 <- extract_all_params(fit_exp_2, df_exp_2, idx)
# pred_data_exp_2 <- resimulate_session(df_exp_2, params_exp_2)
# write_csv(pred_data_exp_2, "../fits/pred_data_exp_2.csv")

pred_data_session_1 <- read_csv("../fits/pred_data_session_1.csv")
pred_data_session_2 <- read_csv("../fits/pred_data_session_2.csv")
pred_data_exp_2 <- read_csv("../fits/pred_data_exp_2.csv")

################################################################################
# EMPIRICAL SUMMARIES
################################################################################
# RT quantiles
quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)

emp_rt_summaries_s1 <- df_session_1 %>%
  group_by(stim_type) %>%
  summarise(
    rt_quantiles = list(quantile(rt, probs = quantiles)),
    .groups = "drop"
  ) %>%
  unnest_wider(rt_quantiles) %>%
  pivot_longer(
    cols = -stim_type,
    names_to = "quantile",
    values_to = "rt"
  ) %>%
  mutate(session = 1)

emp_rt_summaries_s2 <- df_session_2 %>%
  group_by(stim_type) %>%
  summarise(
    rt_quantiles = list(quantile(rt, probs = quantiles)),
    .groups = "drop"
  ) %>%
  unnest_wider(rt_quantiles) %>%
  pivot_longer(
    cols = -stim_type,
    names_to = "quantile",
    values_to = "rt"
  ) %>%
  mutate(session = 2)

emp_rt_summaries_exp1 <- bind_rows(emp_rt_summaries_s1, emp_rt_summaries_s2) %>% 
  mutate(
    quantile = as.numeric(sub("%", "", quantile)) / 100
  )

emp_rt_summaries_exp2 <- df_exp_2 %>%
  group_by(stim_type) %>%
  summarise(
    rt_quantiles = list(quantile(rt, probs = quantiles)),
    .groups = "drop"
  ) %>%
  unnest_wider(rt_quantiles) %>%
  pivot_longer(
    cols = -stim_type,
    names_to = "quantile",
    values_to = "rt"
  ) %>% 
  mutate(
    quantile = as.numeric(sub("%", "", quantile)) / 100
  )

################################################################################
# RT QUANTILES PLOT: EXPERIMENT 1
################################################################################
pred_rt_summaries_s1 <- pred_data_session_1 %>%
  group_by(resim_id, stim_type) %>%
  summarise(
    q10 = quantile(rt, 0.1),
    q30 = quantile(rt, 0.3),
    q50 = quantile(rt, 0.5),
    q70 = quantile(rt, 0.7),
    q90 = quantile(rt, 0.9),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quantile", values_to = "rt") %>%
  group_by(stim_type, quantile) %>%
  summarise(
    mean = mean(rt),
    q_lower = quantile(rt, 0.025),
    q_upper = quantile(rt, 0.975),
    .groups = "drop"
  ) %>%
  mutate(session = 1)

pred_rt_summaries_s2 <- pred_data_session_2 %>%
  group_by(resim_id, stim_type) %>%
  summarise(
    q10 = quantile(rt, 0.1),
    q30 = quantile(rt, 0.3),
    q50 = quantile(rt, 0.5),
    q70 = quantile(rt, 0.7),
    q90 = quantile(rt, 0.9),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quantile", values_to = "rt") %>%
  group_by(stim_type, quantile) %>%
  summarise(
    mean = mean(rt),
    q_lower = quantile(rt, 0.025),
    q_upper = quantile(rt, 0.975),
    .groups = "drop"
  ) %>%
  mutate(session = 2)

pred_rt_summaries <- bind_rows(pred_rt_summaries_s1, pred_rt_summaries_s2) %>% 
  mutate(
    quantile = recode(
      quantile, q10 = 0.1, q30 = 0.3, q50 = 0.5, q70 = 0.7, q90 = 0.9
    )
  )

SMALLER_FONT <- 6
rt_quantile_plot_exp1 <- ggplot(pred_rt_summaries, aes(x = as.factor(quantile), y = mean)) +
  geom_pointrange(
    aes(ymin = q_lower, ymax = q_upper, color = "Re-simulated"),
    linewidth = 1.2, fatten = 3.5, alpha = 0.9
  ) +
  geom_point(
    data = emp_rt_summaries_exp1,
    aes(x = as.factor(quantile), y = rt, color = "Observed"),
    size = 2, alpha = 0.7
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = COLOR_PALETTE[2], "Re-simulated" = COLOR_PALETTE[1])
  ) +
  facet_grid(
    stim_type ~ session,
    labeller = labeller(
      stim_type = c(
        "0" = "New\nstatements",
        "1" = "Repeated\nstatements"
      ),
      session = c(
        "1" = "Session 1\n(10 minute interval)",
        "2" = "Session 2\n(1 week interval)"
      )
    )
  ) +
  labs(
    x = "Quantile",
    y = "Response time (s)"
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
  '../plots/03_rt_quantiles_exp_1.jpeg',
  rt_quantile_plot_exp1,
  device = 'jpeg', dpi = 300,
  width = 8, height = 5
)

################################################################################
# RT QUANTILES PLOT: EXPERIMENT 2
################################################################################
pred_rt_summaries_exp2 <- pred_data_exp_2 %>%
  group_by(resim_id, stim_type) %>%
  summarise(
    q10 = quantile(rt, 0.1),
    q30 = quantile(rt, 0.3),
    q50 = quantile(rt, 0.5),
    q70 = quantile(rt, 0.7),
    q90 = quantile(rt, 0.9),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("q"), names_to = "quantile", values_to = "rt") %>%
  group_by(stim_type, quantile) %>%
  summarise(
    mean = mean(rt),
    q_lower = quantile(rt, 0.025),
    q_upper = quantile(rt, 0.975),
    .groups = "drop"
  ) %>% 
  mutate(
    quantile = recode(
      quantile, q10 = 0.1, q30 = 0.3, q50 = 0.5, q70 = 0.7, q90 = 0.9
    )
  )

SMALLER_FONT <- 6
rt_quantile_plot_exp2 <- ggplot(pred_rt_summaries_exp2, aes(x = as.factor(quantile), y = mean)) +
  geom_pointrange(
    aes(ymin = q_lower, ymax = q_upper, color = "Re-simulated"),
    linewidth = 1.2, fatten = 3.5, alpha = 0.9
  ) +
  geom_point(
    data = emp_rt_summaries_exp2,
    aes(x = as.factor(quantile), y = rt, color = "Observed"),
    size = 2, alpha = 0.7
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = COLOR_PALETTE[2], "Re-simulated" = COLOR_PALETTE[1])
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
  labs(
    x = "Quantile",
    y = "Response time (s)"
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
  '../plots/03_rt_quantiles_exp_2.jpeg',
  rt_quantile_plot_exp2,
  device = 'jpeg', dpi = 300,
  width = 6, height = 4
)






















# Example usage
n_cond <- 4
n_id   <- length(unique(df_session_1$id))

v_array       <- extract_param(draws1, "v", n_cond, n_id)
a_array       <- extract_param(draws1, "a", n_cond, n_id)
bias_array    <- extract_param(draws1, "bias", n_cond, n_id)
ndt_array     <- extract_param(draws1, "ndt", n_cond, n_id)
ndt_var_array <- extract_param(draws1, "ndt_var", n_cond, n_id)

# pick one fit
fit <- fit_session_1  

# v, a, bias saved as 2D arrays [condition, id]
draws_v    <- as_draws_matrix(fit$draws("v"))
draws_a    <- as_draws_matrix(fit$draws("a"))
draws_bias <- as_draws_matrix(fit$draws("bias"))

# t0 (a.k.a ndt) saved as trial-wise
draws_t0   <- as_draws_matrix(fit$draws("t0"))


# Parse column names like "v[cond,id]" into a tibble with indices
.parse_idx <- function(cols) {
  # supports "name[i,j]" OR "name[i]" (we'll fill missing j with NA)
  m <- str_match(cols, "^[^\\[]+\\[(\\d+)(?:,(\\d+))?\\]$")
  tibble(col = cols,
         i = as.integer(m[,2]),
         j = suppressWarnings(as.integer(m[,3])))
}

# Extract a 3D array [cond × id × draw] for parameters saved as name[cond,id]
extract_param_3d <- function(fit, par) {
  dm <- fit$draws(variables = par, format = "draws_matrix")
  mat <- as.matrix(dm)  # draws × columns
  cols <- colnames(mat)
  idx  <- .parse_idx(cols)
  
  if (any(is.na(idx$j))) {
    stop(sprintf("Parameter '%s' does not look 2D (missing second index). Found columns like: %s",
                 par, paste(head(cols, 3), collapse = ", ")))
  }
  
  n_draw <- nrow(mat)
  n_i <- max(idx$i)
  n_j <- max(idx$j)
  
  arr <- array(NA_real_, dim = c(n_i, n_j, n_draw))
  # fill by columns
  for (k in seq_along(cols)) {
    arr[idx$i[k], idx$j[k], ] <- mat[, k]
  }
  arr
}

# Extract a trial-wise matrix [trial × draw] for parameters like t0[trial] / ndt[trial]
extract_trial_by_draw <- function(fit, par) {
  dm <- fit$draws(variables = par, format = "draws_matrix")
  mat <- as.matrix(dm)  # draws × trials
  idx <- .parse_idx(colnames(mat))
  
  if (any(is.na(idx$i))) {
    stop(sprintf("Could not parse trial indices for '%s'.", par))
  }
  
  n_draw <- nrow(mat)
  n_trial <- max(idx$i)
  out <- matrix(NA_real_, nrow = n_trial, ncol = n_draw)  # trial × draw
  
  # fill by trial index
  for (k in seq_along(idx$i)) {
    out[idx$i[k], ] <- mat[, k]
  }
  out
}

# Fallback extractor if ndt/t0 is also [cond,id] instead of [trial]
extract_param_2d_like_ndt <- function(fit, par) {
  # re-use 3D extractor; caller will index [cond,id,draw]
  extract_param_3d(fit, par)
}

## ---------- build posterior arrays from the fit ----------

build_posterior_struct <- function(fit,
                                   par_v   = "v",
                                   par_a   = "a",
                                   par_bias= "bias",
                                   par_t0_trial = c("t0", "ndt"),   # try these (trial-wise)
                                   par_ndt_condid = c("ndt"),       # fallback (cond×id)
                                   verbose = TRUE) {
  if (verbose) message("Extracting v, a, bias as [cond × id × draw]...")
  v_arr    <- extract_param_3d(fit, par_v)
  a_arr    <- extract_param_3d(fit, par_a)
  bias_arr <- extract_param_3d(fit, par_bias)
  
  # Try trial-wise ndt first
  t0_trial_draw <- NULL
  for (cand in par_t0_trial) {
    ok <- try({
      t0_trial_draw <- extract_trial_by_draw(fit, cand)
      TRUE
    }, silent = TRUE)
    if (isTRUE(ok)) {
      if (verbose) message("Found trial-wise ndt as '", cand, "'. Using [trial × draw].")
      break
    }
  }
  
  ndt_condid <- NULL
  if (is.null(t0_trial_draw)) {
    # fallback to cond×id
    for (cand in par_ndt_condid) {
      ok <- try({
        ndt_condid <- extract_param_3d(fit, cand)
        TRUE
      }, silent = TRUE)
      if (isTRUE(ok)) {
        if (verbose) message("Using ndt from '", cand, "' as [cond × id × draw].")
        break
      }
    }
    if (is.null(ndt_condid)) {
      stop("Could not find ndt/t0 either as trial-wise or as cond×id.")
    }
  }
  
  list(v = v_arr, a = a_arr, bias = bias_arr,
       t0_trial_draw = t0_trial_draw,  # trial × draw (or NULL)
       ndt_condid    = ndt_condid)     # cond × id × draw (or NULL)
}




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
# 
# 
# pred_data_s1 <- tibble()
# for (i in 1:length(idx)) {
#   tmp_post <- df_posterior_s1 %>% 
#     filter(draw == idx[i])
#   tmp_pred_data <- tibble()
#   for (j in 1:nrow(df_s1)) {
#     v <- tmp_post$value[
#       tmp_post$parameter == "v" &
#         tmp_post$participant == df_s1$id[j] &
#         tmp_post$condition == df_s1$condition[j]
#     ]
#     a <- tmp_post$value[
#       tmp_post$parameter == "a" &
#         tmp_post$participant == df_s1$id[j] &
#         tmp_post$condition == df_s1$condition[j]
#     ]
#     ndt <- tmp_post$value[
#       tmp_post$parameter == "ndt" &
#         tmp_post$participant == df_s1$id[j] &
#         tmp_post$condition == df_s1$condition[j]
#     ]
#     bias <- tmp_post$value[
#       tmp_post$parameter == "bias" &
#         tmp_post$participant == df_s1$id[j] &
#         tmp_post$condition == df_s1$condition[j]
#     ]
#     x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
#     names(x) <- c("resp", "rt")
#     tmp_pred_data <- tmp_pred_data %>% 
#       bind_rows(x)
#   }
#   tmp_pred_data$id <- df_s1$id
#   tmp_pred_data$condition <- df_s1$condition
#   tmp_pred_data$resim_id <- i
#   tmp_pred_data$stim_type <- df_s1$stim_type
#   tmp_pred_data$factual_truth <- df_s1$factual_truth
#   pred_data_s1 <- pred_data_s1 %>% 
#     bind_rows(tmp_pred_data)
#   print(paste("Resimulation Nr.", i, "finished"))
# }
# 
# write_csv(pred_data_s1, "../data/pred_data_s1.csv")
pred_data_s1 <- read_csv("../data/pred_data_s1.csv")

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
# fit_session_2 <- readRDS("../fits/fit_session_2")
# df_posterior_s2 <- as_tibble(
#   fit_session_2$draws(
#     inc_warmup = FALSE,
#     format = "draws_matrix"
#   )
# ) %>% 
#   mutate(draw = row_number()) %>% 
#   select(
#     draw,
#     matches("^v\\[[1-4],\\s?\\d+\\]$"),
#     matches("^a\\[[1-4],\\s?\\d+\\]$"),
#     matches("^ndt\\[[1-4],\\s?\\d+\\]$"),
#     matches("^bias\\[[1-4],\\s?\\d+\\]$")
#   ) %>% 
#   pivot_longer(
#     cols = -draw,
#     names_to = "param_full",
#     values_to = "value"
#   ) %>% 
#   mutate(
#     parameter = str_extract(param_full, "^[a-z]+"),  # Capture parameter name (v, a, etc.)
#     condition = as.integer(str_extract(param_full, "(?<=\\[)\\d+")),  # Extract condition (1-4)
#     participant = as.integer(str_extract(param_full, "(?<=,\\s?)\\d+")),  # Extract participant number
#     value = as.numeric(value)
#   ) %>%
#   select(draw, parameter, condition, participant, value)
# 
# pred_data_s2 <- tibble()
# for (i in 1:length(idx)) {
#   tmp_post <- df_posterior_s2 %>% 
#     filter(draw == idx[i])
#   tmp_pred_data <- tibble()
#   for (j in 1:nrow(df_s2)) {
#     v <- tmp_post$value[
#       tmp_post$parameter == "v" &
#         tmp_post$participant == df_s2$id[j] &
#         tmp_post$condition == df_s2$condition[j]
#     ]
#     a <- tmp_post$value[
#       tmp_post$parameter == "a" &
#         tmp_post$participant == df_s2$id[j] &
#         tmp_post$condition == df_s2$condition[j]
#     ]
#     ndt <- tmp_post$value[
#       tmp_post$parameter == "ndt" &
#         tmp_post$participant == df_s2$id[j] &
#         tmp_post$condition == df_s2$condition[j]
#     ]
#     bias <- tmp_post$value[
#       tmp_post$parameter == "bias" &
#         tmp_post$participant == df_s2$id[j] &
#         tmp_post$condition == df_s2$condition[j]
#     ]
#     x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
#     names(x) <- c("resp", "rt")
#     tmp_pred_data <- tmp_pred_data %>% 
#       bind_rows(x)
#   }
#   tmp_pred_data$id <- df_s2$id
#   tmp_pred_data$condition <- df_s2$condition
#   tmp_pred_data$resim_id <- i
#   tmp_pred_data$stim_type <- df_s2$stim_type
#   tmp_pred_data$factual_truth <- df_s2$factual_truth
#   pred_data_s2 <- pred_data_s2 %>% 
#     bind_rows(tmp_pred_data)
#   print(paste("Resimulation Nr.", i, "finished"))
# }

# write_csv(pred_data_s2, "../data/pred_data_s2.csv")

pred_data_s2 <- read_csv("../../data/pred_data_s2.csv")

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
fit_exp_2 <- readRDS("../../fits/fit_exp_2_ndt_var.rds")

df_posterior_exp2 <- fit_exp_2$draws(
  inc_warmup = FALSE,
  format = "draws_matrix"
) %>%
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  select(
    draw,
    matches("^v\\[[1-4],\\s?\\d+\\]$"),
    matches("^a\\[[1-4],\\s?\\d+\\]$"),
    matches("^ndt\\[[1-4],\\s?\\d+\\]$"),
    matches("^bias\\[[1-4],\\s?\\d+\\]$"),
    matches("^ndt_s\\[\\d+\\]$")
  ) %>%
  filter(draw %in% idx) %>% 
  pivot_longer(
    cols = -draw,
    names_to = "param_full",
    values_to = "value"
  ) %>%
  mutate(
    parameter = str_extract(param_full, "^[a-z_]+"),
    indices = str_extract_all(param_full, "\\d+"),
    condition = ifelse(
      parameter == "ndt_s",
      NA_integer_,
      as.integer(sapply(indices, `[`, 1))
    ),
    participant = ifelse(
      parameter == "ndt_s",
      as.integer(sapply(indices, `[`, 1)),
      as.integer(sapply(indices, function(x) if(length(x) > 1) x[2] else NA))
    ),
    value = as.numeric(value)
  ) %>%
  select(draw, parameter, condition, participant, value)


df_post_s_exp_2 <- fit_exp_2$draws(
  inc_warmup = FALSE,
  format = "draws_matrix"
) %>%
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  select(draw, starts_with("s[")) %>%
  filter(draw %in% idx) %>% 
  pivot_longer(
    cols = -draw,
    names_to = "param_full",
    values_to = "value"
  ) %>%
  mutate(
    parameter = "s",
    trial = as.integer(str_extract(param_full, "\\d+")),
    value = as.numeric(value)
  ) %>%
  select(draw, parameter, trial, value)
  
pred_data_exp2 <- tibble()
for (i in 1:length(idx)) {
  tmp_post <- df_posterior_exp2 %>%
    filter(draw == idx[i])
  tmp_s_post <- df_post_s_exp_2 %>% 
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
    ndt_s <- tmp_post$value[
      tmp_post$parameter == "ndt_s" &
        tmp_post$participant == df_exp2$id[j]
    ]
    bias <- tmp_post$value[
      tmp_post$parameter == "bias" &
        tmp_post$participant == df_exp2$id[j] &
        tmp_post$condition == df_exp2$condition[j]
    ]
    t0 <- ndt + tmp_s_post$value[j] * ndt_s
    x <- sample_ddm(v = v, a = a, ndt = t0, bias = bias)
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

write_csv(pred_data_exp2, "../../data/pred_data_exp2_ndt_var.csv")

# pred_data_exp2 <- read_csv("../../data/pred_data_exp2.csv")

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
  '../plots/03_rt_quantiles_plot_exp1.jpeg',
  device = 'jpeg', dpi = 300,
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
  '../../plots/03_rt_quantiles_plot_exp2_ndt_var.jpeg',
  device = 'jpeg', dpi = 300,
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
  '../plots/04_resp_prob_plot_exp1.jpeg',
  resp_prob_plot,
  device = 'jpeg', dpi = 300,
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
        text = element_text(size = FONT_SIZE_2),
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
  '../../plots/04_resp_prob_plot_exp2_ndt_var.jpeg',
  resp_prob_plot_exp2,
  device = 'jpeg', dpi = 300,
  width = 12, height = 6
)
