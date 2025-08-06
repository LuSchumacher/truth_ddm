library(tidyverse)
library(magrittr)
library(LaplacesDemon)
library(cmdstanr)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ddm_simulator.R")

softmax <- function(x) {
  exp_x <- exp(x)
  return(exp_x / sum(exp_x))
}

softplus <- function(x) {
  y <- ifelse(x > 20, x + log1p(exp(-x)), log1p(exp(x)))
  return(y)
}

softplus_r <- function(x) {
  if (x > 20) {
    return(x + log1p(exp(-x)))
  } else {
    return(log1p(exp(x)))
  }
}

df <- read_csv('../data/data_session_1.csv')

NUM_SIM <- 50
NUM_SUBS <- 75

PARAM_NAMES <- c(
  "transf_mu_v[1]", "transf_mu_v[2]", "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]", "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]", "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]", "transf_mu_bias[3]", "transf_mu_bias[4]"
)

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18

#------------------------------------------------------------------------------#
# Parameter sampling
#------------------------------------------------------------------------------#
# mu_v <- rnorm(4 * NUM_SIM, 2, 1)
# mu_a <- rnorm(4 * NUM_SIM, 3, 0.5)
# mu_bias <- rnorm(4 * NUM_SIM, 0, 0.2)
# mu_ndt <- rnorm(4 * NUM_SIM, 1, 1)
# 
# sigma_v <- as.vector(apply(matrix(rnorm(4 * NUM_SIM, 0, 1), nrow = 4), 2, softmax))
# sigma_a <- as.vector(apply(matrix(rnorm(4 * NUM_SIM, 0, 1), nrow = 4), 2, softmax))
# sigma_bias <- as.vector(apply(matrix(rnorm(4 * NUM_SIM, 0, 1), nrow = 4), 2, softmax))
# sigma_ndt <- as.vector(apply(matrix(rnorm(4 * NUM_SIM, 0, 1), nrow = 4), 2, softmax))
# 
# sim_id <- rep(1:NUM_SIM, each = 4)
# condition <- rep(1:4, times = NUM_SIM)
# group_params <- tibble(
#   "sim_id" = sim_id,
#   "condition" = condition,
#   "mu_v" = mu_v,
#   "mu_a" = mu_a,
#   "mu_bias" = mu_bias,
#   "mu_ndt" = mu_ndt,
#   "sigma_v" = sigma_v,
#   "sigma_a" = sigma_a,
#   "sigma_bias" = sigma_bias,
#   "sigma_ndt" = sigma_ndt
# )
# 
# group_params$tranf_a <- softplus(group_params$mu_a)
# group_params$tranf_ndt <- softplus(group_params$mu_ndt)
# group_params$tranf_bias <- invlogit(group_params$mu_bias)
# 
# individual_params <- tibble()
# for (i in 1:NUM_SIM) {
#   
#   v <- matrix(0, nrow = 4, ncol = NUM_SUBS)
#   a <- matrix(0, nrow = 4, ncol = NUM_SUBS)
#   bias <- matrix(0, nrow = 4, ncol = NUM_SUBS)
#   ndt <- matrix(0, nrow = 4, ncol = NUM_SUBS)
#   
#   for (j in 1:4) {
#     v[j, ] <- rnorm(
#       NUM_SUBS,
#       group_params$mu_v[group_params$sim_id == i & group_params$condition == j],
#       group_params$sigma_v[group_params$sim_id == i & group_params$condition == j]
#     )
#     a[j, ] <- softplus(rnorm(
#       NUM_SUBS,
#       group_params$mu_a[group_params$sim_id == i & group_params$condition == j],
#       group_params$sigma_a[group_params$sim_id == i & group_params$condition == j]
#     ))
#     bias[j, ] <- invlogit(rnorm(
#       NUM_SUBS,
#       group_params$mu_bias[group_params$sim_id == i & group_params$condition == j],
#       group_params$sigma_bias[group_params$sim_id == i & group_params$condition == j]
#     ))
#     ndt[j, ] <- softplus(rnorm(
#       NUM_SUBS,
#       group_params$mu_ndt[group_params$sim_id == i & group_params$condition == j],
#       group_params$sigma_ndt[group_params$sim_id == i & group_params$condition == j]
#     ))
#   }
#   
#   sim_id <- rep(i, 4*NUM_SUBS)
#   condition <- rep(1:4, times = NUM_SUBS)
#   id <- rep(1:NUM_SUBS, each = 4)
#   tmp_params <- tibble(
#     "sim_id" = sim_id,
#     "id" = id,
#     "condition" = condition,
#     "v" = as.vector(v),
#     "a" = as.vector(a),
#     "bias" = as.vector(bias),
#     "ndt" = as.vector(ndt)
#   )
#   individual_params <- individual_params %>% 
#     bind_rows(tmp_params)
# }
# 
# write_csv(group_params, "../data/param_recovery/group_param_samples_recovery.csv")
# write_csv(individual_params, "../data/param_recovery/individual_params_samples_recovery.csv")
# 
# #------------------------------------------------------------------------------#
# # Data generation
# #------------------------------------------------------------------------------#
# sim_data <- tibble()
# for (i in 1:NUM_SIM) {
#   tmp_sim_data <- tibble()
#   for (j in 1:nrow(df)) {
#     v <- individual_params$v[
#       individual_params$sim_id == i &
#         individual_params$condition == df$condition[j] &
#         individual_params$id == df$id[j]
#     ]
#     a <- individual_params$a[
#       individual_params$sim_id == i &
#         individual_params$condition == df$condition[j] &
#         individual_params$id == df$id[j]
#     ]
#     ndt <- individual_params$ndt[
#       individual_params$sim_id == i &
#         individual_params$condition == df$condition[j] &
#         individual_params$id == df$id[j]
#     ]
#     bias <- individual_params$bias[
#       individual_params$sim_id == i &
#         individual_params$condition == df$condition[j] &
#         individual_params$id == df$id[j]
#     ]
#     
#     x <- sample_ddm(v = v, a = a, ndt = ndt, bias = bias)
#     names(x) <- c("resp", "rt")
#     tmp_sim_data <- tmp_sim_data %>% 
#       bind_rows(x)
#   }
#   tmp_sim_data$id <- df$id
#   tmp_sim_data$condition <- df$condition
#   tmp_sim_data$sim_id <- i
#   tmp_sim_data$stim_type <- df$stim_type
#   tmp_sim_data$factual_truth <- df$factual_truth
#   sim_data <- sim_data %>% 
#     bind_rows(tmp_sim_data)
#   print(paste("Simulation Nr.", i, "finished"))
# }
# 
# write_csv(sim_data, "../data/param_recovery/sim_data_recovery.csv")

#------------------------------------------------------------------------------#
# Model fitting
#------------------------------------------------------------------------------#
# init_fun = function(chains=4){
#   L = list()
#   for (i in 1:chains) {
#     L[[i]] = list(
#       mu_v       = 1.0 + runif(4, -0.5, 0.5),
#       sigma_v    = 0.4 + runif(4, -0.05, 0.05),
#       mu_a       = 2.0 + runif(4, -0.2, 0.2),
#       sigma_a    = 0.4 + runif(4, -0.05, 0.05),
#       mu_bias    = 0.0 + runif(4, -0.05, 0.05),
#       sigma_bias = 0.1 + runif(4, -0.05, 0.05),
#       mu_ndt     = 1 + runif(4, -0.05, 0.05),
#       sigma_ndt  = 0.1 + runif(4, -0.01, 0.01),
#       ndt_s      = 0.05 + runif(1, -0.01, 0.01),
#       z_v        = matrix(runif(4*NUM_SUBS, -1, 1), nrow=4),
#       z_a        = matrix(runif(4*NUM_SUBS, -1, 1), nrow=4),
#       z_ndt      = matrix(runif(4*NUM_SUBS, -1, 1), nrow=4),
#       z_bias     = matrix(runif(4*NUM_SUBS, -1, 1), nrow=4),
#       trel       = matrix(runif(4*NUM_SUBS, 0.01, 0.99), nrow=4)
#     )
#   }
#   return(L)
# }
# 
# m_full <- cmdstan_model(
#   '../model/full_model.stan',
#   cpp_options = list(stan_threads = T)
# )
# 
# for (i in 1:NUM_SIM) {
#   fitting_data <- sim_data[sim_data$sim_id == i, ]
#   stan_data = list(
#     `T`           = nrow(fitting_data),
#     N             = NUM_SUBS,
#     subject_id    = fitting_data$id,
#     resp          = fitting_data$resp,
#     factual_truth = fitting_data$factual_truth,
#     condition     = fitting_data$condition,
#     rt            = fitting_data$rt,
#     minRT         = tapply(
#       fitting_data$rt, list(fitting_data$condition, fitting_data$id), min
#     )
#   )
#   
#   fit_sim_data <- m_full$sample(
#     data = stan_data,
#     init = init_fun(),
#     max_treedepth = 10,
#     adapt_delta = 0.9,
#     refresh = 25,
#     iter_sampling = 1000,
#     iter_warmup = 1000,
#     chains = 4,
#     parallel_chains = 4,
#     threads_per_chain = 2,
#     save_warmup = FALSE
#   )
#   
#   fit_sim_data$save_object(paste0("../fits/fit_sim_data_", i))
#   
#   group_param_estimates <- fit_sim_data$summary(variables = PARAM_NAMES)
#   group_param_estimates <- group_param_estimates %>%
#     select(name = 1, median = 3) %>%
#     mutate(
#       parameter = str_extract(name, "(?<=transf_mu_)[a-z]+"),
#       condition = as.integer(str_extract(name, "\\d+")),
#       sim_id = i
#     ) %>% 
#     select(sim_id, condition, parameter, median) %>% 
#     arrange(sim_id, parameter, condition)
#   
#   write_csv(group_param_estimates, paste0("../data/param_recovery/group_param_estimates_recovery_", i, ".csv"))
#   
#   print(paste("Fitting Nr.", i, "finished"))
# }


#------------------------------------------------------------------------------#
# Evaluation
#------------------------------------------------------------------------------#
true_params <- read_csv("../data/param_recovery/group_param_samples_recovery.csv")
true_params <- true_params %>% 
  select(
    sim_id, condition, mu_v, 
    tranf_a, tranf_ndt, tranf_bias
  ) %>% 
  pivot_longer(
    cols = -c(sim_id, condition),
    names_to = "parameter",
    values_to = "value_true"
  ) %>% 
  mutate(parameter = case_when(
    parameter == "mu_v" ~ "v",
    parameter == "tranf_a" ~ "a",
    parameter == "tranf_ndt" ~ "ndt",
    parameter == "tranf_bias" ~ "bias",
  )) %>% 
  arrange(sim_id, parameter, condition)

true_params$type <- "true"

path <- "../data/param_recovery/"
files <- list.files(path, pattern = "group_param_estimates_recovery")

pred_params <- read_csv(paste0(path, files)) %>% 
  arrange(sim_id, parameter, condition)
  
true_params$value_pred <- pred_params$median
true_params <- true_params %>% 
  mutate(
    parameter = factor(parameter, levels = c("v", "a", "ndt", "bias"))
  )

r2_scores <- true_params %>%
  group_by(parameter) %>%
  summarise(
    r2 = cor(value_true, value_pred)^2
  )
true_params_with_r2 <- true_params %>%
  left_join(r2_scores, by = "parameter")
plots <- true_params_with_r2 %>%
  group_split(parameter) %>%
  imap(~{
    df <- .x
    lims <- range(c(df$value_true, df$value_pred))
    r2_label <- sprintf("R² = %.2f", unique(df$r2))
    
    p <- ggplot(df, aes(x = value_true, y = value_pred)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
      annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1,
               label = r2_label, size = 6, fontface = "italic", color = "black") +
      coord_fixed() +
      xlim(lims) + ylim(lims) +
      ggtitle(unique(df$parameter)) +
      labs(
        x = "True",
        y = if (.y == 1) "Estimated" else NULL
      ) +
      ggthemes::theme_tufte() + 
      theme(
        axis.line = element_line(size = .5, color = "#969696"),
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
        panel.spacing = unit(1., "lines"),
        axis.title.y = if (.y == 1) element_text() else element_blank()
      )
    
    p
  })

combined_plot <- wrap_plots(plots, ncol = 4)
ggsave("../plots/06_parameter_recovery_plot.jpeg", 
       plot = combined_plot, 
       width = 12, height = 8, dpi = 300, 
       device = "jpeg")
