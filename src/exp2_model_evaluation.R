library(tidyverse)
library(magrittr)
library(rstan)
library(bayesplot)
library(cmdstanr)
library(rempsyc)
library(flextable)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fit <- readRDS("../fits/exp2_full_model_without_var.rds")

param_names <- c(
  "transf_mu_v[1]", "transf_mu_v[2]", "transf_mu_v[3]", "transf_mu_v[4]",
  "transf_mu_a[1]", "transf_mu_a[2]", "transf_mu_a[3]", "transf_mu_a[4]",
  "transf_mu_ndt[1]", "transf_mu_ndt[2]", "transf_mu_ndt[3]", "transf_mu_ndt[4]",
  "transf_mu_bias[1]", "transf_mu_bias[2]", "transf_mu_bias[3]", "transf_mu_bias[4]"
)

write_csv(fit$summary(variables=param_names), "../data/exp2_posterior_summaries.csv")

COLOR_PALETTE <- c('#27374D', '#B70404')
FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18

post_samples <- fit$draws(
  variables = param_names,
  inc_warmup = FALSE,
  format = "draws_matrix"
)

post_samples <- as_tibble(post_samples) %>% 
  mutate(draw = 1:8000) %>% 
  pivot_longer(
    cols = starts_with("transf"),
    names_to = c("parameter", "condition"),
    names_pattern = "(mu_\\w+)\\[(\\d+)\\]"
  )

post_samples <- post_samples %>% 
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

plot <- post_samples %>% 
  ggplot(aes(x = factual_truth, y = value, colour = stim_type, fill = stim_type)) +
  geom_violin(alpha = 0.75) +
  facet_wrap(~parameter, scales = "free") +
  theme_classic() +
  scale_color_manual(values = COLOR_PALETTE) +
  scale_fill_manual(values = COLOR_PALETTE)

ggsave("../plots/exp2_param_estimates.png", plot)

posterior_summaries <- fit$summary(variables = param_names) %>%
  as_tibble() %>% 
  mutate_if(is.numeric, round, digits=2)

posterior_summaries$variable <- c(
  "mu_v_1", "mu_v_2", "mu_v_3", "mu_v_4", 
  "mu_a_1", "mu_a_2", "mu_a_3", "mu_a_4", 
  "mu_ndt_1", "mu_ndt_2", "mu_ndt_3", "mu_ndt_4", 
  "mu_bias_1", "mu_bias_2", "mu_bias_3", "mu_bias_4"
)

write_csv(posterior_summaries, "exp2_posterior_summaries.csv")

apa_table <- nice_table(posterior_summaries, title = "Posterior Summaries Experiment 2")
set_flextable_defaults(fonts_ignore=TRUE)
print(apa_table, preview = "pdf")

## Contrasts
contrast <- tibble()
for (name in unique(post_samples$parameter)){
  # main effect repetition
  effect_repetition <- ((post_samples$value[post_samples$parameter == name & post_samples$condition == 2] +
                           post_samples$value[post_samples$parameter == name & post_samples$condition == 4]) / 2) -
    ((post_samples$value[post_samples$parameter == name & post_samples$condition == 1] +
        post_samples$value[post_samples$parameter == name & post_samples$condition == 3]) / 2)
  # main effect factual truth
  effect_truth <- ((post_samples$value[post_samples$parameter == name & post_samples$condition == 1] +
                      post_samples$value[post_samples$parameter == name & post_samples$condition == 2]) / 2) -
    ((post_samples$value[post_samples$parameter == name & post_samples$condition == 2] +
        post_samples$value[post_samples$parameter == name & post_samples$condition == 4]) / 2)
  # interaction
  effect_interaction <- ((post_samples$value[post_samples$parameter == name & post_samples$condition == 1] +
                            post_samples$value[post_samples$parameter == name & post_samples$condition == 4]) / 2) -
    ((post_samples$value[post_samples$parameter == name & post_samples$condition == 2] +
        post_samples$value[post_samples$parameter == name & post_samples$condition == 3]) / 2)
  tmp_df <- tibble(
    "parameter" = rep(name, 8000*3),
    "effect" = rep(c("Repetition status", "Factual truth", "Interaction"), each = 8000),
    "value" = c(effect_repetition, effect_truth, effect_interaction)
  )
  contrast <- rbind(contrast, tmp_df)
}


contrast <- contrast %>%
  mutate(parameter = factor(
    parameter,
    levels = c("mu_v", "mu_a", "mu_ndt", "mu_bias"),  # replace with your actual parameters
    labels = c(expression(mu[v]), expression(mu[a]), expression(mu[ndt]), expression(mu[bias])))
  )

custom_labeller <- labeller(
  parameter = label_parsed,  # Parse only `parameter`
  effect = label_value       # Keep `effect` as strings
)

contrast_plot <- contrast %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = '#27374D', alpha=0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(parameter ~ effect, labeller = custom_labeller) +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.2, 0.2)) + 
  ggtitle("Contrasts Experiment 2") +
  ggthemes::theme_tufte() + 
  ylab("Density")+
  xlab("Value")+
  theme(axis.line = element_line(linewidth = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = FONT_SIZE_3,
                                   angle = 45,
                                   vjust = 0.5),
        axis.text.y = element_text(size = FONT_SIZE_3),
        strip.text.x = element_text(size = FONT_SIZE_2),
        strip.text.y = element_text(size = FONT_SIZE_2, angle = 0),
        text=element_text(size = FONT_SIZE_2),
        plot.title = element_text(size = FONT_SIZE_1,
                                  hjust = 0.5,
                                  face = 'bold'),
        panel.grid = element_line(color = "#969696",
                                  linewidth = 0.2,
                                  linetype = 1),
        legend.spacing.y = unit(0.25, 'cm')) + 
  ggtitle('Group Level Parameter Contrasts')

ggsave("../plots/exp2_contrast_plot.pdf", width=15, contrast_plot)

