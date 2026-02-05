library(WienR)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

FONT_SIZE_1 <- 22
FONT_SIZE_2 <- 20
FONT_SIZE_3 <- 18
SMALLER_FONT <- 2
COLOR_PALETTE <- c('#27374D', '#B70404')

get_ddm_quantiles <- function(q, a, v, w, t0, sv, sw, st0) {
  pars <- list(a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0)
  plow <- 0; phigh <- 1
  low <- 0; high <- 100;
  pd <- WienerCDF(
    100, 2, pars$a, pars$v, pars$w, pars$t0, pars$sv, pars$sw, pars$st0
  )$value
  while (high - low > 0.0001) {
    m <- (high + low ) / 2
    p <- WienerCDF(
      m, 2, pars$a, pars$v, pars$w, pars$t0, pars$sv, pars$sw, pars$st0
      )$value / pd
    if (p > q) high <- m
    else low <- m
  }
  (high + low) / 2
}

get_ddm_quantiles <- Vectorize(get_ddm_quantiles)

v_effect <- function(q, ef, a, v, w, t0, sv, sw, st0) {
  get_ddm_quantiles(q, a, v, w, t0, sv, sw, st0) -
    get_ddm_quantiles(q, a, v + ef, w, t0, sv, sw, st0)
}

w_effect <- function(q, ef,  a, v, w, t0, sv, sw, st0) {
  get_ddm_quantiles(q, a, v, w, t0, sv, sw, st0) -
    get_ddm_quantiles(q, a, v, w + ef, t0, sv, sw, st0)
}

# effect_on_w and effect_on_v chosen such that the proportion of truth ratings
# increase from 0.53 for new items to 0.75 for repeated items as in the data
effect_on_w <- 0.22
effect_on_v <- 0.4775

# Exp. 1, Session 1, True, new statements
pars <- list(a = 2.03, v = 0, w = 0.53, t0 = 1.74, sv = 0, sw = 0, st0 = 2.22) 

WienerCDF(
  100, 2, pars$a, pars$v,pars$w,pars$t0, pars$sv,pars$sw,pars$st0
)$value
WienerCDF(
  100, 2, pars$a, pars$v + effect_on_v, pars$w, pars$t0, pars$sv, pars$sw, pars$st0
)$value
WienerCDF(
  100, 2, pars$a, pars$v, pars$w + effect_on_w, pars$t0, pars$sv, pars$sw, pars$st0
)$value

qs <- seq(0, 0.99, length.out = 200)
df <- tibble(
  q = qs,
  startpoint_effect = w_effect(
    qs, effect_on_w,
    pars$a, pars$v, pars$w, pars$t0, pars$sv, pars$sw, pars$st0
  ),
  drift_effect = v_effect(
    qs, effect_on_v,
    pars$a, pars$v, pars$w, pars$t0, pars$sv, pars$sw, pars$st0
  )
)

df_long <- df %>% 
  pivot_longer(
    cols = c(startpoint_effect, drift_effect),
    names_to = "effect",
    values_to = "value"
  )

ggplot(df_long, aes(x = q, y = value, linetype = effect)) +
  geom_line(linewidth = 1.0) +
  # scale_color_manual(
  #   values = c(
  #     drift_effect = COLOR_PALETTE[1],
  #     startpoint_effect = COLOR_PALETTE[2]
  #   ),
  #   labels = c(
  #     drift_effect = "Effect on drift rate",
  #     startpoint_effect = "Effect on starting point"
  #   )
  # ) +
  scale_linetype_manual(
    values = c(
      drift_effect = "solid",
      startpoint_effect = "dashed"
    ),
    labels = c(
      drift_effect = "Effect on drift rate",
      startpoint_effect = "Effect on starting point"
    )
  ) +
  labs(
    x = "Quantile",
    y = "Decrease in RT for repeated statements (s)",
    color = NULL,
    linetype = NULL
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
    legend.spacing.y = unit(0.2, "cm"),
    legend.key.width = unit(1.5, "cm")
  )

ggsave(
  '../plots/07_qualitative_effect_rt_quantiles.jpeg',
  device = 'jpeg', dpi = 300,
  width = 8.5, height = 5.5
)
