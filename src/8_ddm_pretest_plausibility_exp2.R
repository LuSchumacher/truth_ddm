library(lme4)
library(lmerTest)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read_csv("../data/pretest_plausibility_exp2.csv")

data$Mean_Truth_mc <- data$M - mean(data$M, na.rm = T)

model_exp_2 <- glmer(
  responses ~ Mean_Truth_mc + (1 | id),
  family = binomial(link = "logit"),
  data = data
)

summary(model_exp_2)
