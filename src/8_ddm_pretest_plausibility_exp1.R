library(lme4)
library(lmerTest)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read_csv("../data/pretest_plausibility_exp1.csv")

data$Mean_Truth_mc <- data$M - mean(data$M, na.rm = T)

data_session_1 <- data[which(data$session == '1'), ]
data_session_2 <- data[which(data$session == '2'), ]

model_session_1 <- glmer(
  responses ~ Mean_Truth_mc + (1|id),
  family = binomial(link = "logit"),
  data = data_session_1
)

summary(model_session_1)

model_session_2 <- glmer(
  responses ~ Mean_Truth_mc + (1|id),
  family = binomial(link="logit"),
  data = data_session_2
)

summary(model_session_2)
