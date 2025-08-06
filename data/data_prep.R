library(tidyverse)
library(magrittr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

df <- read.csv2('../data/raw_data.csv')

df %<>%
  select(
    id, trial, session, statement_type,
    factual_truth, responses, dur_trial
  ) %>% 
  rename(
    stim_type = statement_type,
    resp = responses,
    rt = dur_trial
  ) %>% 
  mutate(
    stim_type = as.factor(stim_type),
    factual_truth = ifelse(factual_truth == 'true', 1, 0)
  ) %>% 
  mutate(correct = ifelse(factual_truth == resp, 1, 0)) %>% 
  filter(rt > 0) %>% 
  group_by(id) %>% 
  filter(
    log(rt) > boxplot(log(rt))$stats[1],
    log(rt) < boxplot(log(rt))$stats[5]
  ) %>% 
  ungroup() %>% 
  mutate(
    id = dense_rank(id),
    condition = case_when(
      factual_truth == 1 & stim_type == 0 ~ 1, # true statement, new
      factual_truth == 1 & stim_type == 1 ~ 2, # true statement, repeated
      factual_truth == 0 & stim_type == 0 ~ 3, # false statement, new
      factual_truth == 0 & stim_type == 1 ~ 4  # false statement, repeated
    )
  )

summary <- df %>% 
  group_by(id, session) %>% 
  summarise(
    N = n(),
    min = min(rt),
    max = max(rt),
    mean = mean(rt),
  )

# save session data
df_session_1 <- df %>% 
  filter(session == 1) %>% 
  write_csv(., '../data/data_session_1.csv')

df_session_2 <- df %>% 
  filter(session == 2) %>% 
  write_csv(., '../data/data_session_2.csv')

################################################################################
# EXPERIMENT 2
################################################################################
df <- read_csv2('../data/exp2_raw_data.csv')

df %<>% 
  rename(
    id = VP_Code,
    factual_truth = FactualTruth,
    stim_type = repetition,
    resp = responses,
    rt = responses_rt,
    rt_question = key_resp_question_rt
  ) %>% 
  mutate(
    id = dense_rank(id),
    factual_truth = as.numeric(factual_truth),
    correct = ifelse(resp == factual_truth, 1, 0),
    rt = as.numeric(rt),
    rt_question = as.numeric(rt_question)
  ) %>% 
  select(
    id, stim_type, factual_truth,
    resp, correct, rt, rt_question
  ) %>% 
  arrange(id) %>% 
  group_by(id) %>% 
  mutate(
    rt = ifelse(log(rt) > boxplot(log(rt))$stats[1], rt, NA),
    rt = ifelse(log(rt) < boxplot(log(rt))$stats[5], rt, NA),
    rt_question = ifelse(log(rt_question) > boxplot(log(rt_question))$stats[1], rt_question, NA),
    rt_question = ifelse(log(rt_question) < boxplot(log(rt_question))$stats[5], rt_question, NA),
    rt_question = ifelse(rt_question < 1, NA, rt_question)
  ) %>% 
  ungroup() %>% 
  drop_na(rt, rt_question) %>% 
  ungroup() %>% 
  mutate(
    condition = case_when(
      factual_truth == 1 & stim_type == 0 ~ 1, # true statement, new
      factual_truth == 1 & stim_type == 1 ~ 2, # true statement, repeated
      factual_truth == 0 & stim_type == 0 ~ 3, # false statement, new
      factual_truth == 0 & stim_type == 1 ~ 4  # false statement, repeated
    )
  )

write_csv(df, 'data_exp_2.csv')
