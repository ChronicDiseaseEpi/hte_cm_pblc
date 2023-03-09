# fit models non-cont
library(tidyverse)
library(brms)

dfs <- readRDS("Scratch_data/cfs_non_cont.Rds")
dfs <- map(dfs, ~ .x %>% 
             group_by(model_n) %>% 
             nest())
df_cnt <- map(dfs, ~ .x$data[[1]] %>% distinct(nct_id)) %>% 
  bind_rows(.id = "outcome") %>% 
  mutate(v = 1L) %>% 
  spread(outcome, v)
trials_cnt <- map_int(df_cnt %>% select(-nct_id), sum, na.rm = TRUE)
list2env(dfs, envir = .GlobalEnv)
rm(dfs)

bprior <- c(set_prior("normal(0, 1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "sd"))

## Headache models
mdl <- brm(estimate |se(se) ~ 1 + (1|nct_id), 
           data = hdch$data[hdch$model_n == 0], prior = bprior, chains = 0)
hdch$res <- map(hdch$data, ~ update(mdl, newdata = .x, 
                                    control = list(adapt_delta = 0.999)))
saveRDS(hdch, "Scratch_data/headache_models.Rds")
rm(hdch, mdl)

## osteoporosis models; all three have a placebo/usual care comparator
bpriorfixed <- set_prior("normal(0, 1)", class = "Intercept")
osteop <- osteop %>% 
  unnest(data) %>% 
  ungroup()
osteop <- osteop %>% 
  group_by(model_n, medicine) %>% 
  nest() %>% 
  ungroup()
osteop$n_trials <- map_int(osteop$data, nrow)

mdl <- brm(estimate |se(std.error) ~ 1 , 
           data = osteop$data[[1]], prior = bpriorfixed, chains = 0, backend = "cmdstanr")
osteop$res[osteop$n_trials > 1] <- map(osteop$data[osteop$n_trials > 1], ~ update(mdl, newdata = .x, 
                                        control = list(adapt_delta = 0.999)))
saveRDS(osteop, "Scratch_data/osteop_models.Rds")
rm(osteop, mdl)

## tbe_bleed models
tbe_bleed$data <- map(tbe_bleed$data, ~ .x %>% filter(std.error < 100))
tbe_bleed <- tbe_bleed %>% 
  unnest(data)
tbe_bleed <- tbe_bleed %>% 
  group_by(model_n, condition, cmpr) %>% 
  nest() %>% 
  ungroup()
tbe_bleed$n_trials <- map_int(tbe_bleed$data, nrow)
mdl <- brm(estimate |se(std.error) ~ 1 , 
           data = tbe_bleed$data[[1]], prior = bpriorfixed, chains = 0)
tbe_bleed$res[tbe_bleed$n_trials > 1] <-
  map(tbe_bleed$data[tbe_bleed$n_trials > 1], ~ update(mdl, newdata = .x, control = list(adapt_delta = 0.999)))
saveRDS(tbe_bleed, "Scratch_data/tbe_bleed_models.Rds")
rm(tbe_bleed, mdl)

## tbe thrombolic event models
tbe_dvtpe$data <- map(tbe_dvtpe$data, ~ .x %>% filter(std.error < 100))
tbe_dvtpe <- tbe_dvtpe %>% 
  unnest(data)
tbe_dvtpe <- tbe_dvtpe %>% 
  group_by(model_n, condition, cmpr) %>% 
  nest() %>% 
  ungroup()
tbe_dvtpe$n_trials <- map_int(tbe_dvtpe$data, nrow)
mdl <- brm(estimate |se(std.error) ~ 1, 
           data = tbe_dvtpe$data[[1]], prior = bpriorfixed, chains = 0)
tbe_dvtpe$res[tbe_dvtpe$n_trials > 1] <-
  map(tbe_dvtpe$data[tbe_dvtpe$n_trials > 1], ~ update(mdl, newdata = .x, control = list(adapt_delta = 0.999)))
saveRDS(tbe_dvtpe, "Scratch_data/tbe_dvtpe_models.Rds")
rm(tbe_dvtpe, mdl)
