# fit_models_como_cnt
source("Scripts/01a_prepare_cont_data.R")
library(brms)

## selective FE
cnt$mypriors <- vector(mode = "list", length = nrow(cnt))
cnt$mypriors[cnt$formtype == "simple"] <- map(cnt$mypriors[cnt$formtype == "simple"] , function(x) c(
  set_prior("normal(0, 1)", class = "Intercept")))
cnt$mypriors[!cnt$formtype == "simple"] <- map(cnt$mypriors[!cnt$formtype == "simple"] , function(x) c(
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "sd")))
count <- 1
cnt$res <- pmap(list(cnt$myform,
                     cnt$data,
                     cnt$mypriors), function(myform, cfs, mypriors) {
                       print(count)
                       count <<- count + 1
                       brm(myform,
                           data = cfs,
                           prior = mypriors,
                           control = list(adapt_delta = 0.999),
                           cores = 4,
                           backend = "cmdstanr")
                     })

saveRDS(cnt, "Scratch_data/mdls_cnt.Rds")

## Run as RE for ones currently are FE for appendix table
fe_chs <- cnt %>% 
  filter(formtype == "simple")
re_chs <- cnt %>% 
  filter(!formtype == "simple") %>% 
  select(myform, mypriors) %>% 
  slice(1)
fe_chs <- fe_chs %>% 
  select(-myform, -mypriors) 
fe_chs$myform <- map(1:nrow(fe_chs), ~ re_chs$myform[[1]])
fe_chs$mypriors <- map(1:nrow(fe_chs), ~ re_chs$mypriors[[1]])
count <- 1
fe_chs$res <- pmap(list(fe_chs$myform,
                        fe_chs$data,
                        fe_chs$mypriors), function(myform, cfs, mypriors) {
                       print(count)
                       count <<- count + 1
                       brm(myform,
                           data = cfs,
                           prior = mypriors,
                           control = list(adapt_delta = 0.999),
                           cores = 4,
                           backend = "cmdstanr")
                     })
saveRDS(fe_chs, "Scratch_data/mdls_cnt_re.Rds")
rm(fe_chs)

## Next fit all of these index conditions in a single model.
## Drop outcomes which were not standardised, urate and percent heartburn free days
## Fit a model with ALL in together, note am not asserting is the same,
## just summarising.
## Only do for list of outcomes that were standardised
cnt$res <- NULL
cnt_all <- cnt %>% 
  unnest(data)
cnt_all <- cnt_all %>% 
  filter(!outcome %in% c("Percent heartburn free days", "Urate"))

## Global model ----
mod_glb <- brm(estimate|se(se) ~ 1 + 
                 (1|nct_id) + 
                 (1|cmpr) + 
                 (1|condition_cb), 
               data = cnt_all, 
               control = list(adapt_delta = 0.99),
               cores = 4, backend = "cmdstanr")
summary(mod_glb)
mod_glb_fe <- brm(estimate|se(se) ~ 1 + 
                 (1|cmpr) + 
                 (1|condition_cb), 
               data = cnt_all, 
               control = list(adapt_delta = 0.99),
               cores = 4,  backend = "cmdstanr")
summary(mod_glb_fe)
# similar result for RE and FE for global model
saveRDS(list(re = mod_glb, fe = mod_glb_fe), "Scratch_data/global_models.Rds")
glb <- readRDS("Scratch_data/global_models.Rds")$re
a <- posterior_predict(glb, newdata = data.frame(se = 0, cmpr = "new", condition_hrm = "new", nct_id = "new"), 
                       allow_new_levels = TRUE)
a <- tibble(y = as.vector(a))

## Summarise posteriors using t-distribution ----
smrs_t <- brm(y ~ 1, data = a, family = "student", cores = 4, backend = "cmdstanr")
smrs_t_smry <- summary(smrs_t)
smrs_t_smry <- tibble(m = smrs_t_smry$fixed[1,1],
                      s = smrs_t_smry$spec_pars[1,1], 
                      df = smrs_t_smry$spec_pars[2,1])
a <- bind_rows(a,a,a)
## Sample from t-distribution to compare to raw samples and plot
a$y_t <- rt(12000, df = smrs_t_smry$df)*smrs_t_smry$s + smrs_t_smry$m
a <- a %>% 
  rename(`Posterior` = y,
         `T-distribution` = y_t)
a <- a %>% 
  gather("smpl", "value")
plot1 <- ggplot(a, aes(x = value, colour = smpl)) + 
  geom_density(alpha = 0.2, fill = NA, size = 0.1, ) +
  scale_x_continuous("Comorbidity-count treatment interaction") +
  scale_y_continuous("Density") +
  scale_color_discrete("")
plot1
write_csv(smrs_t_smry, "Outputs/tdist.csv")
