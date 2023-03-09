# como_sq_plot
## Read in comorbidity squared models and plot these

library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)

## Labels 
cond_lbl <- read_csv("Outputs/lookup_cond_lbl.csv")
cmpr_lbl <- read_csv("Outputs/lookup_tx_compar_lbl.csv")
Relbl <- function(x) x %>% 
  inner_join(cond_lbl) %>% 
  select(-condition_cb) %>% 
  inner_join(cmpr_lbl) %>% 
  select(-cmpr, -cmpr2) %>% 
  select(`Index condition` = condition_cb2, `Treatment comparison` = cmpr3, everything()) 
oth <- readRDS("Scratch_data/mdls_mvn_como_sq.Rds")
sng <- readRDS("Scratch_data/mdls_mvn_sngl_trials_como_sq.Rds")
oth$drws <- map(oth$res, ~ .x %>% 
                  as_draws_df() %>% 
                  select(arm = b_termarm, como_cnt = `b_termarm:como_cnt`, como_cnt_sq = `b_termarm:Icomo_cntE2`) %>% 
                  as.matrix())
sng <- sng %>% 
  filter(model_type == "cnt_sq") %>% 
  select(model_type:cmpr, trial_n, cfs, bd) 
sng$bd <- map(sng$bd, ~ as.matrix(.x))
sng$drws <- map2(sng$cfs, sng$bd, ~ mvtnorm::rmvnorm(n = 4000, mean = .x$estimate, sigma = .y))

MakeSmry <- function(drws){
  lvls <- tibble(b0 = 0L, b1 = 0:3, b2 = (0:3)^2) %>% 
   as.matrix() 
  cfs <- t(drws)
  out <- t(lvls %*% cfs)
  tibble(
    como_cnt = 0:3,
    mns = colMeans(out),
    sds = apply(out, 2, sd),
    lci = apply(out, 2, quantile, probs = 0.025),
    uci = apply(out, 2, quantile, probs = 0.975))           
}
oth$smry <- map(oth$drws, MakeSmry)
sng$smry <- map(sng$drws, MakeSmry)

smry1 <- oth %>% 
  select(model_type, condition_cb, cmpr, trial_n, smry) %>% 
  unnest(smry)
smry2 <- sng %>% 
  select(model_type, condition_cb, cmpr, trial_n, smry) %>% 
  unnest(smry)
smry <- bind_rows(smry1, smry2)
rm(smry1, smry2)

smry <- smry %>% 
  mutate(forfac = paste0(condition_cb, 
                         "\n", 
                         cmpr,
                         "\n",
                         trial_n, " trials"))
smry <- smry %>% 
  arrange(condition_cb, cmpr) %>% 
  group_by(condition_cb) %>% 
  mutate(cmpr_n = cumsum(!duplicated(cmpr)) %>% as.factor()) %>% 
  ungroup()

smry <- Relbl(smry)

smry_nst <- smry %>% 
  mutate(forfac = `Index condition`) %>% 
  group_by(`Index condition`) %>% 
  nest() %>% 
  ungroup()
smry_nst$cmpr_ns <- map_int(smry_nst$data, ~ sum(!duplicated(.x$`Treatment comparison`)))
smry_nst$plot <- map2(smry_nst$data, smry_nst$cmpr_ns, ~ {
if (.y >= 2) {
  myplot <- ggplot(.x, aes(x = como_cnt, 
                                            y = mns, 
                                            ymin = lci, 
                                            ymax = uci, 
                                            colour = `Treatment comparison`,
                                            fill = `Treatment comparison`)) +
  geom_line() +
  geom_ribbon(alpha = 0.1, colour = NA) +
  coord_cartesian(ylim = c(-2, 2)) +
    facet_wrap(~forfac)} else {
    
  myplot <-  ggplot(.x, aes(x = como_cnt, 
                                             y = mns, 
                                             ymin = lci, 
                                             ymax = uci)) +
  geom_line() +
  geom_ribbon(alpha = 0.1, colour = NA) +
  coord_cartesian(ylim = c(-2, 2)) +
    facet_wrap(~forfac)
}
  myplot
})

## Note slightly different plot sizes
lng1 <- plot_grid(plotlist = smry_nst$plot[smry_nst$cmpr_ns == 1])
lng2 <- plot_grid(plotlist = smry_nst$plot[smry_nst$cmpr_ns == 2])
lng3 <- plot_grid(plotlist = smry_nst$plot[smry_nst$cmpr_ns == 3])
## note only one here - diabetes - with more than 3 trials
lng4 <- plot_grid(plotlist = smry_nst$plot[smry_nst$cmpr_ns > 3])
pdf("Outputs/como_sq_plot1.pdf", width = 10, height = 10)
lng1
dev.off()
pdf("Outputs/como_sq_plot2.pdf", width = 15, height = 10)
lng2
dev.off()
pdf("Outputs/como_sq_plot3.pdf", width = 20, height = 20)
lng3
dev.off()
pdf("Outputs/como_sq_plot4.pdf", width = 10, height = 10)
lng4
dev.off()

