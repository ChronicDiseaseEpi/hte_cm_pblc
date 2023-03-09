#fit models multivariate within VM box
# similar to "Scripts/04a_fit_models_multivariate.R" just re-runs (on Virtual machine)
## Uses brms to write stancode and stan data. Can then run this within the VM using cmdstanr
# with original as well as wider priors
source("Scripts/01a_prepare_cont_data.R")
library(brms)
library(ggforce)
count <- 1

## Select priors based on https://mc-stan.org/docs/stan-users-guide/multivariate-hierarchical-priors.html
oth$trial_n <- map_int(oth$cfs, ~ sum(!duplicated(.x$nct_id)))
oth_single <- oth %>% 
  filter(trial_n == 1)
oth <- oth %>% 
  filter(trial_n != 1)
oth$mypriors <- vector(mode = "list", length = nrow(oth))

for (prior_orig in c(TRUE,FALSE)) {
  if (prior_orig) {
    oth$mypriors[oth$formtype == "simple"] <- map(oth$mypriors[oth$formtype == "simple"] , function(x) c(
      prior(constant(1), class = "sigma"),
      set_prior("normal(0, 1)", class = "b")))
    oth$mypriors[!oth$formtype == "simple"] <- map(oth$mypriors[!oth$formtype == "simple"] , function(x) c(
      prior(constant(1), class = "sigma"),
      set_prior("normal(0, 1)", class = "b"),
      set_prior("normal(0, 1)", class = "sd"),
      set_prior("lkj_corr_cholesky(1)", class = "L")))
  } else {
    oth$mypriors[oth$formtype == "simple"] <- map(oth$mypriors[oth$formtype == "simple"] , function(x) c(
      prior(constant(1), class = "sigma"),
      set_prior("normal(0, 2)", class = "b")))
    oth$mypriors[!oth$formtype == "simple"] <- map(oth$mypriors[!oth$formtype == "simple"] , function(x) c(
      prior(constant(1), class = "sigma"),
      set_prior("normal(0, 2)", class = "b"),
      set_prior("normal(0, 2)", class = "sd"),
      set_prior("lkj_corr_cholesky(2)", class = "L")))
  }
  brms_make <- function(...){
    list(code = make_stancode(...),
         data = make_standata(...))
  }
  oth$res <- pmap(list(
    oth$myform,
    oth$cfs,
    oth$bd,
    oth$mypriors), function(myform, cfs, bd, mypriors) {
                                 print(count)
                                 count <<- count + 1
                                 brms_make(myform, 
                                     data = cfs, 
                                     data2 = list(covs = bd),
                                     prior = mypriors,
                                     control = list(adapt_delta = 0.999),
                                     cores = 4)  
                               })
  if (!prior_orig) {
    saveRDS(oth, "Scratch_data/mdls_mvn_vmbox.Rds")
  } else {
    saveRDS(oth, "Scratch_data/mdls_mvn_vmbox_prior_orig.Rds")
  }
}
## next export both of these to the VM box

## Imported from the vm box
orig <- readRDS("FromVM/mdls_mvn_vmbox_prior_orig.Rds")
new  <- readRDS("FromVM/mdls_mvn_vmbox.Rds")

res <- bind_rows(orig = orig, new = new, .id = "prior_type")
res <- res %>% 
  select(prior_type, model_type, condition_cb, cmpr, models, res_smry) %>% 
  unnest(res_smry)
res <- res %>% 
  arrange(model_type, condition_cb, cmpr, models, variable, prior_type)
res <- res %>% 
  filter( (variable %>% str_remove_all("[0-9]")) %in% c("b[]"))
a <- read_csv("model_type,variable,lbl
agesex,b[1],arm
agesex,b[2],age
agesex,b[3],sex
bin,b[1],v1
bin,b[2],v2
bin,b[3],v3
bin,b[4],v4
bin,b[5],v5
bin,b[6],v6
cnt_sq,b[1],arm
cnt_sq,b[2],como
cnt_sq,b[3],como_sq
intrvl,b[1],bmi
intrvl,b[2],egfr
intrvl,b[3],fib4
intrvl,b[4],hgb
intrvl,b[5],mbp")
res <- res %>% 
  inner_join(a)
write_csv(res %>% select(prior_type:cmpr, mean:q95), "Outputs/sensitivity_prior.csv")

eachplot <- map(c("agesex", "bin", "cnt_sq", "intrvl"), ~ {
  ggplot(res %>% filter(model_type == .x), aes(x = interaction(cmpr, lbl), y = mean, ymin = q5, ymax = q95, colour = prior_type)) +
    geom_point(position = position_dodge(0.1)) +
    geom_linerange(position = position_dodge(0.1)) +
    coord_flip() +
    facet_wrap(~ condition_cb, scales = "free") +
    ggtitle(.x) +
    scale_x_discrete("Treatment comparison/term") +
    scale_y_continuous("Effect estimate")
})
pdf("Outputs/sensitivity_prior.pdf", height = 15, width = 15)
eachplot
dev.off()



