#fit models multivariate
source("Scripts/01a_prepare_cont_data.R")
library(brms)
count <- 1
## Drop 57 rows where the combined condition/treatment comparisons has only one trial
## from meta-analysis, then save that separately for figure 3 etc
oth$trial_n <- map_int(oth$cfs, ~ sum(!duplicated(.x$nct_id)))
oth_single <- oth %>% 
  filter(trial_n == 1)
oth <- oth %>% 
  filter(trial_n != 1)
oth$mypriors <- vector(mode = "list", length = nrow(oth))
oth$mypriors[oth$formtype == "simple"] <- map(oth$mypriors[oth$formtype == "simple"] , function(x) c(
  prior(constant(1), class = "sigma"),
  set_prior("normal(0, 1)", class = "b")))
oth$mypriors[!oth$formtype == "simple"] <- map(oth$mypriors[!oth$formtype == "simple"] , function(x) c(
  prior(constant(1), class = "sigma"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "sd")))


# Initially ran without comorbidity count squared, so added later and run within QuickRun function
QuickRun <- function(x) {
  x$res <- pmap(list(
    x$myform,
    x$cfs,
    x$bd,
    x$mypriors), function(myform, cfs, bd, mypriors) {
      print(count)
      count <<- count + 1
      brm(myform, 
          data = cfs, 
          data2 = list(covs = bd),
          prior = mypriors,
          control = list(adapt_delta = 0.999),
          cores = 4)  
    })
  x
}
cc2 <- oth %>% 
  filter(model_type == "cnt_sq")
cc2_single <- oth_single %>% 
  filter(model_type == "cnt_sq")
cc2 <- QuickRun(cc2)
saveRDS(cc2, "Scratch_data/mdls_mvn_como_sq.Rds")
saveRDS(cc2_single, "Scratch_data/mdls_mvn_sngl_trials_como_sq.Rds")

oth <- oth %>% 
  filter(!model_type == "cnt_sq")
oth_single <- oth_single %>% 
  filter(!model_type == "cnt_sq")
oth <- QuickRun(oth)
saveRDS(oth, "Scratch_data/mdls_mvn.Rds")
saveRDS(oth_single, "Scratch_data/mdls_mvn_sngl_trials.Rds")


