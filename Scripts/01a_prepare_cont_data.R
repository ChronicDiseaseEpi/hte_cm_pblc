library(tidyverse)
## Bring data together and create formulae and datasets for brms modelling
## also create Table 1

## Read data files ----
rlbl_out <- read_csv("Supporting/relbl_outcomes.csv")
cfs <- read_csv("Data/Model_coefficients.csv.gz")

## Functions ----
# Converts covariance/correlation data from a long format into a
# matrix format needed for model fitting
# Was originally converted into this format prior to export for concision
# Note for CSDR is variance-covariance matrix 
# (checked against standard errors) not a correlation matrix
# whereas YODA extract is a correlation matrix
CnvrtMatrix <- function(a){
  ## a is a dataframe
  ## Duplicate dataframe switching rows and cols
  a <- bind_rows(a, 
                 a %>% 
                   rename(rows = cols, cols = rows)) %>% 
    distinct()
  # convert into matrix format
  a <- a %>% 
    spread(cols, vcov)
  a_rows <- a$rows
  a$rows <- NULL
  a <- as.matrix(a)
  ## warn if any NAs
  if (any(is.na(a))) warning("Missing values in matrix")
  rownames(a) <- a_rows
  a
}

## Read in metadata about files and join to data ----
sheets <- readxl::excel_sheets("Data/review_csdr_yoda_outcome_overview.xlsx")
sheet_choose <- c("bph", "ciu", "diab", "gout", "gord", "htn", 
                  "ibd", "neuro_psych", "osteoporosis", "pul_htn", "resp", 
                  "rhinitis", "ra_pa", "ps_pa", "ank_spond", "sle")
sheets <- map(sheet_choose, ~ readxl::read_excel("Data/review_csdr_yoda_outcome_overview.xlsx", sheet = .x))
names(sheets) <- sheet_choose
sheets <- map(sheets, ~ map2(.x, names(.x), function(mycol, mycolname) if_else(mycol == 1, mycolname, as.character(mycol))) %>% 
                as_tibble())
sheets_cond <- map(sheets, ~ .x %>% select(nct_id, condition_cb =condition_grp, v1:v6))
sheets_cond <- bind_rows(sheets_cond, .id = "condition_grp")

sheets <- map(sheets, ~ {
  commonvars <- c("nct_id", "sponsor_id_sas", "sponsor_id",    "sponsor", "medicine", "condition")
  allvars <- ncol(.x)
  outvars <- allvars - length(commonvars)
  outvars <- paste0("outcome", 1:outvars)
  names(.x)[!names(.x) %in% commonvars] <- outvars
  .x
})
sheets <- bind_rows(sheets, .id = "condition_grp")
## Shows the number of each group of index conditions
sheets %>% 
  count(condition_grp, sort = T) %>% 
  mutate(cn = cumsum(n),
         cp = cn/(sum(n)))
write_csv(sheets, "Outputs/summarise_med_condition_outcome.csv")

## add condition in to cfs ----
cfs %>% 
  anti_join(sheets %>% distinct(nct_id, condition))

cfs <- cfs %>% 
  inner_join(sheets %>% distinct(nct_id, condition))

## read in notes ----
notes <- readxl::read_xlsx("Data/review_csdr_yoda_outcome_overview.xlsx", sheet = "notes")

## Join condition group to each trial to select,
## Note increased rows as have psoriatic arthropathy in inflammatory arthritis and psoriasis categories  
cfs <- cfs %>% 
  inner_join(sheets %>% distinct(condition_grp, nct_id) %>%
               group_by(nct_id) %>% 
               mutate(n_condition_grps = length(condition_grp)) %>% 
               ungroup())
rm(sheets, sheet_choose)

## Select each outcome ----
# select one outcome per trial
outcom_lst <- list(
  ank_spond = c("BASDAI score","mean_basdai"),
  bph = c("IPSS Total Score"),
  ciu = "DLQI Score",
  diab = c("hba1c", "HbAc"),
  gord = "prcntheatburn_free_days",
  gout = "urate",
  htn =  c(
    "Systolic blood pressure",
    "systolicbp"
    ),
  ibd = c("total_cdai", "total_mayo", "CDAI Score", "MAYO Score"),
  osteoporosis = "BMD Total Hip",
  ps_pa = c("PASI Total score", "total_pasi"),
  ra_pa = c("acrn", "ACR numerical"),
  resp = c("FEV", "FEV1", "FVC"),
  rhinitis = c("Total Nasal Symptom Score"),
  sle = c("SLE Disease Activity Index"),
  neuro_psych = c("final_adas_cog_mci_score", "ADAS score", "total_sib", "RLS Symptom Score Total", "UPDRS Total")
)
## Apply selection so one outcome per trial
## Trials in rheumatoid arthritis and psoratic arthritis would appear
## twice if did not exclude as have ACRN and PASI
cfs <- cfs %>% 
  filter(outcome %in% unlist(outcom_lst)) %>% 
  filter(!(outcome %in% outcom_lst$ra_pa & condition_grp == "ps_pa"),
         !(outcome %in% outcom_lst$ps_pa & condition_grp == "ra_pa"))

## Standardise outcome data ----
## Convert trials onto a standard scale for relevant outcomes
cfs <- cfs %>% 
  mutate(estimate_orig = estimate,
         se_orig = se,
         across(c(estimate, se), ~ {case_when(
    outcome %in% "IPSS Total Score" ~ .x/3,
    outcome %in% "Uroflowmetry Qmax" ~ .x/2,
    outcome %in% "DLQI Score" ~ .x/4,
    ## Note convert from 4 to 2.5 as is in % rather than 
    ## in mmol/mol http://www.ngsp.org/ifcc.asp 
    # NGSP = [0.09148 * IFCC] + 2.152
    outcome %in% c("hba1c", "HbAc") ~ .x/2.5,
    outcome %in% c("Diastolic blood pressure",
                   "diastolicbp",
                   "Systolic blood pressure", 
                   "systolicbp") ~ .x/4,
    outcome %in% c("total_cdai","CDAI Score") ~ .x/70,
    outcome %in% c("total_mayo","MAYO Score") ~ .x/3,
    outcome %in% "SLE Disease Activity Index" ~ .x/4,
    outcome %in% c("BASDAI score","mean_basdai") ~ .x/10,
    ## Note PASI is based on taking the pasi50 (ie a 50% reduction) and converting it into an absolute value by multiplying it
    ## by the mean value at baseline
    outcome %in% c("total_pasi", "PASI Total score") ~ .x/7,
    ## Not ACR-N is based on fact that ACR20 is the smallest change we would normally use
    outcome %in% c("acrn", "ACR numerical") ~ .x/10,
    outcome %in% c("final_adas_cog_mci_score", "ADAS score") ~ .x/3,
    outcome %in% "total_sib" ~ .x/5,
    outcome %in% "RLS Symptom Score Total" ~ .x/3,
    outcome %in% "UPDRS Total" ~ .x/5,
    # note MCID is 120 mL and units are litres, so divide by 1000 before calculation
    outcome %in% c("FEV", "FEV1", "FVC") ~ .x/(120/1000),
    outcome %in% "Total Nasal Symptom Score" ~ .x/0.26,
    outcome %in% "BMD Total Hip" ~ .x/0.056,
    TRUE ~ .x)}))

## Relabel outcomes so same name used across CSDR and YODA trials
rlbl_out$lst <- str_split(rlbl_out$outcome_lbl1, patt = ",")
rlbl_out <- rlbl_out %>%
  unnest(lst) %>% 
  mutate(outcome = str_trim(lst)) %>% 
  select(-lst) %>% 
  distinct()
cfs <- cfs %>% 
  inner_join(rlbl_out) %>% 
  rename(outcome_orig = outcome) %>% 
  rename(outcome = outcome_lbl2) %>% 
  select(-outcome_lbl1) %>% 
  distinct()

## Check which trials are not standardised, 
## makes sense - urate and percent heartburn free days
cfs %>% 
  group_by(nct_id) %>% 
  mutate(allsame = estimate == estimate_orig & se == se_orig) %>% 
  ungroup() %>% 
  filter(allsame, model_n == 1, term == "(Intercept)") %>% 
  count(outcome, nct_id)

## Flip values so all in same direction, ie, higher score is worse disease state
cfs <- cfs %>% 
  mutate(estimate = if_else(outcome %in% c(
    "Uroflowmetry Qmax",
    "Percent heartburn free days",
    "ACR numerical",
    "BMD Total Hip",
    # Note not used in final analysis
    "total_sib",
    "FEV1", "FVC"), 
    estimate*-1, estimate))

## Assign treatment comparisons for each condition group ----
cmpr <- list(
  ank_spond = read_csv("Supporting/ank_spond_metadata.csv"),
  # bph - all 6 trials have the same condition and comparison, placebo verss tadalafil
  # ciu - 3 identical trials, all placebo versus omalizumab
  diab = read_csv("Supporting/diab_metadata.csv"),
  # gord - 2 identical trials, all placebo versus Dexlansoprazole MR
  # gout - 2 trials, placebo versus allopurinol and allopurinol versus FEBUXOSTAT. Put all together despite being different
  gout = read_csv("Supporting/gout_metadata.csv"),
  htn = read_csv("Supporting/htn_metadata.csv"),
  ibd = read_csv("Supporting/ibd_metadata.csv") ,
  osteoporosis = read_csv("Supporting/osteoporosis.csv") ,
  ps_pa = read_csv("Supporting/ps_pa_metadata.csv"),
  ra_pa = read_csv("Supporting/ra_pa_metadata.csv"),
  resp = read_csv("Supporting/resp_metadata.csv"),
  rhinitis = read_csv("Supporting/rhinitis_metadata.csv"),
  # sle - 2 identical trials, all placebo versus Belimumab
  neuro_psych = read_csv("Supporting/neuro_psych_metadata.csv"))
cmpr <- map(cmpr, ~ .x %>% select(nct_id, cmpr)) %>% 
  bind_rows(.id = "condition_grp")
cmpr_remain <- cfs %>% 
  distinct(nct_id, condition_grp) %>% 
  anti_join(cmpr %>% select(nct_id)) %>% 
  distinct() %>% 
  arrange(condition_grp, nct_id)
rmn_lkp <- c("bph" = "pu_g04be",
             "ciu" = "pu_r03dx",
             "gord" = "pu_a02bc",
             "sle" = "pu_l04aa")
cmpr_remain <- cmpr_remain %>% 
  mutate(cmpr = rmn_lkp[condition_grp])
cmpr_tot <- bind_rows(cmpr, cmpr_remain)
nrow(cfs)
cfs <- cfs %>% 
  inner_join(cmpr_tot)
nrow(cfs)

## drop trials with solely within class comparisons, 8 trials
drop_within <- cfs %>%  
  filter(cmpr %in% c("m05ba_m05ba",
                     "A10BJ_A10BJ",
                     "N04BC_N04BC",
                     "C09CA_C09CA",
                     "m04aa_m04aa",
                     ## Note the following is metoprolol versus carvedilol - same MOA
                     "C07AB_C07AG")) %>% 
  distinct(nct_id, cmpr)

cfs <- cfs %>% 
  filter(!nct_id %in% drop_within$nct_id)
rm(cmpr, cmpr_remain, cmpr_tot, drop_within)

## Harmonise condition name within condition groups ----
cond_hrm <- read_csv("Supporting/harmonise_condition_names.csv")
cfs <- cfs %>% 
  inner_join(cond_hrm)
# Rename conditions for psoriasis as occurs in two sets (psoriasis and psoriatic arthropathy)
cfs <- cfs %>% 
  inner_join(sheets_cond %>% select(condition_grp, nct_id, condition_cb))
cfs <- cfs %>% 
  mutate(condition_cb = if_else( (condition_hrm == "psor_art" & outcome == "ACR numerical") | 
                                   condition_cb == "rheumatoid arthritis", "inflammatory arthropathy", condition_cb),
         condition_cb = case_when(condition_grp == "htn" ~ "Hypertension",
                                  condition_grp == "rhinitis" ~ "Rhinitis, allergic",
                                  TRUE ~ condition_cb))

## Select models ----
## Restrict to linear models for this analysis
cfs <- cfs %>% 
  filter(outcome_description %in% 
           c("cb", "cc", 
             "lm_models_acrn_cb_update", 
             "lm_models_acrn_cc_update", 
             "lm_models_cb_update", 
             "lm_models_cc_update", 
             "lm_models_hrtbrn_cb_update", 
             "lm_models_hrtbrn_cc_update"))

## Create strings to identify relevant models from larger set of models
cnt_strng <- "bmi|egfr|fib4|hgb"
cnt_strng_nter <- "arm:bmi|arm:egfr|arm:fib4|arm:hgb|arm:mbp"
bin_string <- "v1|v2|v3|v4|v5|v6"
bin_string_nter <- "arm:v1|arm:v2|arm:v3|arm:v4|arm:v5|arm:v6"
age_string <- "arm:age_std"
sex_string <- "arm:sex"
## note which trials are single sex
singlesextrials <- c("NCT00049829", "NCT00384930", "NCT00439244", "NCT00670501", 
                     "NCT00827242", "NCT00848081", "NCT00855582", "NCT00861757", 
                     "NCT00970632")
## Due to FEV"1" present in model description some model numbers 
## mis-assigned. create new model number for these based on 
## intercept (note all linear models here)
cfs <- cfs %>% 
  mutate(model_n_orig = model_n) %>% 
  group_by(outcome_description, nct_id, outcome) %>% 
  mutate(newmod = if_else(term == "(Intercept)", 1L, 0L),
        model_n = cumsum(newmod)) %>% 
  ungroup() 
cfs <- cfs %>% 
  mutate(model_n = if_else(repo == "yoda", as.integer(model_n_orig), as.integer(model_n)))
## Create a variable with all of the coefficient terms for each model within each trial
cfs <- cfs %>%
  group_by(nct_id, outcome_description, outcome, model_n) %>%
  mutate(covs = term %>% unique() %>% sort() %>% paste(collapse = "|")) %>%
  ungroup()

## Split into subsets of data for different ways of modelling comorbidity ----
## binary comorbidities (ie yes no)
bin <- cfs %>%
  filter(outcome_description %in% c("lm_models_acrn_cb_update", "lm_models_cb_update", 
                                    "lm_models_hrtbrn_cb_update", "cb"))
## Comorbidity as a continuous count
cnt_slctd <- cfs %>%
  filter(outcome_description %in% c("lm_models_acrn_cc_update", "lm_models_cc_update", 
                                    "lm_models_hrtbrn_cc_update","cc"))
## comorbidity as interval data, eg egfr, fib4 etc
intrvl_slctd <- cfs %>%
  filter(outcome_description %in% c("lm_models_acrn_cb_update", "lm_models_cb_update", 
                                    "lm_models_hrtbrn_cb_update", "cb"))
# note, is the same dataset as bin; that is, the age/sex models are held within the binary dataset
agesex_slctd <- bin
rm(cfs)

## Create formulae and count trials for each group of conditions; binary comorbidities  ----
## Identify como bin model adjusted for age and sex
bin1 <- bin %>% 
  filter(str_detect(covs, bin_string_nter), 
         !str_detect(covs, cnt_strng),
         str_detect(covs, age_string),
         str_detect(covs, sex_string),
         !nct_id %in% singlesextrials)
## Identify como bin model adjusted for age, single sex trials
bin2 <- bin %>% 
  filter(str_detect(covs, bin_string_nter), 
         !str_detect(covs, cnt_strng),
         str_detect(covs, age_string),
         nct_id %in% singlesextrials)
bin <- bind_rows(bin1, bin2)

## Create file indicating different model types and remove this as allocate model numbers 
trnst <- bin %>% 
  filter(term == "(Intercept)") %>% 
  arrange(nct_id, outcome)  %>% 
  distinct(repo, nct_id, outcome, model_n) %>% 
  count(repo, nct_id, outcome, model_n) %>% spread(model_n, n)
## Note, all of YODA trials are model_n = 25
bin25 <- trnst %>% 
  filter(!is.na(`25`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 25)
trnst <- trnst %>% 
  anti_join(bin25 %>% select(nct_id, outcome))
bin5 <- trnst %>% 
  filter(!is.na(trnst$`5`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 5)
trnst <- trnst %>% 
  anti_join(bin5 %>% select(nct_id, outcome))
bin_slct <- bind_rows(bin25,
                      bin5,
                      .id = "chosen") %>% 
  distinct(nct_id, .keep_all = TRUE)
bin_slctd <- bin %>% 
  semi_join(bin_slct)
rm(bin1, bin2, bin25, bin5, bin, bin_slct, trnst)

## Limit terms to comorbidity (binary) treatment interactions
bin_slctd <- bin_slctd %>% 
  filter(str_detect(term, bin_string_nter))

## Fix transposition of comorbidity for one COPD trial where this was an issue
## - note only need to do this for comorbidity as binary
reverse <- c("arm:v6" = "arm:v5", "arm:v5" = "arm:v6")
bin_slctd <- bin_slctd %>% 
  mutate(term = if_else(nct_id == "NCT01772134" & term %in% reverse, reverse[term], term)) %>% 
  arrange(condition_cb, nct_id, term)

## Create formulae and count trials for each group of conditions; comorbidity count ----
# Limit continuous trials to trials with binary comorbidites
xcld_no_cb <- cnt_slctd %>% 
  filter(!nct_id %in% bin_slctd$nct_id) %>% 
  distinct(nct_id) %>% 
  arrange(nct_id)
cnt_slctd <- cnt_slctd %>% 
  filter(nct_id %in% bin_slctd$nct_id)
## Select models adjusting for age and sex
cnt1 <- cnt_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, "arm:como_cnt"),
         str_detect(covs, age_string),
         str_detect(covs, sex_string),
         !nct_id %in% singlesextrials)
## Select models adjusting for age. Single sex trials
cnt2 <- cnt_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, "arm:como_cnt"),
         str_detect(covs, age_string),
         nct_id %in% singlesextrials)
cnt <- bind_rows(cnt1, cnt2)

## Pull comorbidity count squared models
## Select models adjusting for age and sex
cnt1_sq <- cnt_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, "arm\\:I\\(como_cnt\\^2\\)"),
         str_detect(covs, age_string),
         str_detect(covs, sex_string),
         !nct_id %in% singlesextrials)
## select models adjusting for age, single sex
cnt2_sq <- cnt_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, "arm\\:I\\(como_cnt\\^2\\)"),
         str_detect(covs, age_string),
         nct_id %in% singlesextrials)
cnt_sq <- bind_rows(cnt1_sq, cnt2_sq)

## Create file indicating different model types and remove this as allocate model numbers 
trnst <- cnt %>% 
  filter(term == "arm:como_cnt") %>% 
  arrange(nct_id, outcome)  %>% 
  distinct(repo, nct_id, outcome, model_n) %>% 
  count(repo, nct_id, outcome, model_n) %>% spread(model_n, n)
cnt8 <- trnst %>% 
  filter(!is.na(`8`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 8)
trnst <- trnst %>% 
  anti_join(cnt8 %>% select(nct_id, outcome))
cnt_slct <- bind_rows(cnt8,
                      .id = "chosen") %>% 
  distinct(nct_id, .keep_all = TRUE)
cnt_slctd <- cnt %>% 
  semi_join(cnt_slct)

trnst <- cnt_sq %>% 
  filter(term == "arm:I(como_cnt^2)") %>% 
  arrange(nct_id, outcome)  %>% 
  distinct(repo, nct_id, outcome, model_n) %>% 
  count(repo, nct_id, outcome, model_n) %>% spread(model_n, n)
## all captured by 9
cnt_sq_slctd <- cnt_sq %>% 
  filter(model_n == 9, term %in% c("arm","arm:como_cnt", "arm:I(como_cnt^2)"))
rm(cnt8, cnt1, cnt2, cnt, cnt_slct, trnst)

## Create formulae and count how many trials for each group of conditions; 
## interval measures, egfr, fib4 etc ----
# same 9 excluded with no binary comorbidities
xcld_no_cb2 <- intrvl_slctd %>% 
  filter(!nct_id %in% bin_slctd$nct_id) %>% 
  distinct(nct_id) %>% 
  arrange(nct_id)
intrvl_slctd <- intrvl_slctd %>% 
  filter(nct_id %in% bin_slctd$nct_id)

## Identify interval variables model adjusted for age and sex
intrvl1 <- intrvl_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, cnt_strng_nter),
         str_detect(covs, age_string),
         str_detect(covs, sex_string),
         !nct_id %in% singlesextrials)
## Identify interval variables model adjusted for age. Single sex trials
intrvl2 <- intrvl_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         str_detect(covs, cnt_strng_nter),
         str_detect(covs, age_string),
         nct_id %in% singlesextrials)
intrvl <- bind_rows(intrvl1, intrvl2)
trnst <- intrvl %>% 
  filter(term == "(Intercept)") %>% 
  arrange(nct_id, outcome)  %>% 
  distinct(nct_id, outcome, model_n) %>% 
  count(nct_id, outcome, model_n) %>% spread(model_n, n)
intrvl27 <- trnst %>% 
  filter(!is.na(`27`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 27)
trnst <- trnst %>% 
  anti_join(intrvl27 %>% select(nct_id, outcome))
intrvl_slct <- bind_rows(intrvl27,
                      .id = "chosen")
intrvl_slctd <- intrvl %>% 
  inner_join(intrvl_slct %>% select(-chosen))
rm(intrvl27, intrvl, intrvl_slct, trnst)

## Limits to continuous treatment interactions
intrvl_slctd <- intrvl_slctd %>% 
  filter(str_detect(term, cnt_strng_nter))

## all in bin_slctd
intrvl_slctd %>% 
  filter(!nct_id %in% bin_slctd$nct_id)
## 12 no lab data
bin_slctd %>% 
  filter(!nct_id %in% intrvl_slctd$nct_id) %>% 
  distinct(nct_id, condition_cb) %>% 
  arrange(condition_cb, nct_id)

## Create formulae and count trials for each group of conditions; agesex ----
xcld_no_cb3 <- agesex_slctd %>% 
  filter(!nct_id %in% bin_slctd$nct_id) %>% 
  distinct(nct_id) %>% 
  arrange(nct_id)
# should be TRUE
identical(xcld_no_cb3, xcld_no_cb)
agesex_slctd <- agesex_slctd %>% 
  filter(nct_id %in% bin_slctd$nct_id)
## pull models for trials which are not single sex
agesex1 <- agesex_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         !str_detect(covs, cnt_strng),
         str_detect(covs, age_string),
         str_detect(covs, sex_string),
         !nct_id %in% singlesextrials)
## pull models for trials which are single sex
agesex2 <- agesex_slctd %>% 
  filter(!str_detect(covs, bin_string), 
         !str_detect(covs, cnt_strng),
         str_detect(covs, age_string),
         nct_id %in% singlesextrials)
agesex <- bind_rows(agesex1, agesex2)
trnst <- agesex %>% 
  filter(term == "(Intercept)") %>% 
  arrange(nct_id, outcome)  %>% 
  distinct(repo, nct_id, outcome, model_n) %>% 
  count(repo, nct_id, outcome, model_n) %>% spread(model_n, n)
## Note, all of YODA trials are model_n = 30
agesex30 <- trnst %>% 
  filter(!is.na(`30`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 30)
trnst <- trnst %>% 
  anti_join(agesex30 %>% select(nct_id, outcome))
agesex6 <- trnst %>% 
  filter(!is.na(trnst$`6`)) %>% 
  distinct(nct_id, outcome) %>% 
  mutate(model_n = 6)
trnst <- trnst %>% 
  anti_join(agesex6 %>% select(nct_id, outcome))
agesex_slct <- bind_rows(agesex30,
                         agesex6,
                         .id = "chosen") %>% 
  distinct(nct_id, .keep_all = TRUE)
agesex_slctd <- agesex %>% 
  semi_join(agesex_slct)
rm(agesex30, agesex6, agesex_slct, trnst)

## Limits to arm and arm:sex and arm:age treatment interactions
agesex_slctd <- agesex_slctd %>% 
   filter(term %in% c("arm", "arm:age_std", "arm:sex"))
## all in bin_slctd
agesex_slctd %>% 
  filter(!nct_id %in% bin_slctd$nct_id)
## No missing age ones
bin_slctd %>% 
  filter(!nct_id %in% agesex_slctd$nct_id) %>% 
  distinct(nct_id, condition_cb, repo) %>% 
  arrange(condition_cb, nct_id)

## Drop all other objects and bind model contents together 
rm(age_string, 
   agesex, agesex1, agesex2, 
   bin_string, bin_string_nter,
   cnt_strng, cnt_strng_nter, cond_hrm,
  # condition_grp_name, 
  intrvl1, intrvl2,
  notes, outcom_lst, 
  # outcome_name, 
  reverse, rlbl_out, rmn_lkp,
  sex_string, sheets_cond, singlesextrials, xcld_no_cb,
  xcld_no_cb2, xcld_no_cb3)

## create single dataset with all selected models ----
allcfs <- bind_rows(agesex = agesex_slctd,
                bin  = bin_slctd, 
                cnt = cnt_slctd, 
                cnt_sq = cnt_sq_slctd,
                intrvl = intrvl_slctd,
                .id = "model_type")
rm(agesex_slctd,
   bin_slctd,
   cnt_slctd,
   cnt_sq_slctd,
   cnt_sq, cnt1_sq, cnt2_sq,
   intrvl_slctd)

## Write model formulae based on number of trials ----
## Count trials within each to see which model formulae needed
forform <- allcfs %>% 
  distinct(model_type,
           nct_id,
           condition_cb,
           cmpr) %>%  
  group_by(model_type,
        condition_cb,
        cmpr) %>% 
  summarise(trial = sum(!duplicated(nct_id))) %>% 
  ungroup()
## Each model is for a single treatment comparison/index condition combination.
## for conditions with fewer than 5 trial, treat everything as fixed, otherwise treat it as random
## Otherwise nest trial within treatment comparison
forform <- forform %>% 
  mutate(formtype = if_else(trial < 5, "simple", "complex")) 
forform$myform <- map(forform$formtype, function(formtype_choose) {
 if (formtype_choose == "simple") {
   estimate ~ 0 + term  + fcor(covs)
 } else {
   estimate ~ 0 + term  + (0 + term | nct_id) + fcor(covs) 
 }
})

## Add in covariance matrices ----
## csdr
covs <- map(list("Data/covs_csdr1.csv.gz",
                 "Data/covs_csdr2.csv.gz",
                 "Data/covs_csdr3.csv.gz",
                 "Data/covs_csdr4.csv.gz"), read_csv)
covs <- bind_rows(covs)
covs <- covs %>% 
  filter(outcome_description %in% c(
    "lm_models_acrn_cb_update", 
    "lm_models_acrn_cc_update", 
    "lm_models_cb_update", 
    "lm_models_cc_update", 
    "lm_models_hrtbrn_cb_update", "lm_models_hrtbrn_cc_update"))
covs <- covs %>% 
  mutate(model_n_orig = model_n) %>% 
  group_by(outcome_description, nct_id, outcome) %>% 
  mutate(newmod = if_else(rows == "(Intercept)" & cols == "(Intercept)", 1L, 0L),
         model_n = cumsum(newmod)) %>% 
  ungroup() 
covs <- covs %>% 
  rename(outcome_orig = outcome)
covs <- covs %>% 
  semi_join(allcfs %>% distinct(nct_id, outcome_description, model_n_orig, outcome_orig, model_n)) 
covs <- covs %>% 
  inner_join(allcfs %>% distinct(outcome_orig, outcome)) 

## convert names of rows and columns of covariance matrix so that they match the coefficient terms
Compareterms <- function(){
  print(paste0("In rows not cols ", paste(setdiff(covs$rows, covs$cols), collapse = ", ")))
  print(paste0("In rows not coefficient terms ", setdiff(covs$rows, allcfs$term) %>% paste(collapse = ", ")))
  # print(paste0("In rows not coefficient terms just CSDR ", setdiff(covs$rows, allcfs$term[allcfs$repo == "csdr"])%>% paste(collapse = ", ")))
  print(paste0("In coefficient terms not rows ", setdiff(allcfs$term, covs$rows) %>% paste(collapse = ", ")))
}
Compareterms()
covs <- covs %>% 
  mutate(across(c(rows, cols), ~ .x %>% 
                  str_remove_all("TRUE")),
         across(c(rows, cols), ~  
                case_when(
                  .x == "como_cnt:arm" ~ "arm:como_cnt",
                  .x == "sex:arm" ~ "arm:sex",
                  .x == "arm:anaem" ~ "anaem:arm",
                  .x == "strt" ~ "base",
                  TRUE ~ .x)
  ))
Compareterms()
covs <- covs %>% 
  semi_join(allcfs %>% select(nct_id, outcome_description, model_n, outcome, rows = term)) %>% 
  semi_join(allcfs %>% select(nct_id, outcome_description, model_n, outcome, cols = term))
## YODA
covs2 <- read_csv("Data/covs_yoda.csv.gz")
setdiff(allcfs$outcome_orig[allcfs$repo == "yoda"], covs2$outcome_orig)
intersect(allcfs$outcome_orig[allcfs$repo == "yoda"], covs2$outcome_orig)
covs2 <- covs2 %>% 
  semi_join(allcfs %>% distinct(nct_id, outcome_description, model_n_orig, outcome_orig, outcome, model_n)) 
covs2 <- covs2 %>% 
  inner_join(allcfs %>% distinct(nct_id, outcome_description, model_n_orig, outcome_orig, outcome, model_n)) 
covs2 <- covs2 %>% 
  mutate(across(c(rows, cols), ~ .x %>% 
                  str_replace_all("armtreat", "arm") %>% 
                  str_replace_all("sexTRUE", "sex") %>% 
                  str_remove_all("TRUE$")),
         across(c(rows, cols), ~
                case_when(
                  .x == "como_cnt:arm" ~ "arm:como_cnt",
                  .x == "sex:arm" ~ "arm:sex",
                  .x == "strt" ~ "base",
                  TRUE ~ .x)
  ))

covs2 <- covs2 %>% 
  semi_join(allcfs %>% select(nct_id, model_n, outcome_orig, rows = term)) %>% 
  semi_join(allcfs %>% select(nct_id, model_n, outcome_orig, cols = term))

covs <- covs[ ,intersect(names(covs), names(covs2))]
covs2 <- covs2[ ,intersect(names(covs), names(covs2))]
covs <- bind_rows(csdr = covs,
                  yoda = covs2, .id = "repo") 
rm(covs2)

## Apply correction factor to covariance matrix ----
## Need to apply to variances. Need to square the correction as it is for the
## standard error. 
## small set of model/trial combos which did not converge 
## (model 25, binary predictors, 7 yoda trials)
allcfs <- allcfs %>% 
  filter(!is.na(estimate), !is.na(se))
covs <- covs %>% 
  filter(!is.na(vcov))
crct <- allcfs %>% 
  filter(!estimate == 0) %>% 
  mutate(crct = round(estimate/estimate_orig, 5)) %>% 
  distinct(nct_id, outcome_orig, crct)
covs_out <- covs %>% 
  anti_join(crct)
covs <- covs %>% 
  inner_join(crct)
rm(crct, covs_out)
covs <- covs %>% 
  mutate(vcov = vcov * crct^2)

## Correct continuous variables by dividing by the standard deviation ----
## Multiply by SD to get onto a common scale
## Note that pre-analysis divided BMI by 10 and eGFR, MBP divided by 100. Fib4 and hgb not modified in CSDR
## hgb is is modified in yoda. Therefore group by repo
cnt_sd <- read_csv("Data/Baseline_m_s_csdr.csv")
cnt_sd_yoda <- read_csv("Data/Baseline_m_s_yoda.csv")
arms <- read_csv("Data/arms.csv")
arms <- arms %>% 
  semi_join(allcfs)
cnt_sd <- cnt_sd %>% 
  semi_join(allcfs) %>% 
  inner_join(arms %>% mutate(n = control + treat) %>% select(nct_id, n))
setdiff(names(cnt_sd), names(cnt_sd_yoda))
setdiff(names(cnt_sd_yoda), names(cnt_sd))
cnt_sd_yoda <- cnt_sd_yoda[,names(cnt_sd)]
cnt_sd <- bind_rows(csdr = cnt_sd,
                    yoda = cnt_sd_yoda,
                    .id = "repo")
cnt_sd <- cnt_sd %>% 
  group_by(repo) %>% 
  summarise(across(c(bmi_s, egfr_s, fib4_s, hgb_s, mbp_s), ~ weighted.mean(.x, n, na.rm = TRUE))) %>% 
  ungroup()
names(cnt_sd) <- str_remove(names(cnt_sd),"_s$")
cnt_sd <- cnt_sd %>% 
  gather("term", "crct_x", -repo) %>% 
  mutate(term = paste0("arm:", term))
allcfs <- allcfs %>% 
  left_join(cnt_sd) %>% 
  mutate(crct_x = if_else(is.na(crct_x), 1, crct_x),
         across(c(estimate, se), ~ .x*crct_x))
covs <- covs %>% 
  left_join(cnt_sd %>% rename(rows = term, crct_r = crct_x)) %>% 
  left_join(cnt_sd %>% rename(cols = term, crct_c = crct_x)) %>% 
  mutate(crct_r = if_else(is.na(crct_r), 1, crct_r),
         crct_c = if_else(is.na(crct_c), 1, crct_c),
         vcov = vcov*(crct_r*crct_c))
rm(cnt_sd, cnt_sd_yoda, arms)

## create single dataset with estimates and covariance matrix ----
allcfs_nst <- allcfs %>% 
  select(model_type, repo, outcome_description, term, estimate, se, nct_id, model_n, model_n_orig, outcome, condition, condition_cb, cmpr, condition_hrm) %>% 
  nest(cfs = c(term, estimate, se))
covs_nst <- covs %>% 
  select(-crct, -crct_r, -crct_c) %>% 
  nest(cvs = c(rows, cols, vcov))
rm(covs)

# none with duplicates, should be false
any(map_lgl(covs_nst$cvs, ~ nrow(.x) != nrow(.x %>% distinct(rows, cols))))
## convert dataframe to a matrix
covs_nst$cvs <- map(covs_nst$cvs, CnvrtMatrix)
allcfs_nst <- allcfs_nst %>% 
  left_join(covs_nst)
rm(covs_nst)
allcfs_nst$nocov <- map_lgl(allcfs_nst$cvs, is.null)
## As expected no covariance for yoda with comorbidity count alone (not exported)
## is present for como count squared
allcfs_nst %>% 
  group_by(repo, outcome_description, nocov) %>% 
  summarise(n = n(),
            model_ns = paste(model_n %>% unique() %>% sort(), collapse = ";")) %>% 
  ungroup()
allcfs_nst <- allcfs_nst %>% 
  filter(!outcome_description == "cc" | model_n == 9,
         !(repo == "yoda" & outcome == "ACR numerical" & outcome_description == "cc"))

## Check missingness; should be zero for both
map_int(allcfs_nst$cfs, ~ sum(is.na(.x$estimate) | is.na(.x$se))) %>% table()
map_int(allcfs_nst$cvs, ~ sum(is.na(.x))) %>% table()

# check if any missingness in matrices, should be FALSE as ignores null matrices
map_lgl(allcfs_nst$cvs, ~ any(is.na(.x))) %>% any()

## List variables and check all aligned and that matrix is positive definite, 
## all should be TRUE
allcfs_nst$cfs <- map(allcfs_nst$cfs, ~ .x %>% arrange(term))
map(allcfs_nst$cfs, ~ .x$term) %>% unlist() %>% unique() 
map(allcfs_nst$cvs, rownames) %>% unlist() %>% unique()
map2_lgl(allcfs_nst$cvs, allcfs_nst$cfs, ~ all(.x == .x[.y$term, .y$term])) %>% all()
allcfs_nst$cvs <- map2(allcfs_nst$cvs, allcfs_nst$cfs, ~ .x[.y$term, .y$term])
map_lgl(allcfs_nst$cvs[allcfs_nst$nocov == FALSE], matrixcalc::is.positive.definite) %>% all

## Calculate effect across men and women based on variances and covariance ----
agesex <- allcfs_nst %>% 
  filter(model_type == "agesex")
agesex$sex <- if_else(map_int(agesex$cfs, nrow) == 2, "single", "both")
agesexboth <- agesex
agesexboth$cov <- map_dbl(agesexboth$cvs, ~ if(nrow(.x) ==3 ) {.x["arm:sex", "arm"]} else 0)
agesexboth$cfs <- map(agesexboth$cfs, ~ .x %>% 
                            filter(term %in% c("arm", "arm:sex")))
agesexboth$cfs <- map(agesexboth$cfs, ~ .x %>% 
                        mutate(term = str_remove(term, "\\:")) %>% 
                        pivot_wider(names_from = term, 
                                    values_from = c(estimate, se)))
agesexboth <- agesexboth %>% 
  select(nct_id, outcome, condition_cb, cmpr, cfs, cov) %>% 
  unnest(cfs)
# set estimate to zero where is missing and set se to Infinity;
## coding note cannot do this within pivot wider (above) as is not "missing"
## until unnested. Do this to deals with single sex trials.
agesexboth <- agesexboth %>% 
  mutate(estimate_armsex = if_else(is.na(estimate_armsex), 0, estimate_armsex),
         se_armsex = if_else(is.na(se_armsex), Inf, se_armsex))  %>% 
  arrange(nct_id)
## calculate effects
agesexboth <- agesexboth %>% 
  mutate(est_a = estimate_arm,
         est_b = estimate_arm + estimate_armsex,
         var_a = se_arm^2,
         var_b = se_arm^2 + se_armsex^2 + 2*cov,
         wt_a = 1/var_a,
         wt_b = 1/var_b,
         est_mn = (est_a*wt_a + est_b * wt_b) / (wt_a + wt_b),
         var_mn  = 1/(wt_a + wt_b),
         se_mn = var_mn^0.5)
## examine what happens when only one sex - correctly, as expected, same results as arm alone
agesexboth %>% filter(cov == 0) %>% select(estimate_arm, est_mn, se_arm, se_mn)

## nest data so one row per index condition/comparison combination
agesexboth <- agesexboth %>% 
  select(nct_id, outcome, condition_cb, cmpr, est_mn, se_mn) %>% 
  nest(data = c(nct_id, est_mn, se_mn))
saveRDS(agesexboth, "Scratch_data/bothsexdata.Rds")
rm(agesexboth, agesex)

## Re-organise data into groups of conditions and treatment comparisons ----
allcfs_nst <- allcfs_nst %>% 
  filter(!is.na(condition_cb)) %>% 
  arrange(model_type, condition_cb, cmpr) %>% 
  group_by(model_type, condition_cb, cmpr) %>% 
  nest() %>% 
  ungroup()

## Create block diagonal matrices for fitting within brms ----
CreateDataLong <- function(a) {
  a <- a %>% 
    arrange(nct_id)
  cvs <- a %>% 
    select(nct_id, cvs)
  bd <- Matrix::bdiag(cvs$cvs)
  b <- a %>% 
    select(-cvs) 
  list(cfs = b %>% unnest(cfs),
       bd = bd)
}
## Separate out models with comorbididty as a count and the rest. 
## This is because como as a count only needs the standard error
## as only one parameter of interest for the meta-analysis
allcfs_nst$models <- map_chr(allcfs_nst$data, ~ .x$model_n %>%
                               unique() %>%
                               sort() %>% 
                               paste(collapse = ",") )
cnt <- allcfs_nst %>% 
  filter(model_type == "cnt")
oth <- allcfs_nst %>% 
  filter(!model_type %in% c("cnt"))
oth$data <- map(oth$data, CreateDataLong)
oth$cfs <- map(oth$data, ~ .x[[1]])
oth$bd <- map(oth$data, ~ .x[[2]])
oth$data <- NULL
oth <- oth %>% 
  ungroup()
## check again that block diagonal covariance matrix is positive definite 
## (ie no error in conversion) should be TRUE
map_lgl(oth$bd, ~ matrixcalc::is.positive.definite(as.matrix(.x))) %>% all
## For comorbidity count, only need estimate and se
cnt <- cnt %>% 
  unnest(data)
cnt$cfs <- map(cnt$cfs, ~ .x %>% filter(term == "arm:como_cnt"))
cnt <- cnt %>% 
  unnest(cfs)
cnt <- cnt %>% 
  group_by(model_type, condition_cb, cmpr) %>% 
  nest() %>% 
  ungroup()
## Join formulae list onto data list ----
## need to join using comparison as well as index condition
forformcc <- forform %>% 
  filter(model_type == "cnt")
forformcc$myform[forformcc$formtype == "simple"]  <- map(forformcc$myform[forformcc$formtype == "simple"], function(x) brms::bf(estimate|se(se) ~ 1))
forformcc$myform[!forformcc$formtype == "simple"] <- map(forformcc$myform[!forformcc$formtype == "simple"], function(x) brms::bf(estimate|se(se) ~ 1 + (1|nct_id)))
cnt <- cnt %>% 
  inner_join(forformcc %>% select(model_type, condition_cb, cmpr, formtype, myform))
oth <- oth %>% 
  inner_join(forform %>% select(model_type, condition_cb, cmpr, formtype, myform))
## Drop all objects except model dataframes and formulae and summary tables
rm(allcfs_nst, forform, forformcc,
   CnvrtMatrix, Compareterms, CreateDataLong)

