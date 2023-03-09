# model_post
library(tidyverse)
library(brms)
library(tidybayes)

## continuous outcomes ----
## original and additional processing requested by reviewers
oth_or <- readRDS("Scratch_data/mdls_mvn.Rds")
oth_or$models <- map_chr(oth_or$cfs, ~ paste(.x$model_n %>% sort() %>% unique(), collapse = ";"))
oth_sq <- readRDS("Scratch_data/mdls_mvn_como_sq.Rds")
cnt_or <- readRDS("Scratch_data/mdls_cnt.Rds")
cnt_or$models <- map_chr(cnt_or$data, ~ paste(.x$model_n %>% sort() %>% unique(), collapse = ";"))
cnt_re <- readRDS("Scratch_data/mdls_cnt_re.Rds")
cnt_re$models <- map_chr(cnt_re$data, ~ paste(.x$model_n %>% sort() %>% unique(), collapse = ";"))

## index specific effect estimates ----
slct <- c("models", "model_type", "condition_cb", "cmpr", "formtype", "res")
bth <- bind_rows(oth_or = oth_or %>% select(all_of(slct)), 
                 oth_sq = oth_sq %>% select(all_of(slct)),
                 cnt_or = cnt_or %>% select(all_of(slct)),
                 cnt_re = cnt_re %>% select(all_of(slct)),
                 .id = "model_grp")
rm(oth_or, oth_sq, cnt_or, cnt_re)

## not actually meta-analysis, these are the single trials
oth_sngl <- readRDS("Scratch_data/mdls_mvn_sngl_trials.Rds") %>% 
  mutate(models = "single_trial")

## Labels 
cond_lbl <- read_csv("Outputs/lookup_cond_lbl.csv")
cmpr_lbl <- read_csv("Outputs/lookup_tx_compar_lbl.csv")

## Process single trials ----
oth_sngl <- oth_sngl %>% 
  select(models, model_type, condition_cb, cmpr, cfs) %>% 
  unnest(cfs)
oth_sngl <- oth_sngl %>% 
  mutate(`l-95% CI` = estimate - 1.96*se,
         `u-95% CI` = estimate + 1.96*se) %>% 
  select(models, model_type, condition_cb, cmpr, term, Estimate = estimate, Est.Error = se, `l-95% CI`, `u-95% CI`)

## combine meta and single trials then reformat ----
bth$fixed <- map(bth$res, ~ (summary(.x)$fixed) %>% as_tibble(rownames = "term"))
bth_smry <- bth %>% 
  select(models, model_grp, model_type, condition_cb, cmpr, fixed) %>% 
  unnest(fixed)
bth_smry <- bind_rows(sngl = oth_sngl %>% mutate(model_grp = "single_trial"),
                     meta = bth_smry, .id = "modelled")
bth_smry <- bth_smry %>% 
  mutate(term = str_remove(term, "^term")) %>% 
  rename(Q2.5 = `l-95% CI`,
         Q97.5 = `u-95% CI`) %>% 
  mutate(across(c(Estimate, Q2.5, Q97.5), ~ round(.x, 2) %>% formatC(digits = 2, format = "f", flag = "0"))) %>%
  mutate(res = paste0(Estimate, " (", Q2.5, " to ", Q97.5, ")"),
         res = if_else(str_detect(res, "^\\-"), res, paste0(" ",res))) %>% 
  select(-Rhat, -Bulk_ESS, -Tail_ESS)

## T2 ----
t2 <- bth_smry %>% 
  filter(!model_grp == "oth_sq", !model_grp == "cnt_re") %>% 
  select(-modelled) %>% 
  mutate(res = if_else(Q2.5 > 0 | Q97.5 < 0, paste0(res, "*"), res)) %>% 
  filter(model_type %in% c("agesex", "cnt"),
         term %in% c("arm", "arm:age_std", "arm:sex", "arm:como_cnt", "Intercept"))  %>% 
  select(condition_cb, cmpr, term, res)  %>% 
  spread(term, res, fill = " ") %>% 
  select(-Intercept, -arm) %>% 
  rename(
         `Age treatment interaction` = `arm:age_std`,
         `Sex treatment interaction` = `arm:sex`)

## Binary comorbidities ----
bin <- bth_smry %>% 
  mutate(star = if_else( sign(as.double(Q2.5)) == sign(as.double(Q97.5)), 1L, 0L)) %>% 
  filter(model_type == "bin") %>% 
  mutate(term = str_remove(term, "^arm:")) %>% 
  mutate(res = if_else(star == 1, paste0(res, "*"), res)) %>% 
  select(modelled , condition_cb, cmpr, term, res)
## add comorbidity labels
como_names <- readRDS(file = "../definitive_exports_csdr_yoda/Processed_metadata/comorbidity_top6_csdr.Rds") %>% 
  mutate(condition_cb = str_to_lower(condition_grp))
yoda_names <- readRDS("../definitive_exports_csdr_yoda/Processed_metadata/comorbidity_top6_yoda.Rds")
como_names_lng <- como_names %>% 
  distinct(condition_cb, .keep_all = TRUE) %>% 
  select(condition_cb, v1:v6) %>% 
  gather("term", "como", -condition_cb) 

## Drop inflammatory arthropathy as different 6 top conditions
bin <- bin %>% 
  mutate(condition_cb_orig = condition_cb,
         condition_cb = str_to_lower(condition_cb)) %>% 
  filter(!condition_cb == "inflammatory arthropathy")
setdiff(bin$condition_cb, como_names$condition_cb)
bin <- bin %>% 
  inner_join(como_names_lng) 
bin_tbl <- bin %>% 
 select(condition_cb = condition_cb_orig, cmpr, como, res) %>% 
 spread(como, res, fill = " ")

## continuous predictors ----
cnt_tbl <- bth_smry %>% 
  mutate(star = if_else( sign(as.double(Q2.5)) == sign(as.double(Q97.5)), 1L, 0L)) %>% 
  filter(model_type == "intrvl") %>% 
  mutate(term = str_remove(term, "^arm:") %>% str_to_upper()) %>% 
  mutate(res = if_else(star == 1, paste0(res, "*"), res)) %>% 
  select(condition_cb, cmpr, term, res) %>% 
  spread(term, res, " ")
cnt_tbl

## Relabel conditions and comparisons
Relbl <- function(x) x %>% 
  inner_join(cond_lbl) %>% 
  select(-condition_cb) %>% 
  inner_join(cmpr_lbl) %>% 
  select(-cmpr, -cmpr2) %>% 
  select(`Index condition` = condition_cb2, `Treatment comparison` = cmpr3, everything()) 

t2 <- Relbl(t2)
bin_tbl <- Relbl(bin_tbl)
cnt_tbl <- Relbl(cnt_tbl)

write_csv(t2, "Outputs/table2.csv", na = " ")
write_csv(bin_tbl, "Outputs/table_binary_predictors.csv", na = " ")
write_csv(cnt_tbl, "Outputs/table_continous_predictors.csv", na = " ")

## Non continuous outcomes ----
dfs <- readRDS("Scratch_data/cfs_non_cont.Rds")
dfs <- map(dfs, ~ .x[, intersect(names(.x), c("nct_id", "outcome", "medicine", "control","cmpr", "condition"))] %>% distinct())
dfs$tbe_bleed$outcome <- "bleeding"
dfs$tbe_dvtpe$outcome <- "dvt/dvt or pe"
dfs$osteop <- dfs$osteop %>% 
  mutate(cmpr = case_when(medicine == "Teriparatide" ~ "pu_H05AA02", 
                          medicine == "zoledronic acid" ~ "pu_M05BA08"),
         control = "placebo")
dfs$hdch <- dfs$hdch %>% 
  mutate(cmpr = "pu_N03AX11",
         control = "placebo")
dfs$osteop$outcome <- "Fracture of VERT"
dfs$osteop$condition <- "Osteoporosis"

dfs <- map(dfs, ~ .x %>% 
  group_by(outcome, cmpr, condition, medicine, control) %>% 
  summarise(n_trials = sum(!duplicated(nct_id))) %>% 
  ungroup())

# if treated lower risk of headache, comorbidity and comorbidity-treatment interaction null
hdch <- readRDS("Scratch_data/headache_models.Rds")
hdch$smry <- map(hdch$res, ~ (summary(.x)$fixed) %>% as_tibble(rownames = "term")) 
hdch_smry <- hdch %>% 
  ungroup() %>% 
  select(model_n, smry) %>% 
  unnest(smry)
hdch_smry <- bind_cols(hdch_smry, 
  dfs$hdch) 
rm(hdch)
## if treated lower risk of fracture.
oste <- readRDS("Scratch_data/osteop_models.Rds")
oste$smry <- map2(oste$res, oste$data, ~ {
  if (nrow(.y) == 1) {
    .y %>% mutate(term = "Intercept") %>%  select(term, Estimate = estimate, Est.Error = std.error, `l-95% CI`= conf.low, `u-95% CI` = conf.high)
    } else{
    (summary(.x)$fixed) %>% as_tibble(rownames = "term")
  }
})
oste_smry <- oste %>% 
  ungroup() %>% 
  select(model_n, smry, medicine) %>% 
  unnest(smry)
oste_smry <- oste_smry %>% 
  inner_join(dfs$osteop)
rm(oste)
## tbe bleed
bld <- readRDS("Scratch_data/tbe_bleed_models.Rds")
bld$smry <- map2(bld$res, bld$data, ~ {
  if (nrow(.y) == 1) {
    .y %>% mutate(term = "Intercept",
                  `l-95% CI` = estimate - 1.96*std.error,
                  `u-95% CI` = estimate + 1.96*std.error) %>%
      select(term, Estimate = estimate, Est.Error = std.error, `l-95% CI`, `u-95% CI`)
  } else{
    (summary(.x)$fixed) %>% as_tibble(rownames = "term")
  }
})

bld_smry <- bld %>% 
  ungroup() %>% 
  select(model_n, smry, cmpr, condition) %>% 
  unnest(smry)
bld_smry <- bld_smry %>% 
  inner_join(dfs$tbe_bleed)
rm(bld)

## tbe pe/dvt
pedvt <- readRDS("Scratch_data/tbe_dvtpe_models.Rds")
pedvt$smry <- map2(pedvt$res, pedvt$data, ~ {
  if (nrow(.y) == 1) {
    .y %>% mutate(term = "Intercept",
                  `l-95% CI` = estimate - 1.96*std.error,
                  `u-95% CI` = estimate + 1.96*std.error) %>%
      select(term, Estimate = estimate, Est.Error = std.error, `l-95% CI`, `u-95% CI`)
  } else{
    (summary(.x)$fixed) %>% as_tibble(rownames = "term")
  }
})

pedvt_smry <- pedvt %>% 
  ungroup() %>% 
  select(model_n, smry, cmpr, condition) %>% 
  unnest(smry)
pedvt_smry <- pedvt_smry %>% 
  inner_join(dfs$tbe_dvtpe)
rm(pedvt)

## all_non_cnt
non_cnt <- bind_rows(tbe_pedvt = pedvt_smry,
                     tbe_bld = bld_smry,
                     hdch = hdch_smry,
                     oste = oste_smry,
                     .id = "cond_grp")


non_cnt <- non_cnt %>% 
  mutate(cmpr = if_else(cmpr == "pu_B01AB_plac", "pu_B01AB", cmpr),
         across(c(Estimate, `l-95% CI`, `u-95% CI`), function(x) {
    exp(x) %>% 
      round(2) %>% 
      formatC(digits = 2, format = "f")
  }))
model_n_lkp <- c(
  m0 = "a_arm",
  m1 = "b_como_cnt",
  m3 = "c_arm_como_cnt_unad",
  m8 = "d_arm_como_cnt_adj")
non_cnt2 <- non_cnt %>%
  mutate(res = paste0(Estimate, " (",`l-95% CI`, "-",`u-95% CI`, ")"),
         model_n2 = model_n_lkp[paste0("m", model_n)]) %>% 
  select(cond_grp, outcome, condition, intervention = medicine, comparator = control, cmpr, model_n2, n_trials,  res) %>% 
  spread(model_n2, res) %>% 
  separate(cmpr, into = c("atc1", "atc2"), sep = "_", extra = "merge") %>% 
  select(cond_grp, outcome, condition, intervention, atc2, #intervention_atc = atc2, 
          comparator, 
         atc1,
         n_trials, a_arm, b_como_cnt, c_arm_como_cnt_unad, d_arm_como_cnt_adj) %>% 
  arrange(cond_grp)

non_cnt3 <- non_cnt2 %>% 
  mutate(comparator = case_when(
    comparator %in% c("plac", "placebo") ~ "Placebo",
    comparator == "warf" ~ "Warfarin",
    comparator == "lmwh" ~ "LMWH"),
    intervention = str_to_sentence(intervention),
    intervention = paste0(intervention, " (", atc2, ")"),
    comparator = if_else(comparator == "Placebo", "Placebo", paste0(comparator, " (", atc1, ")"))) %>% 
  select(condition, outcome, intervention, comparator, n_trials,
         a_arm, b_como_cnt, c_arm_como_cnt_unad, d_arm_como_cnt_adj) %>% 
  arrange(condition, outcome, intervention, comparator)

non_cnt4 <- non_cnt3 %>% 
  select(-a_arm, -b_como_cnt, -c_arm_como_cnt_unad)
names(non_cnt4) <- str_to_sentence(names(non_cnt4))
non_cnt4 <- non_cnt4 %>% 
  rename(`Comorbidity treatment interaction` = D_arm_como_cnt_adj) %>% 
  group_by(Condition) %>% 
  mutate(Outcome = if_else(duplicated(Outcome), "", Outcome)) %>% 
  ungroup() %>% 
  mutate(Condition = if_else(duplicated(Condition), "", Condition),
         Condition = if_else(Condition == "treat", "Treatment", str_to_sentence(Condition)),
         Outcome = str_replace_all(Outcome, "dvt", "DVT") %>% 
           str_replace_all("pe", "PE") %>% 
           str_replace_all("VERT", "vertebrae") %>% 
           str_replace_all("^no\\.", "No. of ") %>% 
           str_replace("^b", "B")) 

## Note for semuloparin model did not converge in original trial IPD so not included (SE >2000)
write_csv(non_cnt3, "Outputs/res_non_cnt_out.csv")
write_csv(non_cnt4, "Outputs/table3.csv")

## random effects for appendix ----
cnt_re <- readRDS("Scratch_data/mdls_cnt_re.Rds")
cnt_re <- cnt_re %>% 
  select(model_type, condition_cb, cmpr, data) %>% 
  unnest(data) %>% 
  distinct(model_type, condition_cb, cmpr, nct_id) %>% 
  count(condition_cb, cmpr) %>% 
  filter(n %in% 2:4)

t2_sa <- bth_smry %>% 
  filter(model_grp %in% c("cnt_or", "cnt_re")) %>% 
  select(model_grp, condition_cb, cmpr, res)  %>% 
  inner_join(cnt_re) %>% 
  spread(model_grp, res) %>% 
  rename("Fixed effects" = cnt_or, "Random effects" = cnt_re,
         `N trials` = n)
t2_sa <- Relbl(t2_sa)
write_csv(t2_sa, "Outputs/Comorbidity_count_re.csv")
