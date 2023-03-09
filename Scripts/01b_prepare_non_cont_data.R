# prepare_non_cont_models
library(tidyverse)

## Metadata ----
csdr <- read_csv("Data/trial_metadata_csdr.csv")
yoda <- read_csv("Data/trial_metadata_yoda.csv")
como_names_csdr <- read_csv("Data/top6_csdr.csv")
como_names_yoda <- read_csv("Data/top6_yoda.csv")

## Headache trials ----
## 5 headache trials all have comorbidity count; these were Poisson models, 
## all Topiramate for migraine
hdch <- read_csv("Data/Model_coefficients_yoda_hdch.csv")
hdch <- yoda %>% 
  inner_join(hdch)
hdch <- hdch %>% 
  mutate(term = if_else(term == "como_cnt:armtreat", "armtreat:como_cnt", term))
hdch %>% 
  filter(model_n %in% c(1, 3, 8), main_exposure == "cc") %>% 
  count(model_n, nct_id, term) %>% 
  spread(term, n, fill = "-")
hdch <- hdch %>% 
  filter((model_n == 1 & term == "como_cnt") |
           (model_n == 3 & term %in% c("armtreat","armtreat:como_cnt")) |
           (model_n == 8 & term == "armtreat:como_cnt")) %>% 
  mutate(model_n = if_else(model_n == 3 & term == c("armtreat"), 0, model_n)) %>% 
  arrange(model_n)

## Cox models, various ----
cox <- read_csv("Data/Model_coefficients_cox.csv.gz")
cox <- cox %>% 
  filter(outcome_description %in% c("coxph_models_cc_update")) %>% 
  mutate(term = case_when(
    term == "como_cnt:arm" ~ "arm:como_cnt",
    term == "sex:arm" ~ "arm:sex",
    TRUE ~ term))
## must have comorbidity count (ie comorbidity data). As expected this drop 4 trials 
## (all cardiovascular which did not record comorbidities) and 11 outcomes
cox <- cox %>% 
  group_by(nct_id) %>% 
  mutate(anycomo = any(term == "como_cnt")) %>% 
  ungroup() 
cox <- cox %>% 
  filter(anycomo) %>% 
  select(-anycomo)
cox <- csdr %>% 
  inner_join(cox)

## examine binary
cox_cb <- read_csv("Data/Model_coefficients_cox.csv.gz") %>% 
  filter(outcome_description %in% c("coxph_models_cb_update")) %>% 
  semi_join(cox %>% distinct(nct_id))
  
## limit to relevant models. Note pulling model with comorbidity interaction 
## (but not age and sex twice, once to get treatment effect in model with arm:como interaction)
cox <- cox %>% 
  filter((model_n == 1 & term == "como_cnt") |
         (model_n == 3 & term %in% c("arm","arm:como_cnt")) |
         (model_n == 8 & term == "arm:como_cnt")) %>% 
  mutate(model_n = if_else(model_n == 3 & term == c("arm"), 0, model_n))

## Osteoporosis trials (subset of Cox) ----
osteop <- cox %>% 
  filter(nct_id %in%  c("NCT00670501", "NCT00049829", "NCT00046254"),
         outcome %in% c("Fracture of VERT", "new vertebral fractures")) %>% 
  arrange(model_n)
osteop_cb <- cox_cb %>% 
  filter(nct_id %in% unique(osteop$nct_id))
osteop_cb %>% 
  anti_join(como_names_csdr %>% filter(condition_grp == "Osteoporosis"))
osteop_cb <- osteop_cb %>% 
  mutate(term = str_remove(term, "TRUE$"))
osteop_lk <- como_names_csdr %>% filter(condition_grp == "Osteoporosis") %>% 
  select(v1:v6) %>% 
  distinct() %>% 
  gather("term", "new_term")
osteop_lk2 <- osteop_lk %>% 
  mutate(term = paste0("arm:", term))
osteop_lk <- bind_rows(osteop_lk,
                       osteop_lk2)
rm(osteop_lk2)
osteop_cb <- osteop_cb %>% 
  left_join(osteop_lk) %>% 
  mutate(term = if_else(!is.na(new_term), new_term, term)) %>% 
  select(-new_term)
# drop from remaining 
cox <- cox %>% 
  filter(!nct_id %in% osteop$nct_id)
# drop one diabetes and one erosive esoph trial as cannot fit meta-analysis to these
other <- cox %>% 
  filter(nct_id %in% c("NCT00968708", "NCT00321737"))
## Thromboembolic or bleeding (subset of Cox) ----
tbe_bleed <- cox %>% 
  filter(!nct_id %in% c("NCT00968708", "NCT00321737"))
rm(cox, other)
tbe_bleed_tx <- read_csv("Supporting/tbe_drugs.csv")
tbe_bleed_cond <- read_csv("Supporting/tbe_condition.csv")
tbe_bleed <- tbe_bleed %>% 
  select(nct_id, term, estimate, std.error, outcome, model_n) %>% 
  left_join(tbe_bleed_tx) %>% 
  left_join(tbe_bleed_cond)
rm(tbe_bleed_tx, tbe_bleed_cond)
tbe_bleed_cb <- cox_cb %>% 
  filter(nct_id %in% unique(tbe_bleed$nct_id))
## All in data, all except diabetes trial have same conditions
como_names_csdr %>% 
  filter(nct_id %in% tbe_bleed_cb$nct_id)
## Bleeding (subset of tbe_bleed trials) ---- 
## 8 of 9 trials have a bleeding measure, the one post op (of two post op ones) does not NCT00168818
bleed <- tbe_bleed %>% 
  filter(!nct_id %in% "NCT00168818",
         outcome %in% c("bleeding", "CLINICALLY RELEVANT NON MAJOR BLEEDING", "Inv reported-All bleeding events"))  %>% 
  arrange(model_n)

## DVT or PE (subset of tbe trials) ----
## 8 of 9 trials have dvt or pe, choose most similar outcomes across trials. 
## Trial without DVT or PE is a stroke trial
`US DVT left leg` <- c("NCT00152971", "NCT00168818")
`DVT or PE` <- c("NCT00291330", "NCT00680186", "NCT00694382", "NCT00657150")
DVT <- c("NCT00329238", "NCT00558259")
tbe <- tbe_bleed %>% 
  filter(!nct_id %in% "NCT00262600") %>% 
  filter( (nct_id %in% DVT & outcome == "DVT") |
          (nct_id %in% `DVT or PE` & outcome == "DVT or PE") |
          (nct_id %in% `US DVT left leg` & outcome == "US DVT left leg")) %>% 
  arrange(model_n)

rm(`DVT or PE`, `US DVT left leg`,
   como_names_csdr, como_names_yoda, cox_cb, csdr, DVT, 
   osteop_cb, osteop_lk, tbe_bleed, tbe_bleed_cb, yoda)

## Save outputs and create summaries ----
savlst <- list(tbe_bleed = bleed, tbe_dvtpe = tbe, hdch = hdch, osteop = osteop)
saveRDS(savlst, 
        "Scratch_data/cfs_non_cont.Rds")
tbe_out <-  savlst[c("tbe_dvtpe", "tbe_bleed")] %>%
  bind_rows(.id = "out_broad") %>% 
  distinct(nct_id, condition, cmpr, medicine, control, out_broad) %>% 
  group_by(nct_id, condition, cmpr, medicine, control) %>% 
  summarise(outcomes = paste(sort(out_broad), collapse = " & ")) %>% 
  ungroup()
hdch <- savlst$hdch %>% 
  distinct(nct_id, condition, outcome, medicine) %>% 
  mutate(cmpr = "N03AX11")
osteo_nct_id <- savlst$osteop %>% 
  distinct(nct_id) %>% 
  summarise(nct_ids = paste(nct_id %>% sort() %>% unique(), collapse = ";")) %>% 
  mutate(`Index conditions` = "Osteoporosis") 
osteo <- savlst$osteop %>% 
  distinct(nct_id, medicine, condition, outcome) %>% 
  mutate(cmpr = case_when(medicine == "Teriparatide" ~ "H05AA",
                          medicine == "zoledronic acid" ~ "M05BA"),
         outcome = "Vertebral fracture",
         condition = "Osteoporosis") %>% 
  count(condition, outcome, cmpr) %>% 
  mutate(trials = sum(n),
         n_tx_compar = sum(!duplicated(cmpr))) %>% 
  group_by(condition, outcome, trials, n_tx_compar) %>% 
  summarise(cmpr = paste(cmpr, n) %>% paste(collapse = "; ")) %>% 
  rename(`Index conditions` = condition, Outcome = outcome, Trials = trials, `N Tx compar` = n_tx_compar, `Tx compar` = cmpr) %>% 
  mutate(Nesting = "simple") %>% 
  ungroup()
osteo <- osteo %>% 
  inner_join(osteo_nct_id)
tbe_out_nct_id <- tbe_out %>% 
  distinct(nct_id) %>% 
  summarise(nct_ids = paste(nct_id %>% sort() %>% unique(), collapse = ";")) %>% 
  mutate(`Index conditions` = "Thromboembolic") 
tbe_out2 <- tbe_out %>% 
  transmute(nct_id = nct_id,
            condition = "Thromboembolic",
         cmpr = cmpr %>% str_remove("pu_") %>% 
           str_remove("_plac") %>% 
           str_replace("B01AE07", "B01AE") %>% 
           str_replace("B01AA03", "B01AA") %>% 
           str_replace("_", " vs "),
         outcomes = case_when(
           outcomes == "tbe_dvtpe" ~ "DVT or PE",
           outcomes == "tbe_bleed" ~ "Bleeding",
           TRUE ~ "DVT or PE and Bleeding")) %>% 
  group_by(condition, cmpr) %>% 
  mutate(cmpr_n = sum(!duplicated(nct_id))) %>% 
  group_by(condition, outcomes) %>% 
  mutate(out_n = sum(!duplicated(nct_id))) %>% 
  ungroup() %>%
  mutate(trials = sum(!duplicated(nct_id))) %>% 
  select(-nct_id) %>% 
  distinct() %>% 
  mutate(n_cmprs = sum(!duplicated(cmpr)),
         cmpr = paste(cmpr, cmpr_n),
         outcomes = paste(outcomes, out_n)) %>% 
  distinct(condition, trials, n_cmprs, cmpr, outcomes) %>% 
  group_by(condition, trials, n_cmprs) %>% 
  summarise(across(c(cmpr, outcomes), ~ .x %>% sort() %>% unique() %>% paste(collapse = "; "))) %>% 
  ungroup()
tbe_out2 <- tbe_out2 %>% 
  rename(`Index conditions` = condition, Outcome = outcomes, Trials = trials, `N Tx compar` = n_cmprs, `Tx compar` = cmpr) %>% 
  mutate(Nesting = "complex")
tbe_out_conditions <- tbe_out %>% 
  transmute(condition = condition %>% 
  str_replace("treat", "Treatment") %>% 
  str_to_title()) %>% 
  count(condition) %>% 
  mutate(condition = paste(condition, n, " trials")) %>% 
  pull(condition) %>% 
  paste(collapse = "; ")
tbe_out2 <- tbe_out2 %>% 
  inner_join(tbe_out_nct_id)
hdch_nct_id <- hdch %>% 
  distinct(nct_id) %>% 
  summarise(nct_ids = paste(nct_id %>% sort() %>% unique(), collapse = ";")) %>% 
  mutate(`Index conditions` = "Migraine") 
hdch2 <- hdch %>% 
  mutate(outcome = str_replace(outcome, "no\\.", "No. ")) %>% 
  count(condition, outcome, cmpr) %>% 
  mutate(cmpr = paste(cmpr, n),
         Trials = n) %>% 
  rename(`Index conditions` = condition, Outcome = outcome, Trials = Trials, `N Tx compar` = n, `Tx compar` = cmpr) %>% 
  mutate(Nesting = "simple")
hdch2 <- hdch2 %>% 
  inner_join(hdch_nct_id)
## Create categorical component of table 1
t1_catg <- bind_rows(tbe_out2 %>% slice(0),
                     hdch2,
                     osteo,
                     tbe_out2)

t1_catg2 <- t1_catg %>% 
  mutate(`Tx compar` = `Tx compar` %>% 
           str_replace("N03AX11", "topiramate (N03AX11)") %>% 
           str_replace("H05AA", "Parathyroid hormones and analogues (H05AA)") %>% 
           str_replace("M05BA", "Bisphosphonates (M05BA)") %>% 
           str_replace("B01AA", "Vitamin K antagonists (B01AA)") %>% 
           str_replace("B01AB", "Heparin group (B01AB)") %>% 
           str_replace("B01AE", "Direct thrombin inhibitors (B01AE)") %>% 
           str_replace_all("\\b([0-9]\\b)", "\\[\\1\\]") %>% 
           stringi::stri_replace_last_fixed(";", " &") %>% 
           str_replace_all(";", ","),
         Outcome =  str_replace_all(Outcome, "\\b([0-9]\\b)", "\\[\\1\\]")) %>% 
  select(`Index conditions`, Outcome, `Treatment comparisons` = `Tx compar`, Trials, nct_ids)

write_csv(t1_catg2, "Outputs/table1_catg.csv")
write_lines(tbe_out_conditions, "Outputs/note_on_tbe.txt")
