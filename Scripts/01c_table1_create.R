library(tidyverse)

## run code and read data ----
source("Scripts/01a_prepare_cont_data.R")
rm(oth, cnt)
t1_catg <- read_csv("Outputs/table1_catg.csv")
csdr_yoda <- read_csv("Data/age_sex_como_count.csv")

## Functions ----
## Can calculate combined standard deviation of groups as per
## https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formula_for_combining_groups.htm
CombMean <- function(n1, m1, n2, m2){
  top <- n1*m1 + n2*m2
  top/(n1 + n2)
}

CombSd <- function(n1, m1, n2, m2, s1, s2){
  top = (n1 - 1)*s1^2 + (n2 - 1)*s2^2 + (n1*n2/(n1 + n2)) * (m1^2 + m2^2 - 2*m1*m2)
  bottom = n1 + n2 - 1
  res = (top/bottom)^0.5
  res
}

CombSdVectorised <- function(n, m, s){
  myrows <- length(n)
  myrows_now <- myrows
  while (myrows_now >= 2) {
    ## select first two values
    n1 <- n[1]
    n2 <- n[2]
    s1 <- s[1]
    s2 <- s[2]
    m1 <- m[1]
    m2 <- m[2]
    ## replace first value with combination of first two values and drop second value
    new_s1 <- CombSd(n1, m1, n2, m2, s1, s2)
    new_m1 <- CombMean(n1, m1, n2, m2)
    new_n1 <- n1 + n2
    s[1] <- new_s1
    m[1] <- new_m1
    n[1] <- new_n1
    s <- s[-2]
    m <- m[-2]
    n <- n[-2]
    # print(s)
    # print(m)
    ## recalculate the length
    myrows_now <- length(s)
  }
  # print(sum(n))
  s
}

CollapseSex <- function(x) {
  tibble(
    age_m = weighted.mean(x$age_m, x$n),
    age_s = CombSdVectorised(x$n, x$age_m, x$ages_s),
    como_cnt = weighted.mean(x$como_cnt, x$n),
    male = sum(x$male*x$n),
    female = sum((1 - x$male)*x$n),
    tot = sum(x$n))}

## Summary statistics; Table 1 continuous ----
allcfs_smry <- allcfs %>% 
  distinct(model_type,
           nct_id,
           condition_cb,
           cmpr,
           outcome) 
## remove underline from il
allcfs_smry <- allcfs_smry %>%
  mutate(cmpr2 = str_replace_all(cmpr, "_il", "-IL") %>% 
           str_remove("pu_") %>% 
           str_to_upper() %>% 
           str_replace("_", " vs "))
tx_cmpar_lkp <- allcfs_smry %>% 
  distinct(cmpr, cmpr2)
allcfs_smry <- allcfs_smry %>% 
  select(-cmpr) %>% 
  rename(cmpr = cmpr2)

## convert ATC codes to names
atc <- read_csv("Supporting/relevant_atc_codes.csv")
allcfs_smry_lbl <- allcfs_smry %>% 
  group_by(model_type, condition_cb, cmpr) %>% 
  summarise(Outcome = paste(outcome %>% unique() %>% sort(), collapse = "; "),
            Trials = length(nct_id),
            nct_ids = paste(nct_id %>% unique() %>% sort(), collapse = ";")) %>% 
  ungroup() %>% 
  select(model_type, condition_cb, Outcome, cmpr,Trials, nct_ids)
allcodes <- allcfs_smry_lbl$cmpr %>% 
  unique()
allcodes <- allcodes %>% 
  str_split(pattern = "\\b") %>% 
  unlist() %>% 
  unique()
setdiff(allcodes, atc$code)
allcodes <- intersect(allcodes, atc$code) %>% sort()
allcodes <- atc %>% 
  filter(code %in% allcodes)
allcodes_lkp <- paste0(allcodes$nm, " (", allcodes$code, ")")
names(allcodes_lkp) <- allcodes$code
allcodes_lkp <- c(allcodes_lkp, "L01EX" = "Other protein kinase inhibitors (L01EX)")
allcfs_smry_lbl$cmpr2 <- allcfs_smry_lbl$cmpr
for (i in seq_along(allcodes_lkp)) {
  allcfs_smry_lbl$cmpr2 <- str_replace_all(allcfs_smry_lbl$cmpr2, names(allcodes_lkp)[i], allcodes_lkp[i]) %>% 
    print()
}
tx_cmpar_lkp <- tx_cmpar_lkp %>% 
  inner_join(allcfs_smry_lbl %>% 
               select(cmpr3 = cmpr2, cmpr2 = cmpr) %>% 
               distinct())
allcfs_smry_lbl <- allcfs_smry_lbl %>% 
  arrange(condition_cb)  %>% 
  mutate(condition_cb2 = case_when(
    condition_cb == "Ank spond" ~ "Ankylosing Spondylitis",
    condition_cb == "ibd" ~ "IBD",
    condition_cb == "inflammatory arthropathy" ~ "Inflammatory Arthropathy",
    condition_cb == "Psoriasisps_pa" ~ "Psoriasis",
    condition_cb == "Psoriasisra_pa" ~ "Psoriatic arthropathy",
    TRUE ~ condition_cb),
    cmpr2 = if_else(cmpr == "L04AB", "Tumor necrosis factor alpha (TNF-alpha) inhibitors (L04AB)", cmpr2))
allcfs_smry_lbl
trials_cmpr <- allcfs_smry %>% 
  distinct(nct_id, cmpr) 
trials_cmpr2 <- trials_cmpr %>% 
  inner_join(allcfs_smry_lbl %>% distinct(cmpr, cmpr2)) %>% 
  select(nct_id, dc_detailed = cmpr2)
trials_cmpr3 <- trials_cmpr %>% 
  inner_join(allcfs_smry_lbl %>% distinct(cmpr, cmpr2)) %>% 
  inner_join(allcfs_smry %>% distinct(nct_id, condition_cb)) %>% 
  inner_join(allcfs_smry_lbl %>% distinct(condition_cb, condition_cb2)) %>% 
  select(nct_id, cmpr = cmpr2, cond = condition_cb2) %>% 
  distinct()
write_csv(trials_cmpr3, "Outputs/lookup_tx_compar_trial_lvl.csv")
## create lookup for conditions and treatment comparisons
write_csv(tx_cmpar_lkp,
          "Outputs/lookup_tx_compar_lbl.csv")
write_csv(allcfs_smry_lbl %>% 
            distinct(condition_cb, condition_cb2),
          "Outputs/lookup_cond_lbl.csv")
write_csv(trials_cmpr2,
          "Outputs/lookupdrugclassesdetailed.csv")
## Create table 1
t1_cont <- allcfs_smry_lbl %>% 
  select(-cmpr, -condition_cb) %>% 
  rename(`Index conditions` = condition_cb2, `Treatment comparisons` = cmpr2) %>% 
  filter(!model_type == "intrvl") %>% 
  select(-model_type) %>%
  distinct() %>% 
  arrange(`Index conditions`, `Treatment comparisons`) 
t1_cont <- t1_cont %>% 
  mutate(`Treatment comparisons` = paste0(`Treatment comparisons`, " [", Trials, "]")) %>% 
  group_by(`Index conditions`, Outcome) %>% 
  summarise(`Treatment comparisons` = paste(`Treatment comparisons`, collapse = ", "),
            Trials = sum(Trials),
            nct_ids = paste(nct_ids, collapse = ";")) %>% 
  ungroup()
t1_cont <- t1_cont %>% 
  mutate(`Treatment comparisons` = stringi::stri_replace_last_fixed(`Treatment comparisons`, ',', ' &'))
rm(allcfs_smry, allcfs_smry_lbl, allcodes, atc, 
   trials_cmpr, trials_cmpr2, trials_cmpr3,
   tx_cmpar_lkp, i, allcodes_lkp)

## read categorical table 1  ----
t1w <- bind_rows(cont = t1_cont,
                 catg = t1_catg,
                 .id = "outcome_type")
t1w <- t1w %>% 
  mutate(t1id = seq_along(outcome_type))
t1 <- t1w %>% 
  select(t1id, outcome_type, nct_ids) 
t1$nct_id <- str_split(t1w$nct_ids, pattern = ";")
t1 <- t1 %>% 
  select(-nct_ids) %>% 
  unnest(nct_id) %>% 
  distinct()
rm(t1_catg, t1_cont)

## Calculate total in each
csdr_yoda %>% 
  semi_join(t1 %>% 
              filter(outcome_type == "cont")) %>% 
  summarise(trials = length(.data$nct_id),
            participants = sum(tot))
csdr_yoda %>% 
  semi_join(t1 %>% 
              filter(outcome_type == "catg")) %>% 
  summarise(trials = length(.data$nct_id),
            participants = sum(tot))
csdr_yoda %>% 
  semi_join(t1 %>% 
              filter(outcome_type == "catg")) %>% 
  semi_join(t1 %>% 
              filter(outcome_type == "cont")) %>% 
  summarise(trials = length(.data$nct_id),
            participants = sum(tot))
csdr_yoda %>% 
  semi_join(t1) %>% 
  summarise(trials = length(.data$nct_id),
            participants = sum(tot))

## aggregate by outcome-type/condition/treatment
t1_agg <- t1 %>% 
  left_join(csdr_yoda) 
t1_agg <- t1_agg %>% 
  group_by(t1id) %>% 
  nest() %>% 
  ungroup()
# summarise means within these
t1_agg$data2 <- map(t1_agg$data, ~ .x %>% 
                      summarise(age_s = CombSdVectorised(tot, age_m, age_s),
                                age_m = weighted.mean(age_m, tot),
                                male = sum(male),
                                female = sum(female),
                                tot = sum(tot),
                                como_cnt = weighted.mean(como_cnt),
                                nct_id = paste(nct_id %>% sort() %>% unique(), collapse = ";")) %>% 
                      mutate(male = paste0(round(100*male/tot, 1), "%"),
                             como_cnt = round(como_cnt, 2),
                             age = paste0(round(age_m, 1), " (", round(age_s, 1), ")")) %>% 
                      select(nct_id, participants = tot, age, male, como_cnt))
t1_agg_smry <- t1_agg %>% 
  select(t1id, data2) %>% 
  unnest(data2)
t1w2 <- t1w %>% 
  left_join(t1_agg_smry) %>% 
  select(-nct_ids, -t1id)
## calculate % with each comorbidity count 
t1w2$res <- map(t1w2$como_cnt, ~ round(100*dpois(0:2, .x),1) %>% as.list())
t1w2$res <- map(t1w2$res, ~ {
  names(.x) <- c("Zero", "One", "Two")
  as_tibble(.x) %>% 
    mutate("Three+" = round(100 - (Zero + One + Two),1))})
t1w2 <- t1w2 %>% 
  unnest(res)
t1w2 <- t1w2 %>% 
  mutate(across(c( Zero, One, Two, `Three+`), ~ if_else(is.na(.x), "", paste0(.x, "%"))))
write_csv(t1w2, "Outputs/table1_with_age_sex_como.csv")

## breakdown
t1w2_per_out <- t1w2 
t1w2_per_out$nct_id <- str_split(t1w2_per_out$nct_id, patt = ";")
t1w2_per_out <- t1w2_per_out %>% 
  unnest(nct_id)
t1w2_per_out <- t1w2_per_out %>% 
  distinct(outcome_type, nct_id, participants)
t1w2_per_out %>% 
  group_by(outcome_type) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            participants = sum(participants[!duplicated(nct_id)])) %>% 
  ungroup()

a <- t1w2_per_out %>% 
  distinct(nct_id, participants)

# %>% 
#   summarise(trials = sum(!duplicated(nct_id)),
#             participants = sum(participants[!duplicated(nct_id)])) 
