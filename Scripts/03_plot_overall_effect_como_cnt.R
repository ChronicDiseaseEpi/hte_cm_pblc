library(tidyverse)
library(brms)
library(ggthemes)

## read in treatmetn effect data data ---- 
bth <- readRDS(file = "Scratch_data/bothsexdata.Rds")
bth$trials <- map_int(bth$data, nrow)
bth_sngl <- bth %>% 
  filter(trials == 1)
bth_mlt <- bth %>% 
  filter(!trials == 1)
## Run models for overall (at standardised mean age) treatment effects ----
## Only run this if changing models, is time consuming
myfx <- function(){
  count <- 1
  ## Random effect
  bth_mlt$res <- map(bth_mlt$data, function(mydata) {
    print(count)
    count <<- count + 1
    brm(est_mn|se(se_mn) ~ 
          1 + (1|nct_id)
        , 
        data = mydata, 
        # prior = bprior,
        control = list(adapt_delta = 0.99),
        cores = 4, backend = "cmdstanr")  
  })
  saveRDS(bth_mlt, file = "Scratch_data/bothsexdatamdls_brm_re.Rds")
  count <- 1
  ## Fixed effect
  bth_mlt$res <- map(bth_mlt$data, function(mydata) {
    print(count)
    count <<- count + 1
    brm(est_mn|se(se_mn) ~ 
          1
        , 
        data = mydata, 
        # prior = bprior,
        control = list(adapt_delta = 0.999),
        cores = 4, backend = "cmdstanr")  
  })
  saveRDS(bth_mlt, file = "Scratch_data/bothsexdatamdls_brm_fe.Rds")
}
# myfx()
bth_brm_re <- readRDS(file = "Scratch_data/bothsexdatamdls_brm_re.Rds")
bth_brm_fe <- readRDS(file = "Scratch_data/bothsexdatamdls_brm_fe.Rds")

ExtractRename <- function(mdls, suffix = "re"){
  mdls$fixed <- map(mdls$res, ~ (summary(.x)$fixed) %>% as_tibble(rownames = "term"))
  mdls <- mdls %>% 
    unnest(fixed)
  mdls <- mdls %>% 
    mutate(term = str_remove(term, "^term")) %>% 
    select(-Rhat, -Bulk_ESS, -Tail_ESS)
  mdls <- mdls %>% 
    rename(est = Estimate, se = Est.Error, lci = `l-95% CI`, uci = `u-95% CI`)
  mdls <- mdls %>% 
    select(outcome, condition_cb, cmpr, trials, est, se, lci, uci)
  mdls
}
bth_brm_re <- ExtractRename(bth_brm_re)
bth_brm_fe <- ExtractRename(bth_brm_fe, suffix = "fe")

## Choose based on trial N, cut off is 5
bth_brm <- bind_rows(bth_brm_fe %>% filter(trials <  5),
                     bth_brm_re %>% filter(trials >= 5)) %>% 
  select(-outcome)
bth_all <- bind_rows(
  raw = bth %>% 
    select(condition_cb, cmpr, data, trials) %>% 
    unnest(data) %>% 
    rename(est = est_mn, se = se_mn) %>% 
    mutate(lci = est - 1.96*se,
           uci = est + 1.96*se),
  meta = bth_brm %>% 
    mutate(nct_id = " Meta-analysis"),
  .id = "res_type") %>% 
  arrange(condition_cb, cmpr, desc(nct_id))

## add in comorbidity count effect estimates from individual trials for plotting
cnt <- readRDS("Scratch_data/mdls_cnt.Rds")
cnt_raw <- cnt %>% 
  select(condition_cb, cmpr, data) %>% 
  unnest(data) %>% 
  select(condition_cb, cmpr, nct_id, est = estimate, se) %>% 
  mutate(lci = est - 1.96*se,
         uci = est + 1.96*se)
cnt$fixed <- map(cnt$res, ~ (summary(.x)$fixed) %>% as_tibble(rownames = "term"))

## add in comorbidity count meta-analyses for plotting
cnt_ma <- cnt %>% 
  select(condition_cb, cmpr, fixed) %>% 
  unnest(fixed) %>% 
  select(condition_cb, cmpr, est = Estimate, se = Est.Error,
         lci = `l-95% CI`,
         uci = `u-95% CI`) %>% 
  mutate(nct_id = " Meta-analysis")

## combine all trials and meta-analyses for comorbidity counts into a single table
cnt_bth <- bind_rows(
  raw = cnt_raw,
  meta = cnt_ma,
  .id = "res_type") %>% 
  arrange(condition_cb, cmpr, desc(nct_id))

## count the number of trials per condition/treatment comparison
cnt_bth <- cnt_bth %>% 
  group_by(condition_cb, cmpr) %>% 
  mutate(trials = sum(str_detect(nct_id, "NCT"))) %>% 
  ungroup()

## combine overall treatment effects and comorbidity count effects into a single table
bth_all <- bind_rows(`Treatment` = bth_all,
                     `Treatment interaction` = cnt_bth,
                     .id = "result")

## Assign colours shape and alpha
bth_all2 <- bth_all %>% 
  mutate(mycolour = case_when(
           result == "Treatment" ~ "black",
           result == "Treatment interaction" ~ "red"),
         myshape = if_else(res_type == "raw", 1, 18),
         myalpha = if_else(res_type == "raw", 0.2, 1))

bth_all3 <- bth_all2 %>% 
  mutate(nct_id = case_when(
    !res_type == "raw" & result == "Treatment" ~ " Meta-analysis",
    !res_type == "raw" & result == "Treatment interaction" ~ " Meta-analysis",
    TRUE ~ nct_id))
bth_all <- bth_all3
rm(bth_all2, bth_all3)
rm(bth, bth_brm, bth_brm_fe, bth_brm_re, bth_mlt, bth_sngl, cnt, cnt_bth, cnt_ma, cnt_raw)

## convert ATC codes to upper case etc
bth_all <- bth_all %>% 
  separate(cmpr, into = c("one", "two"), sep = "_", extra = "merge") %>% 
  mutate(one = str_to_upper(one),
         two = str_to_upper(two),
         cmpr = if_else(one == "PU", two, paste0(one, " vs ", two)),
         cmpr = str_replace_all(cmpr, "_", "-")) %>% 
  select(-one, -two)
write_csv(bth_all, "Outputs/main_eff_nter_eff.csv")

## Plot code for meta-analysis ----
## relabel facets
CleanDrugCode <- function(x) {
  y <- x %>% 
    str_split_fixed(pattern = "\\(", n = 2)
  y1 <- y[,1, drop = TRUE]
  y1 <- str_trim(y1)
  y1 <- case_when(
    y1 == "Ank spond" ~ "Ankylosing spondylitis",
    y1 == "ibd" ~ "IBD",
    y1 == "Systemic Lupus Erythematosus" ~ "SLE",
    y1 == "inflammatory arthropathy" ~ "Inflammatory arthropathy",
    y1 == "Pulmonary Disease, Chronic Obstructive" ~ "COPD",
    y1 == "Parkinson Disease" ~ "Parkinson's disease",
    TRUE ~ y1)
  
  y2 <- y[,2, drop = TRUE]
  y2 <- y2 %>%  
    str_replace("pu_", "") %>%
    str_replace_all("_il", "-il") %>%
    str_replace("_", " VS ") %>%
    str_to_upper()
  paste0(y1, " (", y2)
}

## Add drug class detailed labels (so same as in tables)
dc_long <- read_csv("Outputs/lookupdrugclassesdetailed.csv")
## joining on NCT ID then propagating within comparison
bth_all <- bth_all %>% 
  left_join(dc_long) %>% 
  group_by(cmpr) %>% 
  mutate(dc_detailed = if_else(is.na(dc_detailed), dc_detailed[1], dc_detailed)) %>% 
  ungroup()
## change insulin and analogues to sentence case (although in CAPS in ATC)
bth_all <- bth_all %>% 
  mutate(dc_detailed = str_replace(dc_detailed, "INSULINS AND ANALOGUES", "Insulins and analogues")) 

bth_all <- bth_all %>% 
  mutate(condition_cb2 =  case_when(
    condition_cb == "Ank spond" ~ "Ankylosing spondylitis",
    condition_cb == "ibd" ~ "IBD",
    condition_cb == "Systemic Lupus Erythematosus" ~ "SLE",
    condition_cb == "inflammatory arthropathy" ~ "Inflammatory arthropathy",
    condition_cb == "Pulmonary Disease, Chronic Obstructive" ~ "COPD",
    condition_cb == "Parkinson Disease" ~ "Parkinson's disease",
    TRUE ~ condition_cb))

## Order alphabetically to allow splitting at 50%
alphab <- bth_all %>% 
  filter(!trials == 1) %>% 
  arrange(condition_cb2, dc_detailed) %>% 
  group_by(condition_cb2, dc_detailed) %>% 
  summarise(n_trials = sum(!duplicated(nct_id[str_detect(nct_id, "^NCT")]))) %>% 
  ungroup() %>% 
  mutate(m = n_trials/sum(n_trials),
         cumm = cumsum(m)) 
alphab$condition_cb2[alphab$cumm < 0.5]  %>% unique() %>% dput()
splitma <- c("Ankylosing spondylitis", "Asthma", "BPH", "CIU", "COPD", "Dementia", 
             "Diabetes", "GORD")
bth_all <- bth_all %>% 
  mutate(splitma = if_else(condition_cb2 %in% splitma, 1L, 2L))
# Make facet labelsplit DC detailed
LabelCondDDc <- function(x) {
  x %>% 
    mutate(faclabel = paste0(
      condition_cb2,
      "\n",
      str_replace(dc_detailed, " vs ", " vs\n")))
}
bth_all <- LabelCondDDc(bth_all )

## Alternative plots as per reviewer comments
bth_nst <- bth_all %>% 
  group_by(condition_cb2, dc_detailed) %>% 
  nest() %>% 
  ungroup()
bth_nst$plot <- pmap(list(bth_nst$data, bth_nst$condition_cb2, bth_nst$dc_detailed), function(.x, .y, .z) {
  ggplot(.x, aes(x = nct_id, 
                 y = est, 
                 ymin = lci, 
                 ymax = uci, 
                 shape = myshape,
                 colour = mycolour)) +
    geom_point() +
    geom_linerange() +
    coord_flip(ylim = c(-1.5, 1.5)) +
    geom_hline(yintercept = 0, colour = "grey") +
    scale_colour_identity() +
    scale_shape_identity() +
    scale_x_discrete("") +
    scale_y_continuous("") +
    facet_wrap(~mycolour) +
    ggtitle(.y, subtitle = .z) + 
    ggthemes::theme_tufte() +
    theme(strip.text.x = element_blank(),
          panel.spacing.x = unit(-5, "lines"))
})
bth_nst$plot[[1]]
bth_nst$n_trials <- map_int(bth_nst$data, ~ sum(!duplicated(.x$nct_id[str_detect(.x$nct_id, "^NCT")])))
bth_meta <- bth_nst %>% 
  filter(!n_trials == 1) %>% 
  mutate(grps = c(rep(letters[1:4], each = 6), "e", "e")) %>% 
  select(grps, plot) %>% 
  nest(data = plot)
bth_meta$cp <- map(bth_meta$data, ~ cowplot::plot_grid(plotlist = .x$plot, ncol = 2))
map2(bth_meta$grps, bth_meta$cp, ~ {
  tiff(filename = paste0("Outputs/Figure2", .x, ".tiff"), 
       compression = "lzw", res = 300,
       height = 12,
       width = 12,
       unit = "in")
  print(.y)
  dev.off()
  })

## plot single trials
sngl <- bth_all %>% 
  filter(trials == 1, !res_type == "meta", str_detect(nct_id, "^NCT")) %>% 
  mutate(mycolour = case_when(
    result == "Treatment" ~ "black",
    result == "Treatment interaction" ~ "red"))

sngl_plot <- ggplot(sngl, aes(x = faclabel, y = est, ymin = lci, ymax = uci,
                              colour = mycolour)) +
  geom_point(position = position_dodge(-0.5), shape = 1) +
  geom_linerange(position = position_dodge(-0.5)) +
  coord_flip(ylim = c(-3.5, 3.5)) +
  scale_colour_identity() +
  ggthemes::theme_tufte() +
  geom_hline(yintercept = 0, alpha = 0.7) +
  scale_x_discrete("", limits = rev) +
  scale_y_continuous("") +
  theme(axis.ticks.y = element_blank())
sngl_plot

tiff("Outputs/Figure3.tiff", res = 600, height = 20, width = 15, unit = "cm", compression = "lzw")
sngl_plot
dev.off()

a <- magick::image_read("Outputs/Figure3.tiff", density = 600)
a <- magick::image_annotate(a, "Better efficacy", size = 50, gravity = "southwest", location = "+2200+90", style = "italic")
a <- magick::image_annotate(a, "Worse efficacy",  size = 50, gravity = "southwest", location = "+3000+90", style = "italic")
magick::image_write(a, path = "Outputs/Figure3_lbl.tiff", format = "tiff")
