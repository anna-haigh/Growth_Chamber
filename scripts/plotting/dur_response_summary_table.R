library(tidyverse)
library(lubridate)
library(dplyr)
library(knitr)
library(kableExtra)
library(emmeans)
library(nlme)
library(flextable)
library(officer)

# Load data
emmdcrit      <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmdcrit.csv")
emmd15        <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmd15.csv")
emmd50        <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmd50.csv")
dd1resilemm   <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_d1resilemm.csv")
dresisemm     <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_resisemm.csv")
demm_rec_rate <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_emm_rec_rate.csv")
demm_mod_day  <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_emm_mod_day.csv")

head(demm_rec_rate)
# Function for big table
tidy_emm <- function(summary_df, metric_name, group_vars = c("Elev", "Prev_trmt", "Needle_age")) {
  df <- as.data.frame(summary_df)
  
  if ("emmean" %in% names(df)) df <- df %>% rename(estimate = emmean)
  if (!"p.value" %in% names(df)) df$p.value <- NA_real_
  
  df %>%
    dplyr::select(dplyr::all_of(group_vars), estimate, SE, lower.CL, upper.CL, p.value) %>%
    dplyr::mutate(Metric = metric_name, .before = 1)
}


dcrit_tab       <- tidy_emm(emmdcrit, "Dcrit")
d50_tab         <- tidy_emm(emmd50, "D50")
d15_tab         <- tidy_emm(emmd15, "D15")
dresilience_tab  <- tidy_emm(dd1resilemm, "Resilience")
dresistance_tab  <- tidy_emm(dresisemm, "Resistance")
drate_tab        <- tidy_emm(demm_rec_rate, "Recovery rate coefficient")
drecdate_tab     <- tidy_emm(demm_mod_day, "Recovery date")


all_metricsd <- bind_rows(
  dcrit_tab, d15_tab, d50_tab, dresilience_tab, dresistance_tab, drate_tab, drecdate_tab
)

head(all_metricsd)

table_wided <- all_metricsd %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control", "mod", "high"), labels = c("control","moderate","high"))) %>%
  mutate(
    margin = upper.CL - estimate,
    estimate_fmt = sprintf("%.2f ± %.2f", estimate, margin)
  ) %>%
  select(Prev_trmt, Elev, Needle_age, Metric, estimate_fmt) %>%
  pivot_wider(
    names_from = Metric,
    values_from = estimate_fmt
  ) %>%
  arrange(Prev_trmt, Elev, Needle_age)

table_wided

table_wided <- table_wided %>%
  rename(
    `Previous Treatment` = Prev_trmt,
    Elevation    = Elev,
    `Needle Age` = Needle_age,
  )

names(table_wided)

ft2 <- flextable(table_wided) %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(part = "header") %>%
  merge_v(j = c("Needle Age", "Elevation", "Previous Treatment")) %>%
  valign(valign = "top", part = "body") %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ft2 

doc2 <- read_docx() %>%
  body_add_flextable(ft2) %>%
  body_end_section_landscape()

print(doc2, target = "~/Masters/Growth_Chamber/results/tables/duration_recovery_table.docx")
