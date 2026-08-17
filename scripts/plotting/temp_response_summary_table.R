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
emmtcrit     <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmtcrit.csv")
emmt15       <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmt15.csv")
emmt50       <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmt50.csv")
d1resilemm   <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_d1resilemm.csv")
resisemm     <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_resisemm.csv")
emm_rec_rate <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_emm_rec_rate.csv")
emm_mod_day  <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_emm_mod_day.csv")

head(emmtcrit)
# Function for big table
tidy_emm <- function(summary_df, metric_name, group_vars = c("Elev", "Micro", "Needle_age")) {
  df <- as.data.frame(summary_df)
  
  if ("emmean" %in% names(df)) df <- df %>% rename(estimate = emmean)
  if (!"p.value" %in% names(df)) df$p.value <- NA_real_
  
  df %>%
    dplyr::select(dplyr::all_of(group_vars), estimate, SE, lower.CL, upper.CL, p.value) %>%
    dplyr::mutate(Metric = metric_name, .before = 1)
}


tcrit_tab       <- tidy_emm(emmtcrit, "Tcrit")
t50_tab         <- tidy_emm(emmt50, "T50")
t15_tab         <- tidy_emm(emmt15, "T15")
resilience_tab  <- tidy_emm(d1resilemm, "Resilience")
resistance_tab  <- tidy_emm(resisemm, "Resistance")
rate_tab        <- tidy_emm(emm_rec_rate, "Recovery rate coefficient")
recdate_tab     <- tidy_emm(emm_mod_day, "Recovery date")


all_metrics <- bind_rows(
  tcrit_tab, t15_tab, t50_tab, resilience_tab, resistance_tab, rate_tab, recdate_tab
)

head(all_metrics)


table_wide <- all_metrics %>%
  mutate(
    margin = upper.CL - estimate,
    estimate_fmt = sprintf("%.2f ± %.2f", estimate, margin)
  ) %>%
  dplyr::select(Elev, Micro, Needle_age, Metric, estimate_fmt) %>%
  tidyr::pivot_wider(
    names_from = Metric,
    values_from = estimate_fmt
  ) %>%
  dplyr::arrange(Elev, Micro, Needle_age)

table_wide

table_pub <- table_wide %>%
  dplyr::rename(
    Elevation    = Elev,
    Microbiome   = Micro,
    `Needle Age` = Needle_age,
  )

table_pub <- table_pub %>%
  dplyr::rename(
    `Tcrit (°C)`        = Tcrit_estimate_fmt,
    `Tcrit p`            = Tcrit_p_fmt,
    `T50 (°C)`           = T50_estimate_fmt,
    `T50 p`              = T50_p_fmt,
    `T85 (°C)`           = T85_estimate_fmt,
    `T85 p`              = T85_p_fmt,
    `Resilience`         = Resilience_estimate_fmt,
    `Resilience p`       = Resilience_p_fmt,
    `Resistance`         = Resistance_estimate_fmt,
    `Resistance p`       = Resistance_p_fmt,
    `Recovery rate`      = `Recovery rate_estimate_fmt`,
    `Recovery rate p`    = `Recovery rate_p_fmt`,
    `Recovery date (days)` = `Recovery date_estimate_fmt`,
    `Recovery date p`      = `Recovery date_p_fmt`
  )


names(table_wide)

ft <- flextable(table_pub) %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(part = "header") %>%
  merge_v(j = c("Needle Age", "Elevation")) %>%
  valign(valign = "top", part = "body") %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()

ft

doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_end_section_landscape()

print(doc, target = "~/Masters/Growth_Chamber/results/tables/thermal_recovery_table.docx")
