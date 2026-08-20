library(tidyverse)
library(lubridate)
library(dplyr)

dur <- read.csv("~/Masters/Growth_chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_treatment2_full.csv")
head(dur)

dur<- dur %>%
  mutate(.,
         Seedling_ID = factor(Seedling_ID),
         Elev = factor(Elev),
         Genotype = factor(Genotype),
         Micro = factor(Micro),
         Needle_age = factor(Needle_age),
         Measurement = as.factor(Measurement),
         Date = mdy(Date),
         Fo = as.numeric(Fo),
         Fv = as.numeric(Fv),
         Fm = as.numeric(Fm),
         Fv.Fm=as.numeric(Fv.Fm))

dur$Seed_ID <- paste(dur$Seedling_ID, dur$ID2)

head(dur)

# Check if mean initial Fv/Fm differs between groups
initial <- dur %>%
  filter(Measurement == "0")

initialmod <-lme(
  Fv.Fm ~ Micro * Needle_age * Elev * Prev_trmt,
  random = ~1 | Seedling_ID,
  data = initial
)
anova(initialmod)

drrr_df <- dur %>%
  filter(Measurement == 0 | Measurement == 7 | Measurement == "rec1" |Measurement == "rec2" | Measurement == "rec7") %>%
  filter(Trmt_2 == "heat") %>%
  group_by(Seed_ID, Seedling_ID, Prev_trmt, Needle_age, Micro, Elev, Genotype, Measurement) %>%
  summarize(mean_FvFm = mean(Fv.Fm, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Measurement, values_from = mean_FvFm, names_prefix = "M")

head(drrr_df)

# Resilience(Ingrish) = (Recovery - Stress) / (Initial - Stress)

drrr_df <- drrr_df %>%
  mutate(Resistance = M7 / M0) %>%
  mutate(Day1_Recovery = Mrec1 / M7) %>%
  mutate(Day7_Recovery = Mrec7 / M7) %>%
  mutate(Day1_Resilience = (Mrec1 - M7)/ (M0 - M7)) %>%
  mutate(Day2_Resilience = (Mrec2 - M7)/ (M0 - M7)) %>%
  mutate(Day7_Resilience = (Mrec7 - M7)/ (M0 - M7)) 

head(drrr_df)

write.csv(drrr_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_rrr.csv", row.names=FALSE)

# Recovery only for curve fitting
drec_long <- dur %>%
  filter(Measurement %in% c("7", "rec1", "rec2", "rec7")) %>%
  filter(Trmt_2 %in% ("heat")) %>%
  mutate(
    time = case_when(
      Measurement == "7"   ~ 0,
      Measurement == "rec1" ~ 1,
      Measurement == "rec2" ~ 2,
      Measurement == "rec7" ~ 7
    )
  )

head(drec_long)

write.csv(drec_long, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_rec_long.csv", row.names=FALSE)

# Control summary for recovery curve fitting
durcontrol_sum <- dur %>%
  group_by(Elev, Needle_age, Prev_trmt, Measurement, Trmt_2) %>%
  summarise(control_FvFm = mean(Fv.Fm, na.rm = TRUE), .groups = "drop") %>%
  filter(Measurement == "rec7", Trmt_2 == "control") %>%
  select(Elev, Needle_age, Prev_trmt, control_FvFm)

write.csv(durcontrol_sum, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_control_summary.csv", row.names=FALSE)

