library(tidyverse)
library(lubridate)
library(dplyr)

flu<-read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_full_final.csv")
head(flu)

flu<- flu %>%
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

flu$Seed_ID <- paste(flu$Seedling_ID, flu$ID2)

head(flu)

# Check if mean initial Fv/Fm differs between groups
initial <- flu %>%
  filter(Measurement == "25")

initialmod <-lme(
  Fv.Fm ~ Micro * Needle_age * Elev,
  random = ~1 | Seedling_ID,
  data = initial
)
anova(initialmod)


rrr_df <- flu %>%
  filter(Measurement == 25 | Measurement == 45 | Measurement ==50 | Measurement == "rec1" |Measurement == "rec2" | Measurement == "rec7") %>%
  filter(Treatment == "high") %>%
  group_by(Seed_ID, Seedling_ID, Needle_age, Micro, Elev, Genotype, Measurement) %>%
  summarize(mean_FvFm = mean(Fv.Fm, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Measurement, values_from = mean_FvFm, names_prefix = "M")

head(rrr_df)

# Resilience(Ingrish) = (Recovery - Stress) / (Initial - Stress)

rrr_df <- rrr_df %>%
  mutate(Resistance = M50 / M25) %>%
  mutate(Resistance45 = M45 / M25) %>%
  mutate(Day1_Recovery = Mrec1 / M50) %>%
  mutate(Day2_Recovery = Mrec2 / M50) %>%
  mutate(Day7_Recovery = Mrec7 / M50) %>%
  mutate(Day1_Resilience = (Mrec1 - M50)/ (M25 - M50)) %>%
  mutate(Day2_Resilience = (Mrec2 - M50)/ (M25 - M50)) %>%
  mutate(Day7_Resilience = (Mrec7 - M50)/ (M25 - M50)) 


write.csv(rrr_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_rrr.csv", row.names=FALSE)


# Recovery only for curve fitting

rec_long <- flu %>%
  filter(Measurement %in% c("50", "rec1", "rec2", "rec7")) %>%
  filter(Treatment %in% ("high")) %>%
  mutate(
    time = case_when(
      Measurement == "50"   ~ 0,
      Measurement == "rec1" ~ 1,
      Measurement == "rec2" ~ 2,
      Measurement == "rec7" ~ 7
    )
  )

head(rec_long)

write.csv(rec_long, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_rec_long.csv", row.names=FALSE)

# Control summary for recovery curve fitting
control_sum <- flu %>%
  group_by(Elev, Micro, Needle_age, Measurement, Treatment) %>%
  summarise(control_FvFm = mean(Fv.Fm, na.rm = TRUE), .groups = "drop") %>%
  filter(Measurement == "rec7", Treatment == "control") %>%
  select(Elev, Micro, Needle_age, control_FvFm)

write.csv(control_sum, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_control_summary.csv", row.names=FALSE)
