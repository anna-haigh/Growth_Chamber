library(tidyverse)
library(here)
library(lubridate)
library(car)

# Read the data
flur2<- read.csv("~/Masters/Growth_chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_treatment2_full.csv")

# Clean data for duration model
# Logistic function where asymptotes exist at Fv/Fm = 0.8, Fv/Fm = 0

Fv.Fm = 0.8 / (1 + exp(time))

dur <- flur2 %>%
  filter(Measurement %in% c("0", "3", "4", "5", "6", "7")) %>%
  filter(Trmt_2 %in% ("heat")) %>%
  mutate(
    time = case_when(
      Measurement == "0" ~ 0,
      Measurement == "3" ~ 3,
      Measurement == "4" ~ 4,
      Measurement == "5" ~ 5,
      Measurement == "6" ~ 6,
      Measurement == "7" ~ 7
    )
  )

dur$Seed_ID <- paste(dur$Seedling_ID, dur$ID2)

head(dur)

# Remove rows with missing values in key columns
cat("Original rows:", nrow(dur), "\n")
dur <- dur[!is.na(dur$Fv.Fm) & !is.na(dur$time) & 
             !is.na(dur$Needle_age) & !is.na(dur$Seed_ID), ]
cat("After removing NAs:", nrow(dur), "\n\n")

# Fix Fv.Fm values that are exactly 0 or 1 (causes logit errors)
extreme_zero <- dur$Fv.Fm == 0
extreme_one <- dur$Fv.Fm == 1

if(sum(extreme_zero, na.rm = TRUE) > 0){
  cat("Replacing", sum(extreme_zero, na.rm = TRUE), "values of 0 with 0.001\n")
  dur$Fv.Fm[extreme_zero] <- 0.001
}

if(sum(extreme_one, na.rm = TRUE) > 0){
  cat("Replacing", sum(extreme_one, na.rm = TRUE), "values of 1 with 0.999\n")
  dur$Fv.Fm[extreme_one] <- 0.999
}

# Remove values outside valid range (if any)
out_of_range <- dur$Fv.Fm < 0 | dur$Fv.Fm > 1
if(sum(out_of_range, na.rm = TRUE) > 0){
  cat("Removing", sum(out_of_range, na.rm = TRUE), "rows with Fv.Fm outside 0-1 range\n")
  dur <- dur[!out_of_range, ]
}

# Check individuals with insufficient data points
individual_counts <- aggregate(Fv.Fm ~ Needle_age + Seed_ID, data = dur, FUN = length)
names(individual_counts)[3] <- "n_points"
low_count_individuals <- individual_counts[individual_counts$n_points < 5, ]


# Duration model
source("~/Masters/Growth_Chamber/scripts/functions/duration_function.R")

resultsdur <- duration_tol(
  fv_fm = Fv.Fm, 
  time = time,
  individual = Seed_ID,
  data = dur, 
  print_graph = FALSE  # Will display graphs but not save them
)

head(resultsdur)

dur_df <- as.data.frame(resultsdur)

dur_df <- dur_df %>%
  rename(., Seed_ID = ID)
head(dur_df)

seedling_sum <- dur %>%
  subset(., select = -c(Measurement,Date,time,Time,Fo,Fv,Fm,Fv.Fo,Fv.Fm))

seedling_sum <- distinct(seedling_sum)

dur_dat<-merge(dur_df, seedling_sum, by = ("Seed_ID")) %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control", "mod", "high")))

View(dur_dat)

write.csv(dur_dat, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/duration_curve_data.csv", row.names=FALSE)

