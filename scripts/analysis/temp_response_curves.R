library(tidyverse)
library(lubridate)
library(dplyr)
library(car)

# Read in data
flu<-read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_full_final.csv")

head(flu)

flu<- flu %>%
  filter(Measurement %in% c("25","30","35","40","45","50")) %>%
  filter(Seedling_ID %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")) %>%
  mutate(., 
         Elev = factor(Elev),
         Genotype = factor(Genotype),
         Micro = factor(Micro),
         Needle_age = factor(Needle_age),
         Measurement = as.numeric(Measurement),
         Date = mdy(Date),
         Fo = as.numeric(Fo),
         Fv = as.numeric(Fv),
         Fm = as.numeric(Fm)
  )

flu$Seed_ID <- paste(flu$Seedling_ID, flu$ID2)

head(flu)

# Check and process data to run through function
print(colSums(is.na(flu[, c("Fv.Fm", "Measurement", "Needle_age", "Seed_ID")])))

# Check species levels
cat("Unique needle ages:", paste(unique(flu$Needle_age), collapse = ", "), "\n")
cat("Unique seed IDs:", paste(sort(unique(flu$Seed_ID)), collapse = ", "), "\n\n")

# Check data for each individual
table_counts <- table(flu$Needle_age, flu$Seed_ID)
print(table_counts)

# Check for individuals with too few data points (<5)
individual_counts <- aggregate(Fv.Fm ~ Needle_age + Seed_ID, data = flu, FUN = length)
names(individual_counts)[3] <- "n_points"
problem_individuals <- individual_counts[individual_counts$n_points < 5, ]
if(nrow(problem_individuals) > 0){
  print(problem_individuals)
} else {
  cat("None - all individuals have sufficient data points\n")
}

# Remove rows with missing values
cat("Original rows:", nrow(flu), "\n")
flu <- flu[!is.na(flu$Fv.Fm) & !is.na(flu$Measurement) & 
             !is.na(flu$Needle_age) & !is.na(flu$Seed_ID), ]
cat("After removing NAs:", nrow(flu), "\n\n")

# Fix Fv.Fm values that are exactly 0 or 1 (causes logit errors)
extreme_zero <- flu$Fv.Fm == 0
extreme_one <- flu$Fv.Fm == 1

if(sum(extreme_zero, na.rm = TRUE) > 0){
  cat("Replacing", sum(extreme_zero, na.rm = TRUE), "values of 0 with 0.001\n")
  flu$Fv.Fm[extreme_zero] <- 0.001
}

if(sum(extreme_one, na.rm = TRUE) > 0){
  cat("Replacing", sum(extreme_one, na.rm = TRUE), "values of 1 with 0.999\n")
  flu$Fv.Fm[extreme_one] <- 0.999
}


# Remove values outside valid range
out_of_range <- flu$Fv.Fm < 0 | flu$Fv.Fm > 1
if(sum(out_of_range, na.rm = TRUE) > 0){
  cat("Removing", sum(out_of_range, na.rm = TRUE), "rows with Fv.Fm outside 0-1 range\n")
  flu <- flu[!out_of_range, ]
}


# Run data through temperature response function

source("~/Masters/Growth_Chamber/scripts/functions/temp_response_function.R")


results <- Thermal_tol2(
  fv_fm = Fv.Fm, 
  temperature = Measurement,
  individual = Seed_ID,
  data = flu, 
  print_graph = FALSE  # Will display graphs but not save them
)

thermal_df <- results[, c("ID", 
                          "ctmax", "ctmax.lci", "ctmax.uci",
                          "tcrit", "tcrit.lci", "tcrit.uci",
                          "t15", "t15.lci", "t15.uci")]

head(thermal_df)

thermal_df <- thermal_df %>%
  rename(., Seed_ID = ID)
head(thermal_df)

seed_sum <- flu %>%
  subset(., select = -c(Measurement,Date,Time,Fo,Fv,Fm,Fv.Fo,Fv.Fm))

seed_sum <- distinct(seed_sum)
head(seed_sum)

curve_dat<-merge(thermal_df, seed_sum, by = c("Seed_ID"))

write.csv(curve_dat, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/response_curve_data.csv", row.names=FALSE)



#### END SCRIPT












