library(tidyverse)
library(lubridate)
library(dplyr)
library(knitr)
library(kableExtra)
library(emmeans)
library(nlme)
library(car)

# Read in data
drrr_df <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_rrr.csv")
drec_long <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_rec_long.csv")
dcontrol_sum <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_control_summary.csv")

# Linear mixed effects modeling
dresistmodel<-lme(
  Resistance ~ Needle_age * Prev_trmt* Elev,
  random = ~ 1 | Seed_ID,
  data=drrr_df,
  na.action = na.exclude)

anova(dresistmodel)

dd1resil <- lme(
  Day1_Resilience ~ Needle_age *Prev_trmt * Elev,
  random = ~ 1 | Seed_ID,
  data=drrr_df,
  na.action = na.exclude
)

anova(dd1resil)

dd7resil <- lme(
  Day7_Resilience ~ Needle_age * Prev_trmt * Elev,
  random = ~ 1 | Seed_ID,
  data=drrr_df,
  na.action = na.exclude
)

anova(dd7resil)

# Emmeans analysis
dresisemm <- emmeans(dresistmodel, ~  Needle_age * Prev_trmt * Elev)
write.csv(dresisemm, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_resisemm.csv", row.names=FALSE)


dd1resilemm <- emmeans(dd1resil, ~  Needle_age * Prev_trmt * Elev)
write.csv(dd1resilemm, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_d1resilemm.csv", row.names=FALSE)

dd7resilemm <- emmeans(dd7resil, ~  Needle_age * Prev_trmt * Elev)
write.csv(dd7resilemm, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_d7resilemm.csv", row.names=FALSE)




# Log linear recovery
source("~/Masters/Growth_Chamber/scripts/functions/log_linear_recovery_function.R")
head(drec_long)

dres <- log_linear_recovery(
  data = drec_long,
  group_vars = c("Seed_ID", "Needle_age", "Prev_trmt", "Elev", "Genotype"),
  time_var = "time",
  response_var = "Fv.Fm",
  control_target = dcontrol_sum,
  control_join_vars = c("Elev", "Needle_age", "Prev_trmt")
)

dres$target_tab %>% filter(is.finite(days_to_recovery))

head(dres$fit_tab)        
head(dres$fit_curves)     
View(dres$target_tab)     

converged_ids <- dres$target_tab %>%
  filter(is.finite(days_to_recovery), days_to_recovery <= 100) %>%
  pull(Seed_ID)   

fit_curves_final <- dres$fit_curves %>%
  filter(Seed_ID %in% converged_ids)

fit_curves_plot <- fit_curves_final %>%
  mutate(Elev_Micro = interaction(Elev, Micro, sep = ":"))

mean_curves <- fit_curves_final %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control","mod","high"), labels = c("control","moderate","high"))) %>%
  group_by(Needle_age, Prev_trmt, Elev, time) %>%
  summarise(
    meanfitted = mean(fitted, na.rm = TRUE),
    se_fitted = sd(fitted, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

control_ref <- dres$target_tab %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control","mod","high"), labels = c("control","moderate","high"))) %>%
  distinct(Elev, Prev_trmt, Needle_age, target) 

ggplot(mean_curves, aes(x = time, y = meanfitted, color = Elev_Micro, fill = Elev_Micro)) +
  geom_ribbon(aes(ymin = meanfitted - se_fitted, ymax = meanfitted + se_fitted),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  facet_grid(Needle_age ~ Prev_trmt) +
  geom_hline(data = control_ref, aes(yintercept = target, color = Elev_Micro),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  labs(x = "Days since heat treatment", y = "Fv/Fm",
       color = "Elevation : Microbiome", fill = "Elevation : Microbiome") +
  ylim(0,0.9) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom")

ggplot(mean_curves, aes(x = time, y = meanfitted, color = Prev_trmt, fill = Prev_trmt)) +
  geom_ribbon(aes(ymin = meanfitted - se_fitted, ymax = meanfitted + se_fitted),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(data = control_ref, aes(yintercept = target, color = Prev_trmt),
             linetype = "dashed", linewidth = 0.7, show.legend = FALSE) + 
  facet_grid(Needle_age ~ Elev,
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  labs(x = "Days since heat treatment", y = "Fv/Fm",
       color = "Previous Treatment", fill = "Previous Treatment") +
  scale_color_manual(values = c("high" = "red", "moderate" = "orange", "control" = "green")) +
  scale_fill_manual(values = c("high" = "red", "moderate" = "orange", "control" = "green")) +
  theme_bw(base_size = 13) +
  ylim(0,0.85) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position="bottom") +
  ggsave("~/Masters/Growth_Chamber/results/figures/dur_response_recovery.png")


# Recovery metric models
dres$target_tab %>%
  summarise(
    n_na   = sum(is.na(days_to_recovery)),
    n_inf  = sum(is.infinite(days_to_recovery)),
    n_nan  = sum(is.nan(days_to_recovery))
  )

dres$target_tab <- dres$target_tab %>%
  filter(is.finite(days_to_recovery), days_to_recovery <= 100)

dmod_day <- lme(
  days_to_recovery ~ Needle_age * Prev_trmt * Elev,
  random = ~1 | Seed_ID,
  data = dres$target_tab,
  na.action = na.exclude
)

anova(dmod_day)

aov_tab <- anova(dmod_day) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  rename(
    numDF  = numDF,
    denDF  = denDF,
    `F-value` = `F-value`,
    `p-value` = `p-value`
  ) %>%
  mutate(
    `p-value` = case_when(
      `p-value` < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(`p-value`, 3))
    ),
    `F-value` = round(`F-value`, 2)
  )

tab <- aov_tab %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated day of recovery",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE)

tab


drec_rate <- lme(
  B ~ Needle_age * Prev_trmt * Elev,
  random = ~1 | Seed_ID,
  data = dres$fit_tab,
  na.action = na.exclude
)

anova(drec_rate)

aov_tab2 <- anova(rec_rate) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  rename(
    numDF  = numDF,
    denDF  = denDF,
    `F-value` = `F-value`,
    `p-value` = `p-value`
  ) %>%
  mutate(
    `p-value` = case_when(
      `p-value` < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(`p-value`, 3))
    ),
    `F-value` = round(`F-value`, 2)
  )

tab2 <- aov_tab2 %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated recovery rate coefficient",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE)

tab2


# Emmeans analysis
demm_mod_day <- emmeans(dmod_day, ~ Needle_age *  Prev_trmt * Elev)
demm_mod_day
write.csv(demm_mod_day, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_emm_mod_day.csv", row.names=FALSE)

demm_rec_rate <- emmeans(drec_rate, ~ Needle_age * Prev_trmt * Elev)
write.csv(demm_rec_rate, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/dur_emm_rec_rate.csv", row.names=FALSE)


# Pairwise contrasts
contrast(demm_mod_day, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()

demm_mod_day_df <- demm_mod_day %>% as_tibble()

# Plot EMMs
ggplot(demm_mod_day_df, aes(x = Elev, y = (emmean), 
                           color = Prev_trmt, group = Prev_trmt)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = (lower.CL), ymax = (upper.CL)),
                width = 0.2, position = position_dodge(0.3)) +
  facet_wrap(~Needle_age)+
  labs(
    x     = "Elevation",
    y     = "Days to Recovery ± 95% CI",
    color = "Previous Treatment"
  ) +
  theme_bw() +
  theme(legend.position = "top")




