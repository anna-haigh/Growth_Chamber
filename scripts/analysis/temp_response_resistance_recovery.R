library(tidyverse)
library(lubridate)
library(dplyr)
library(knitr)
library(kableExtra)
library(emmeans)
library(nlme)
library(car)

rrr_df <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_rrr.csv")
rec_long <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_rec_long.csv")
control_sum <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_control_summary.csv")


# Linear mixed effects modeling
resistmodel<-lme(
  Resistance ~ Needle_age * Elev * Micro,
  random = ~ 1 | Seedling_ID,
  data=rrr_df,
  na.action = na.omit)

anova(resistmodel)

aov_res <- anova(resistmodel) %>%
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

tab_res <- aov_res %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated resistance to heat treatment",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/temp_resistance_anova.pdf")


resist45model<-lme(
  Resistance45 ~ Needle_age * Elev * Micro,
  random = ~ 1 | Seedling_ID,
  data=rrr_df,
  na.action = na.omit)

anova(resist45model)

d1resil <- lme(
  Day1_Resilience ~ Needle_age * Elev * Micro,
  random = ~ 1 | Seedling_ID,
  data=rrr_df,
  na.action = na.omit
)

anova(d1resil)

aov_resil <- anova(d1resil) %>%
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

tab_resil <- aov_resil %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated resilience to heat treatment",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/temp_resilience_anova.pdf")


d7resil <- lme(
  Day7_Resilience ~ Needle_age * Elev * Micro,
  random = ~ 1 | Seed_ID,
  data=rrr_df,
  na.action = na.omit
)

anova(d7resil)


# Emmeans analysis
resisemm <- emmeans(resistmodel, ~  Needle_age * Elev * Micro)
write.csv(resisemm, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_resisemm.csv", row.names=FALSE)

View(contrast(resisemm, method = "pairwise") %>%
       summary(infer = TRUE) %>%
       as_tibble())

resisemm_df <- as.data.frame(resisemm)

ggplot(resisemm_df, aes(x = Micro, y = emmean, color = Needle_age)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Elev)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Mean resistance to heat treatment"), 
       color = "Microbiome") +
  geom_hline(yintercept= 1, lty="dashed", color="gray34")+
  theme_bw() + 
  theme(legend.position = c(.1,.1))

d1resilemm <- emmeans(d1resil, ~  Needle_age * Elev * Micro)
write.csv(d1resilemm, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_d1resilemm.csv", row.names=FALSE)

d1resilemm_df <- as.data.frame(d1resilemm)

ageresil <- emmeans(d1resil, ~ Needle_age)
ageresil

contrast(ageresil, method = "pairwise") %>%
       summary(infer = TRUE) %>%
       as_tibble()

elevresil <- emmeans(d1resil, ~ Elev)
elevresil

contrast(elevresil, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()

agemicroresil <- emmeans(d1resil, ~ Needle_age * Micro)
agemicroresil

contrast(agemicroresil, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()


ggplot(d1resilemm_df, aes(x = Elev, y = emmean, color = Micro)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Mean resilience to heat treatment 1 day post-stress"), 
       color = "Microbiome") +
  geom_hline(yintercept= 1, lty="dashed", color="gray34")+
  theme_bw() + 
  theme(legend.position = c(.1,.1))


d7resilemm <- emmeans(d7resil, ~  Needle_age * Elev * Micro)
d7resilemm

d7resilemm_df <- as.data.frame(d7resilemm)

ggplot(d7resilemm_df, aes(x = Elev, y = emmean, color = Micro)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Mean resilience to heat treatment 7 days post-stress"), 
       color = "Microbiome") +
  geom_hline(yintercept= 1, lty="dashed", color="gray34")+
  theme_bw() + 
  theme(legend.position = c(.1,.1))



# Log linear recovery
source("~/Masters/Growth_Chamber/scripts/functions/log_linear_recovery_function.R")
head(rec_long)

res <- log_linear_recovery(
  data = rec_long,
  group_vars = c("Seed_ID","Seedling_ID", "Needle_age", "Micro", "Elev", "Genotype"),
  time_var = "time",
  response_var = "Fv.Fm",
  control_target = control_sum,
  control_join_vars = c("Elev", "Micro", "Needle_age")
)

res$target_tab %>% filter(is.finite(days_to_recovery))

head(res$fit_tab)        
head(res$fit_curves)     
head(res$target_tab)     

converged_ids <- res$target_tab %>%
  filter(is.finite(days_to_recovery))   

fit_curves_final <- res$fit_curves %>%
  filter(Seed_ID %in% converged_ids)

fit_curves_plot <- fit_curves_final %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
  mutate(Elev_Micro = interaction(Elev, Micro, sep = " + "))

mean_curves <- fit_curves_plot %>%
  group_by(Needle_age, Elev_Micro, time) %>%
  summarise(
    meanfitted = mean(fitted, na.rm = TRUE),
    se_fitted = sd(fitted, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

head(mean_curves)

control_ref <- res$target_tab %>%
  distinct(Elev, Micro, Needle_age, target) %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
  mutate(Elev_Micro = interaction(Elev, Micro, sep = " + "))

rec_figure <- ggplot(mean_curves, aes(x = time, y = meanfitted, color = Elev_Micro, fill = Elev_Micro)) +
  geom_ribbon(aes(ymin = meanfitted - se_fitted, ymax = meanfitted + se_fitted),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Needle_age) +
  geom_hline(data = control_ref, aes(yintercept = target, color = Elev_Micro),
             linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
  labs(x = "Days Since Heat Treatment", y = "Fv/Fm",
       color = "Elevation + Microbiome", fill = "Elevation + Microbiome") +
  ylim(0,0.9) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = c(0.8,0.2),
        panel.grid.minor = element_blank())

rec_figure
ggsave("temp_response_recovery.png", rec_figure, path = "~/Masters/Growth_Chamber/results/figures")


# Recovery metric models
mod_day <- lme(
  days_to_recovery ~ Needle_age * Micro * Elev,
  random = ~1 | Seedling_ID,
  data = res$target_tab,
  na.action = na.omit
)

anova(mod_day)

aov_tab <- anova(mod_day) %>%
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
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/temp_recovery_day_anova.pdf")

tab


rec_rate <- lme(
  B ~ Needle_age * Micro * Elev,
  random = ~1 | Seedling_ID,
  data = res$fit_tab,
  na.action = na.omit
)

anova(rec_rate)

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
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/temp_recovery_rate_anova.pdf")

tab2

# Emmeans analysis
emm_mod_day <- emmeans(mod_day, ~ Needle_age *  Micro * Elev)
write.csv(emm_mod_day, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_emm_mod_day.csv", row.names=FALSE)

View(contrast(emm_mod_day, method = "pairwise") %>%
       summary(infer = TRUE) %>%
       as_tibble())

agemodday<-emmeans(mod_day, ~ Needle_age)
agemodday
contrast(agemodday, method = "pairwise") %>%
    summary(infer = TRUE) %>%
    as_tibble()

emm_rec_rate <- emmeans(rec_rate, ~ Needle_age * Micro * Elev)
write.csv(emm_rec_rate, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/temp_emm_rec_rate.csv", row.names=FALSE)

View(contrast(emm_rec_rate, method = "pairwise") %>%
       summary(infer = TRUE) %>%
       as_tibble())

# Pairwise contrasts
contrast(emm_mod_day, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()

emm_mod_day_df <- emm_mod_day %>% as_tibble()

# Plot EMMs
ggplot(emm_mod_day_df, aes(x = Elev, y = (emmean), 
                           color = Micro, group = Micro)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = (lower.CL), ymax = (upper.CL)),
                width = 0.2, position = position_dodge(0.3)) +
  facet_wrap(~Needle_age)+
  labs(
    x     = "Elevation",
    y     = "Days to Recovery ± 95% CI",
    color = "Microbial Inoculation"
  ) +
  theme_bw() +
  theme(legend.position = "top")


