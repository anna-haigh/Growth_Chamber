library(tidyverse)
library(lubridate)
library(dplyr)
library(knitr)
library(kableExtra)
library(webshot2)
library(emmeans)
library(nlme)

# Read in data
dur_dat <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/duration_curve_data.csv")

ggplot(dur_dat, aes(x=R2)) + 
  geom_density()  
mean(dur_dat$R2, na.rm=TRUE)


# Linear mixed effects models
dcritmod <- lme(
  Dcrit ~ Prev_trmt*Needle_age*Elev,
  random = ~ 1 | Seed_ID,
  data=dur_dat,
  na.action = na.exclude
)
anova(dcritmod)

d15mod <- lme(
  D15 ~ Prev_trmt*Needle_age*Elev,
  random = ~ 1 | Seed_ID,
  data=dur_dat,
  na.action = na.exclude
)
anova(d15mod)

d50mod <- lme(
  D50 ~ Prev_trmt*Needle_age*Elev,
  random = ~ 1 | Seed_ID,
  data=dur_dat,
  na.action = na.exclude
)
anova(d50mod)

# Residual checks 
dur_dat<- dur_dat %>%
  mutate(res1=residuals(d50mod, type="pearson"),
         res2=residuals(dcritmod, type="pearson"),
         res3=residuals(d15mod, type="pearson"))

ggplot(dur_dat, aes(Needle_age, res1))+
  geom_boxplot()
ggplot(dur_dat, aes(Elev, res1))+
  geom_boxplot()
ggplot(dur_dat, aes(x=res1)) + 
  geom_density()

ggplot(dur_dat, aes(Needle_age, res2))+
  geom_boxplot()
ggplot(dur_dat, aes(Elev, res2))+
  geom_boxplot()
ggplot(dur_dat, aes(x=res2)) + 
  geom_density()

ggplot(dur_dat, aes(Needle_age, res3))+
  geom_boxplot()
ggplot(dur_dat, aes(Elev, res3))+
  geom_boxplot()
ggplot(dur_dat, aes(x=res3)) + 
  geom_density()


# ANOVA tables
aov_table <- anova(d50mod) %>%
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

table1 <- aov_table %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated D{50}",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/D50_anova.pdf")


aov_table2 <- anova(dcritmod) %>%
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

table2 <- aov_table2 %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated Dcrit",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/Dcrit_anova.pdf")

aov_table3 <- anova(d15mod) %>%
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

table3 <- aov_table3 %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for estimated D15",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/D15_anova.pdf")

# Emmeans analysis

emmd50 <- emmeans(d50mod, ~  Prev_trmt * Needle_age * Elev)
emmd50

emmd50_df <- as.data.frame(emmd50)

write.csv(emmd50_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmd50.csv", row.names=FALSE)


ggplot(emmd50_df, aes(x = Elev, y = emmean, color = Prev_trmt)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean D "[50]*""), 
       color = "Previous treatment") +
  theme_bw() + 
  theme(legend.position = c(.1,.1))


emmdcrit <- emmeans(dcritmod, ~  Prev_trmt * Needle_age * Elev)
emmdcrit

emmdcrit_df <- as.data.frame(emmdcrit)

write.csv(emmdcrit_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmdcrit.csv", row.names=FALSE)


ggplot(emmdcrit_df, aes(x = Elev, y = emmean, color = Prev_trmt)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean T "[crit]*""), 
       color = "Previous treatment") +
  theme_bw() + 
  theme(legend.position = c(.1,.9))

emmd15 <- emmeans(d15mod, ~  Needle_age * Elev * Prev_trmt)
emmd15

emmd15_df <- as.data.frame(emmd15)

write.csv(emmd15_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmd15.csv", row.names=FALSE)


ggplot(emmd15_df, aes(x = Elev, y = emmean, color = Prev_trmt)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean T "[15]*""), 
       color = "Previous Treatment") +
  theme_bw() + 
  theme(legend.position = c(.1,.9))
