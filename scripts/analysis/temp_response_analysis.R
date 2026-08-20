library(tidyverse)
library(lubridate)
library(dplyr)
library(knitr)
library(kableExtra)
library(webshot2)
library(emmeans)
library(nlme)

curve_dat <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/response_curve_data.csv")
head(curve_dat)

# Linear mixed effects models

critmod <- lme(
  tcrit ~ Needle_age*Elev*Micro,
  random = ~ 1 | Seedling_ID,
  data=curve_dat
)
anova(critmod)


t50mod <- lme(
  ctmax ~ Needle_age*Elev*Micro,
  random = ~ 1 | Seedling_ID,
  data=curve_dat
)
anova(t50mod)


t15mod <- lme(
  t15 ~ Needle_age*Elev*Micro,
  random = ~ 1 | Seedling_ID,
  data=curve_dat
)
anova(t15mod)


# Residual checks 

curve_dat<- curve_dat %>%
  mutate(res1=residuals(t50mod, type="pearson"),
         res2=residuals(critmod, type="pearson"),
         res3=residuals(t15mod, type="pearson"))

ggplot(curve_dat, aes(Needle_age, res1))+
  geom_boxplot()
ggplot(curve_dat, aes(Micro, res1))+
  geom_boxplot()
ggplot(curve_dat, aes(Elev, res1))+
  geom_boxplot()
ggplot(curve_dat, aes(x=res1)) + 
  geom_density()

ggplot(curve_dat, aes(Needle_age, res2))+
  geom_boxplot()
ggplot(curve_dat, aes(Micro, res2))+
  geom_boxplot()
ggplot(curve_dat, aes(Elev, res2))+
  geom_boxplot()
ggplot(curve_dat, aes(x=res2)) + 
  geom_density()

ggplot(curve_dat, aes(Needle_age, res3))+
  geom_boxplot()
ggplot(curve_dat, aes(Micro, res3))+
  geom_boxplot()
ggplot(curve_dat, aes(Elev, res3))+
  geom_boxplot()
ggplot(curve_dat, aes(x=res3)) + 
  geom_density()

# ANOVA tables
aov_table <- anova(t50mod) %>%
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
    caption = "Linear mixed-effects model ANOVA table for estimated T50",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/T50_anova.pdf")


aov_table2 <- anova(critmod) %>%
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
    caption = "Linear mixed-effects model ANOVA table for estimated Tcrit",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/Tcrit_anova.pdf")

aov_table3 <- anova(t15mod) %>%
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
    caption = "Linear mixed-effects model ANOVA table for estimated T15",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE) %>%
  save_kable("~/Masters/Growth_Chamber/results/tables/T15_anova.pdf")


# Emmeans analysis

emmt50 <- emmeans(t50mod, ~  Needle_age * Elev * Micro)
emmt50

emmt50_df <- as.data.frame(emmt50)

write.csv(emmt50_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmt50.csv", row.names=FALSE)


ggplot(emmt50_df, aes(x = Elev, y = emmean, color = Micro)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean T "[50]*""), 
       color = "Microbiome") +
  theme_bw() + 
  theme(legend.position = c(.1,.1))


emmtcrit <- emmeans(critmod, ~  Needle_age * Elev * Micro)
emmtcrit

emmtcrit_df <- as.data.frame(emmtcrit)

write.csv(emmtcrit_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmtcrit.csv", row.names=FALSE)


ggplot(emmtcrit_df, aes(x = Elev, y = emmean, color = Micro)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean T "[crit]*""), 
       color = "Microbiome") +
  theme_bw() + 
  theme(legend.position = c(.1,.9))

emmt15 <- emmeans(t15mod, ~  Needle_age * Elev * Micro)
emmt15

emmt15_df <- as.data.frame(emmt15)

write.csv(emmt15_df, "~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/emmt15.csv", row.names=FALSE)


ggplot(emmt15_df, aes(x = Elev, y = emmean, color = Micro)) +
  geom_point(size = 3, position = position_dodge(0.7)) +
  facet_wrap(~Needle_age)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(0.7)) +
  labs(x = "Seed Source Elevation", y = expression("Estimated mean T "[15]*""), 
       color = "Microbiome") +
  theme_bw() + 
  theme(legend.position = c(.1,.9))

# Pairwise based on significant stats
etcrit<- emmeans(critmod, ~ Needle_age)
etcrit

contrast(etcrit, method = "pairwise") %>%
    summary(infer = TRUE) %>%
    as_tibble()

et15<- emmeans(t15mod, ~ Needle_age)
et15
contrast(et15, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()
