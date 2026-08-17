library(tidyverse)
library(here)
library(lubridate)
library(lmerTest)
library(emmeans)
library(nlme)
library(knitr)
library(kableExtra)

Amax <- read.csv("haigh_growth_chamber/full_labeled.csv")
dur<-read.csv(here("Chlorophyll_fluorescence/chl_fluor_treatment2_full.csv"))

Amax <- Amax %>% rename_at('seedling_id', ~'Seedling_ID')
Amax <- Amax %>%
  subset(., amax_corrected >=0) %>%
  subset(., gsw_corrected >=0) %>%
  subset(., select = c(Seedling_ID, treatment_1, treatment_2, measurement,leaf_area,amax_corrected, gsw_corrected))
head(Amax)

dur<- dur %>%
  subset(., select = c(Seedling_ID, Elev, Genotype, Micro))
merged <- merge(Amax, dur, by="Seedling_ID")
head(merged)

Amax_sum <- merged %>%
  group_by(Seedling_ID, treatment_1, treatment_2, measurement, amax_corrected, gsw_corrected, Elev, Genotype,
           Micro) %>%
  summarise(
    leaf_area=mean(leaf_area),
    amax_corrected=mean(amax_corrected),
    gsw_corrected=mean(gsw_corrected),
    .groups="drop"
  )

Amax_sum <- Amax_sum %>%
  mutate(.,
         Seedling_ID = factor(Seedling_ID),
         Elev = factor(Elev),
         Genotype = factor(Genotype),
         Micro = factor(Micro),
         measurement = factor(measurement),
         treatment_1 = factor(treatment_1),
         treatment_2 = factor(treatment_2),
         amax_corrected = as.numeric(amax_corrected),
         gsw_corrected = as.numeric(gsw_corrected))
head(Amax_sum)


ggplot(Amax_sum, aes(x=measurement, y=amax_corrected, fill=Micro)) +
  geom_boxplot() + 
  ylim(-5,15)

june<-Amax_sum %>%
  filter(measurement %in% c("06pre","06rec1", "06rec7","07pre"))

june <- june %>%
  mutate(measurement = factor(
    measurement,
    levels = c("06pre","06rec1", "06rec7","07pre"),
    labels = c("pre", "rec01", "rec07", "rec14")
  ))

head(june)

ggplot(june, aes(x=measurement, y=amax_corrected, fill=treatment_1))+
  geom_boxplot()+
  ylim(-5,15) + 
  scale_fill_manual(values = c("high" = "red", "mod" = "orange", "control" = "green"))

ggplot(june, aes(x=measurement, y=gsw_corrected, fill=treatment_1))+
  geom_boxplot()+
  ylim(0,1.5) + 
  scale_fill_manual(values = c("high" = "red", "mod" = "orange", "control" = "green"))

##JULY IS NESTED
july <- Amax_sum %>%
  filter(measurement %in% c("07pre", "07rec1","07rec7"))
head(july)

july$treatment_1 <- factor(july$treatment_1, levels=c("control","mod","high"))

ggplot(july, aes(x=measurement, y=amax_corrected, fill=treatment_2))+
  geom_boxplot()+
  ylim(-5,15) +
  scale_fill_manual(values = c("heat" = "red","control" = "white"))

ggplot(july, aes(x=measurement, y=amax_corrected, fill=treatment_1))+
  geom_boxplot()+
  ylim(-2,12) +
  scale_fill_manual(values = c("high" = "red", "mod" = "orange", "control" = "green"))+
  labs(y=expression(A[max]~"("*mu*"mol CO"[2]~m^{-2}~s^{-1}*")"),
       x="Measurement",
       fill="Previous heat treatment")+
  theme_bw()+
  facet_wrap(~Elev)+
  theme(legend.position="top")
###########################################
# June Contrasts relative to control
june_sub <- june %>%
  filter(measurement %in% c("pre", "rec07")) %>%
  droplevels() %>%
  mutate(treatment_1 = factor(treatment_1, levels = c("control", "mod", "high")))

model_june <- lme(
  amax_corrected ~ measurement * treatment_1 * Elev * Micro,
  random = ~ 1 | Seedling_ID,
  data = june_sub,
  na.action = na.omit
)

summary(model_june)
anova(model_june)

emm_june <- emmeans(model_june, ~ treatment_1 * measurement | Elev * Micro)

diD_june <- contrast(
  emm_june,
  interaction = c("trt.vs.ctrl1", "revpairwise"),
  by = c("Elev", "Micro")
)

summary(diD_june, infer = c(TRUE, TRUE)) 
diD_june_df <- as.data.frame(summary(diD_june, infer = c(TRUE, TRUE)))

ggplot(diD_june_df, aes(x = treatment_1_trt.vs.ctrl1, y = estimate, color = treatment_1_trt.vs.ctrl1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
  facet_grid(Elev ~ Micro, labeller = labeller(.rows = label_both, .cols = label_both)) +
  scale_color_manual(values = c("mod - control" = "steelblue", "high - control" = "red")) +
  labs(
    x = "Treatment vs. Control",
    y = "Contrasts between treatment and control day 7 and initial Amax",
    color = "Contrast"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "none")
###############################
# July
july_sub <- Amax_sum %>%
  filter(measurement %in% c("07pre", "07rec1", "07rec7")) %>%
  droplevels() %>%
  mutate(
    measurement = factor(measurement,
                         levels = c("07pre", "07rec1", "07rec7"),
                         labels = c("pre", "rec01", "rec07")),
    treatment_1 = factor(treatment_1, levels = c("control", "mod", "high")),
    treatment_2 = factor(treatment_2, levels = c("control", "heat"))
  )

model_july <- lme(
  amax_corrected ~ measurement + treatment_2 * treatment_1 * Elev * Micro,
  random = ~ 1 | Seedling_ID,
  data = july_sub,
  na.action = na.omit
)

summary(model_july)

emm_july <- emmeans(model_july, ~ treatment_2 * measurement | treatment_1)

diD_july <- contrast(
  emm_july,
  interaction = c("trt.vs.ctrl1", "revpairwise"),
  by = "treatment_1"
)

summary(diD_july, infer = c(TRUE, TRUE))

diD_july_df <- as.data.frame(summary(diD_july, infer = c(TRUE, TRUE)))

diD_july_df

ggplot(diD_july_df, aes(x = measurement_revpairwise, y = estimate, color = measurement_revpairwise)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15,
                position = position_dodge(width = 0.4)) +
  facet_wrap(~ treatment_1) +
  scale_color_manual(values = c("rec01 - pre" = "goldenrod", "rec07 - pre" = "firebrick")) +
  labs(
    x = "Timepoint comparison",
    y = "Difference-in-differences estimate\n(heat change − control change)",
    title = "July Amax: heat effect on recovery relative to control, by previous treatment",
    color = "Comparison"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "none")
##########################################################
##model fitting june

head(june)

juner <- june %>%
  filter(measurement == "pre" | measurement == "rec01" | measurement == "rec07" | measurement == "rec14") %>%
  group_by(Seedling_ID, Micro, Elev, treatment_1, Genotype, measurement) %>%
  summarize(mean_amax = mean(amax_corrected, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(names_from = measurement, values_from = mean_amax, names_prefix = "M")

View(juner)

juner <- juner %>%
  mutate(Resistance = Mrec01 / Mpre) %>%
  mutate(Day7_Recovery = Mrec07 / Mrec01) %>%
  mutate(Day1_Resilience = Mrec01 / Mpre) %>%
  mutate(Day7_Resilience = Mrec07 / Mpre) %>%
  mutate(Resilience = (Mrec07 - Mrec01) / (Mpre - Mrec07)) 

junemod <- lme(
  Resistance ~ treatment_1*Elev*Micro,
  random = ~ 1 | Seedling_ID,
  data=juner,
  na.action = na.omit
)
anova(junemod)
coef(junemod)

emmjun <- emmeans(junemod, ~  treatment_1 * Elev * Micro)
emmjun

contrast(emmjun, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()

emmjun_df <- emmjun %>% as_tibble()

ggplot(emmjun_df, aes(x= treatment_1, y= emmean, color=Elev, group=Elev))+
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = (lower.CL), ymax = (upper.CL)),
                width = 0.2, position = position_dodge(0.3)) +
  facet_wrap(~Micro)+
  geom_hline(yintercept = 1, lty="dashed", color="gray34")+
  labs(
    x     = "Treatment",
    y     = "Resistance",
    color = "Elevation"
  ) +
  scale_color_manual(values =c("high"="red", "low"="lightblue"))+
  theme_bw() +
  theme(legend.position = "top")

junemod2 <- lme(
  Resilience ~ treatment_1*Elev,
  random = ~ 1 | Seedling_ID,
  data=juner,
  na.action = na.omit
)
anova(junemod2)
coef(junemod2)

emmjun2 <- emmeans(junemod2, ~  treatment_1 * Elev)
emmjun2

contrast(emmjun2, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble()

emmjun2_df <- emmjun2 %>% as_tibble()

ggplot(emmjun2_df, aes(x= treatment_1, y= emmean, color=Elev, group=Elev))+
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = (lower.CL), ymax = (upper.CL)),
                width = 0.2, position = position_dodge(0.3))+
  geom_hline(yintercept = 1, lty="dashed", color="gray34")+
  labs(
    x     = "Treatment",
    y     = "Resistance",
    color = "Elevation"
  ) +
  scale_color_manual(values =c("high"="red", "low"="lightblue"))+
  theme_bw() +
  theme(legend.position = "top")




juner7<- emmjun %>%
  contrast(
    method = list(
      "Change in Amax Control" = contr7 - contpre,
      "Change in Amax Moderate" = modr7 - modpre,
      "Change in Amax High" = highr7 - highpre
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
juner7


ggplot(juner7, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="Difference in pre-treatment and 7 day recovery",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)


contrelrec7 <- contr7 - contpre
modrelrec7 <- modr7 - modpre
highrelrec7<- highr7 - highpre



junerelrec1<- emmjun %>%
  contrast(
    method = list(
      "High temp. relative to control" =  highrelrec1 - contrelrec1,
      "Mod. temp. relative to control" = modrelrec1 - contrelrec1)
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
junerelrec1

ggplot(junerelrec1, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="Mean Change in Amax Before and After Treatment Relative to Control",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)


junerelrec7<- emmjun %>%
  contrast(
    method = list(
      "High temp. recovery" = highrelrec7 - contrelrec7,
      "Mod. temp. recovery" = modrelrec7 - contrelrec7
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
junerelrec7


ggplot(junerelrec7, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="7 day relative recovery (compared with control)",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)

junerelrec14<- emmjun %>%
  contrast(
    method = list(
      "High temp. recovery" = highrelrec14 - contrelrec14,
      "Mod. temp. recovery" = modrelrec14 - contrelrec14
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
junerelrec14


ggplot(junerelrec14, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="14 day relative recovery (compared with control)",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)

ggplot(juner14, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="Difference in pre-treatment and 14 day recovery",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)



juner14<- emmjun %>%
  contrast(
    method = list(
      "Change in Amax Control" = contr14 - contpre,
      "Change in Amax Moderate" = modr14 - modpre,
      "Change in Amax High" = highr14 - highpre
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
juner14

contrelrec14 <- contr14 - contpre
modrelrec14 <- modr14- modpre
highrelrec14<- highr14 - highpre

juner714<- emmjun %>%
  contrast(
    method = list(
      "Change in Amax Control" = contr14 - contr7,
      "Change in Amax Moderate" = modr14 - modr7,
      "Change in Amax High" = highr14 - highr7
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
juner714


ggplot(juner714, aes(x= contrast, y= estimate))+
  geom_errorbar(
    width =.15, lwd= .75,
    aes(ymin= lower.CL, ymax= upper.CL)
  )+
  geom_point(size= 3)+
  labs(
    y ="Difference in 14 and 7 day recovery",
    x =NULL
  )+
  theme_bw()+
  geom_hline(yintercept= 0, colour="grey34",lty= 2)

##Figure out how to combine relative recovery into one figure
junerelrec1 <- emmjun %>%
  contrast(
    method = list(
      "High temp. relative to control" = highrelrec1 - contrelrec1,
      "Mod. temp. relative to control" = modrelrec1 - contrelrec1
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble() %>%
  mutate(timepoint = "Post-treatment (Day 1)")

junerelrec7 <- emmjun %>%
  contrast(
    method = list(
      "High temp. relative to control" = highrelrec7 - contrelrec7,
      "Mod. temp. relative to control" = modrelrec7 - contrelrec7
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble() %>%
  mutate(timepoint = "7-Day Recovery")

junerelrec14 <- emmjun %>%
  contrast(
    method = list(
      "High temp. relative to control" = highrelrec14 - contrelrec14,
      "Mod. temp. relative to control" = modrelrec14 - contrelrec14
    )
  ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble() %>%
  mutate(timepoint = "14-Day Recovery")

# --- Combine into one data frame ---

june_all <- bind_rows(junerelrec1, junerelrec7, junerelrec14) %>%
  mutate(timepoint = factor(timepoint, 
                            levels = c("Post-treatment (Day 1)", 
                                       "7-Day Recovery", 
                                       "14-Day Recovery")))

# --- Combined faceted plot ---

ggplot(june_all, aes(x = contrast, y = estimate, colour = contrast)) +
  geom_errorbar(
    width = .15, lwd = .75,
    aes(ymin = lower.CL, ymax = upper.CL)
  ) +
  geom_point(size = 3) +
  facet_wrap(~ timepoint, ncol = 3) +
  labs(
    y = expression("Mean change in A"[max]*" relative to control"),
    x = NULL,
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_hline(yintercept = 0, colour = "grey34", lty = 2)





emmjunel <- emmeans(junemod, ~  measurement * Elev * Micro * treatment_1)
emmjunel

helincontpre<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helincontr1<- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helincontr7<- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helincontr14<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
lelincontpre<-c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
lelincontr1<- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
lelincontr7<- c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
lelincontr14<-c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helhighpre<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helhighr1<- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
helhighr7<- c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
helhighr14<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
lelhighpre<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
lelhighr1<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
lelhighr7<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
lelhighr14<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
helmodpre<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
helmodr1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)
helmodr7 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
helmodr14<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0)
lelmodpre<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
lelmodr1 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
lelmodr7 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
lelmodr14<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)


juner1elev<- emmjunel %>%
  contrast(
    method = list(
           "Change in Amax Control (high)" = helcontr1 - helcontpre,
           "Change in Amax Moderate (high)" = helmodr1 - helmodpre,
           "Change in Amax High (high)" = helhighr1 - helhighpre,
           "Change in Amax Control (low)" = lelcontr1 - lelcontpre,
           "Change in Amax Moderate (low)" = lelmodr1 - lelmodpre,
           "Change in Amax High (low)" = lelhighr1 - lelhighpre
         )
    ) %>%
  summary(infer = c(T, F)) %>%
  as_tibble()
juner1elev


head(july)



july %>% count(treatment_1, treatment_2)

# Check sample size and balance
nrow(july)
july %>% count(treatment_1, treatment_2, Elev, Micro, measurement)

# Check for any single-level factors
july %>% summarise(across(c(treatment_1, treatment_2, Elev, Micro, measurement), n_distinct))


julymod <- lme(
  amax_corrected ~ treatment_1 / treatment_2 +  Elev + Micro + measurement -1,
  random = ~ 1 | Seedling_ID,
  data = july,
  na.action = na.omit
)

anova(julymod)


aov2_table <- anova(julymod) %>%
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

aov2_table %>%
  kbl(
    caption = "Linear mixed-effects model ANOVA table for July A$_{max}$",
    col.names = c("Fixed Effect", "num df", "den df", "F", "p"),
    align = c("l", "c", "c", "c", "c"),
    booktabs = TRUE   # clean lines for LaTeX/PDF
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE)

emm_full <- emmeans(julymod, ~ treatment_2 | treatment_1 * measurement)
summary(emm_full)

emm_full_df <- as.data.frame(emm_full)

emm_full_df %>%
  filter(treatment_2 == "heat") %>%
  ggplot(aes(x = measurement, y = emmean, 
             color = treatment_1, group = treatment_1)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  scale_x_discrete(labels = c("07pre" = "Pre", "07rec1" = "Day 1", "07rec7" = "Day 7")) +
  scale_color_manual(values = c("control" = "steelblue", "high" = "firebrick", "mod" = "goldenrod"),
                     labels  = c("control" = "Control", "high" = "High", "mod" = "Moderate")) +
  labs(
    x     = "Measurement",
    y     = expression(A[max]~"("*mu*"mol CO"[2]~m^{-2}~s^{-1}*")"),
    color = "Elevation Treatment"
  ) +
  theme_bw() +
  theme(legend.position = "top")

