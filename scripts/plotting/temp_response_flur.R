library(tidyverse)
library(lubridate)
library(dplyr)
library(car)

flu <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_full_final.csv")
curve_dat <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/response_curve_data.csv")

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

flu$ID <- paste(flu$Seedling_ID, flu$ID2)

head(flu)

flu_sum <- flu %>%
  group_by(Needle_age, Elev, Micro, Measurement)%>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
  summarise(
    n=n(),
    meanFvFm=mean(Fv.Fm, na.rm=TRUE),
    sdFvFm=sd(Fv.Fm, na.rm=TRUE),
    minFvFm=min(Fv.Fm, na.rm=TRUE),
    maxFvFm=max(Fv.Fm, na.rm=TRUE),
    .groups="drop")

head(flu_sum)

sum_dat <- curve_dat %>%
  group_by(Needle_age, Elev, Micro) %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
  summarize(
    n=n(),
    ctmax=mean(ctmax, na.rm=TRUE),
    tcrit=mean(tcrit, na.rm=TRUE),
    t15=mean(t15, na.rm=TRUE),
    .groups="drop"
  )

fit_curves_temp <- flu_sum %>%
  group_by(Needle_age, Elev, Micro) %>%
  group_modify(~ {
    dat <- .x
    cof <- tryCatch(
      coef(lm(car::logit(meanFvFm) ~ Measurement, data = dat)),
      error = function(e) c(0, 0)
    )
    m <- tryCatch(
      nls(meanFvFm ~ theta1 / (1 + exp(-(theta2 + theta3 * Measurement))),
          data = dat,
          start = list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
          control = nls.control(maxiter = 1000, tol = 1e-3)),
      error = function(e) NULL
    )
    if (is.null(m)) return(data.frame(Measurement = NA_real_, pred = NA_real_))
    temp_seq <- seq(min(dat$Measurement), max(dat$Measurement), length.out = 100)
    data.frame(Measurement = temp_seq, pred = predict(m, newdata = data.frame(Measurement = temp_seq)))
  }) %>%
  ungroup()


ref_lines <- sum_dat %>%
  select(Needle_age, Elev, Micro,tcrit, t15, ctmax) %>%
  pivot_longer(cols = c(tcrit, t15, ctmax),
                      names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("tcrit", "t15", "ctmax"), labels = c("Tcrit","T15","T50"))
)

figure <- ggplot() +
  geom_point(data = flu_sum, aes(x = Measurement, y = meanFvFm),
             color = "black", alpha = 0.7, size = 2) +
  geom_line(data = fit_curves_temp, aes(x = Measurement, y = pred),
            color = "gray34", linewidth = 1) +
  geom_vline(data = ref_lines, aes(xintercept = value, color = metric, linetype = metric),
             linewidth = 0.7) +
  facet_grid(Needle_age ~ Elev + Micro,
             labeller = labeller(label_value)) +
  scale_linetype_manual(values = c("Tcrit" = "dashed", "T15" = "dashed", "T50" = "dashed"),
                        guide = "none") +
  geom_hline(yintercept=0.4, linetype="dashed",color="gray34")+
  scale_color_manual(name = "Reference value", values = c("Tcrit" = "steelblue", "T15" = "orange", "T50" = "red")) +
  labs(x = "Temperature (°C)", y = expression(F[v]/F[m])) +
  ylim(0,0.83)+
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
figure
ggsave("temp_response_fig.png", figure, path = "~/Masters/Growth_Chamber/results/figures")
