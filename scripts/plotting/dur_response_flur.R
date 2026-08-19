library(tidyverse)
library(lubridate)
library(dplyr)
library(car)

# Read in data
dur_dat <- read.csv("~/Masters/Growth_Chamber/data/Chlorophyll_fluorescence/processed/duration_curve_data.csv")
flur2<- read.csv("~/Masters/Growth_chamber/data/Chlorophyll_fluorescence/raw/chl_fluor_treatment2_full.csv")

dur <- flur2 %>%
  filter(Measurement %in% c("0", "3", "4", "5", "6", "7")) %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control", "mod", "high"), labels = c("Control","Moderate","High"))) %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
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

dur$ID <- paste(dur$Seedling_ID, dur$ID2)
head(dur)

dur_sum <- dur %>%
  group_by(Needle_age, Prev_trmt, Elev, Micro, time)%>%
  summarise(
    n=n(),
    meanFvFm=mean(Fv.Fm, na.rm=TRUE),
    sdFvFm=sd(Fv.Fm, na.rm=TRUE),
    minFvFm=min(Fv.Fm, na.rm=TRUE),
    maxFvFm=max(Fv.Fm, na.rm=TRUE),
    .groups="drop")

head(dur_sum)

sum_dat <- dur_dat %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control", "mod", "high"), labels = c("control","moderate","high"))) %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  mutate(Micro = factor(Micro, levels = c("inoc","uninoc"), labels = c("Inoculated","Uninoculated"))) %>%
  group_by(Needle_age, Prev_trmt, Elev, Micro) %>%
  summarize(
    n=n(),
    D50=mean(D50, na.rm=TRUE),
    Dcrit=mean(Dcrit, na.rm=TRUE),
    D15=mean(D15, na.rm=TRUE),
    .groups="drop"
  )

fit_curves_dur <- dur_sum %>%
  group_by(Needle_age, Prev_trmt, Elev, Micro) %>%
  group_modify(~ {
    dat <- .x
    cof <- tryCatch(
      coef(lm(car::logit(meanFvFm) ~ time, data = dat)),
      error = function(e) c(0, 0)
    )
    m <- tryCatch(
      nls(meanFvFm ~ theta1 / (1 + exp(-(theta2 + theta3 * time))),
          data = dat,
          start = list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
          control = nls.control(maxiter = 1000, tol = 1e-3)),
      error = function(e) NULL
    )
    if (is.null(m)) return(data.frame(time = NA_real_, pred = NA_real_))
    time_seq <- seq(min(dat$time), max(dat$time + 1), length.out = 100)
    data.frame(time = time_seq, pred = predict(m, newdata = data.frame(time = time_seq)))
  }) %>%
  ungroup()

ref_lines <- sum_dat %>%
  select(Needle_age, Prev_trmt, Elev, Micro, Dcrit, D15, D50) %>%
  pivot_longer(cols = c(Dcrit, D15, D50),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Dcrit", "D15", "D50"))
  )

figure <- ggplot() +
  geom_point(data = dur_sum, aes(x = time, y = meanFvFm, color = Prev_trmt),
             size = 2) +
  geom_line(data = fit_curves_dur, aes(x = time, y = pred, color = Prev_trmt),
            linewidth = 1, alpha = 0.6) +
  facet_grid(Needle_age ~ Elev + Micro,
             labeller = labeller(label_value)) +
  scale_color_manual(values = c("Control" = "#21AF29", "Moderate" = "#D4AB7B", "High" = "#A7144C")) +
  labs(x = "Duration at 45°C (hours)", y = expression(F[v]/F[m]),
       color = "Previous Treatment") +
  geom_hline(yintercept = 0.4, lty = "dashed", color = "gray34") + 
  geom_hline(yintercept = 0.8, lty = "dashed", color = "gray70") +
  xlim(0,8) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom",
        panel.grid.minor=element_blank())

figure
ggsave("dur_response_fig.png", figure, path = "~/Masters/Growth_Chamber/results/figures")

# No microbe
dur_sum2 <- dur %>%
  group_by(Needle_age, Prev_trmt, Elev, time)%>%
  summarise(
    n=n(),
    meanFvFm=mean(Fv.Fm, na.rm=TRUE),
    sdFvFm=sd(Fv.Fm, na.rm=TRUE),
    minFvFm=min(Fv.Fm, na.rm=TRUE),
    maxFvFm=max(Fv.Fm, na.rm=TRUE),
    .groups="drop")


sum_dat2 <- dur_dat %>%
  mutate(Prev_trmt = factor(Prev_trmt, levels = c("control", "mod", "high"), labels = c("Control","Moderate","High"))) %>%
  mutate(Needle_age = factor(Needle_age, levels = c("new","old"), labels=c("Developing Foliage","Mature Foliage"))) %>%
  mutate(Elev = factor(Elev, levels = c("low","high"), labels = c("Low Elevation", "High Elevation"))) %>%
  group_by(Needle_age, Prev_trmt, Elev) %>%
  summarize(
    n=n(),
    D50=mean(D50, na.rm=TRUE),
    Dcrit=mean(Dcrit, na.rm=TRUE),
    D15=mean(D15, na.rm=TRUE),
    .groups="drop"
  )

fit_curves_dur2 <- dur_sum2 %>%
  group_by(Needle_age, Prev_trmt, Elev) %>%
  group_modify(~ {
    dat <- .x
    cof <- tryCatch(
      coef(lm(car::logit(meanFvFm) ~ time, data = dat)),
      error = function(e) c(0, 0)
    )
    m <- tryCatch(
      nls(meanFvFm ~ theta1 / (1 + exp(-(theta2 + theta3 * time))),
          data = dat,
          start = list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
          control = nls.control(maxiter = 1000, tol = 1e-3)),
      error = function(e) NULL
    )
    if (is.null(m)) return(data.frame(time = NA_real_, pred = NA_real_))
    time_seq <- seq(min(dat$time), max(dat$time + 3), length.out = 100)
    data.frame(time = time_seq, pred = predict(m, newdata = data.frame(time = time_seq)))
  }) %>%
  ungroup()

ref_lines2 <- sum_dat2 %>%
  select(Needle_age, Prev_trmt, Elev, Dcrit, D15, D50) %>%
  pivot_longer(cols = c(Dcrit, D15, D50),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Dcrit", "D15", "D50"))
  )

figure2 <- ggplot() +
  geom_point(data = dur_sum2, aes(x = time, y = meanFvFm),
             alpha = 0.6, size = 2) +
  geom_line(data = fit_curves_dur2, aes(x = time, y = pred, color = Prev_trmt),
            linewidth = 1.1, alpha =0.6) +
  geom_vline(data = ref_lines2,
             aes(xintercept = value, color = Prev_trmt, linetype = metric),
             linewidth = 0.7) +
  facet_grid(Prev_trmt ~ Elev + Needle_age,
             labeller = labeller(label_value)) +
  scale_color_manual(values = c("Control" = "#21AF29", "Moderate" = "#D4AB7B", "High" = "#A7144C")) +
  scale_linetype_manual(values = c("Dcrit" = "dashed", "D15" = "dotted", "D50" = "solid")) +
  labs(x = "Duration at 45°C (hours)", y = expression(F[v]/F[m]),
       color = "Previous Treatment", linetype = "Reference value") +
  geom_hline(yintercept = 0.4, lty = "dashed", color = "gray34") + 
  geom_hline(yintercept = 0.8, lty = "dashed", color = "gray70") +
  xlim(0,8.5) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom",
        panel.grid.minor=element_blank())

figure2
ggsave("dur_response_fig2.png", figure2, path = "~/Masters/Growth_Chamber/results/figures")

