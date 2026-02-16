library(tidyverse)
library(emmeans)
library(lmerTest)
library(ggplot2)
library(patchwork)

##Load RDS file (script from assignment 2 with cleaned data)
AllAf_data_clean <- read_rds("AllAf_data_clean.rds")
str(AllAf_data_clean)

##Select for triazole class and hemolysin values
hemolysin <- AllAf_data_clean|> 
  select(Strain, Generation, TriazoleClass, ZOC.CD.24h., ZOC.CD.48h., ZOC.CD.72h.) |>
  slice(326:437)

##Pivot longer and make strain names unique and change time names
hemolysin_long <- hemolysin |>
  rename_with(~ str_replace_all(.x, "[ZOC.CD.;]", "")) |>
  mutate(Strain = make.unique(as.character(Strain))) |> #Creates unique names
  pivot_longer(cols = contains("h"), names_to = "Time", values_to = "ZOC_CD")

##Model
hemo_mod <- lmer(
  log1p(ZOC_CD) ~ Triazolelass * Time + (1 | Strain),
  data = hemolysin_long)

summary(hemo_mod)

##Extract EMMs per time and 95% CIs
emm_hemo <- emmeans(
  hemo_mod,
  ~ Triazolelass | Time)

emm_hemo_df <- as.data.frame(emm_hemo)

##Pairwise differences between Triazole classes (WITH 95% CI)
hemo_contrasts <- contrast(
  emm_hemo,
  method = "pairwise",
  adjust = "tukey")

hemo_contrasts_df <- confint(hemo_contrasts, level = 0.95) |>
  as.data.frame()

##Plot residuals
residuals_hemo <- plot(hemo_mod)
qqplot_hemo <- qqnorm(resid(hemo_mod))
qqplot_hemo <- qqline(resid(hemo_mod))

emm_hemo_df <- as.data.frame(emm_hemo)

##Plot showing 95% CI for difference in group means (crossing 0 means no significant difference between resistant and susceptible groups)
ggplot(hemo_contrasts_df,
       aes(x = Time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2,
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Time",
    y = "Difference in hemolysin (Resistant − Susceptible, log scale)")
