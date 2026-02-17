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
## BMB: shouldn't slice by row number, should filter based on
## information in data set

##Pivot longer and make strain names unique and change time names
hemolysin_long <- hemolysin |>
  ## BMB: modify rename so it doesn't get rid of C in TriazoleClass
  rename_with(~ str_replace_all(.x, "([.]|ZOC.CD)", "")) |>
  mutate(Strain = make.unique(as.character(Strain))) |> #Creates unique names
  pivot_longer(cols = contains("h"), names_to = "Time", values_to = "ZOC_CD")

##Model
hemo_mod <- lmer(
  ## BMB: why are you using log1p? You don't have any zero values here.
  ## log() would be more appropriate, easier to interpret
  log1p(ZOC_CD) ~ TriazoleClass * Time + (1 | Strain),
  data = hemolysin_long)
## BMB: note that (1 | Strain) is probably too specific, may need
## to allow (Time | Strain) or (Triazole
## if these are counts we may want to use a GLM instead?
summary(hemo_mod)

hemo_mod2 <- update(hemo_mod, log(ZOC_CD) ~ .)
## changes very little because data don't vary much, proportionally

##Extract EMMs per time and 95% CIs
emm_hemo <- emmeans(
  hemo_mod,
  ~ TriazoleClass | Time)

## BMB: did you ever show a prediction plot?
plot(emm_hemo)

## BMB: why convert to data frame?
emm_hemo_df <- as.data.frame(emm_hemo)

##Pairwise differences between Triazole classes (WITH 95% CI)
hemo_contrasts <- contrast(
  emm_hemo,
  method = "pairwise", ## BMB: only one pair here anyway???
  adjust = "tukey")

hemo_contrasts_df <- confint(hemo_contrasts, level = 0.95) |>
  as.data.frame()

##Plot residuals
residuals_hemo <- plot(hemo_mod)
qqplot_hemo <- qqnorm(resid(hemo_mod))
qqplot_hemo <- qqline(resid(hemo_mod))
## BMB: note these have a heavy left tail, thin right tail

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

## BMB: mark 2
