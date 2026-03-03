library(tidyverse)
library(emmeans)
library(lmerTest)

##Load RDS file
AllAf_data_clean <- read_rds("AllAf_data_clean.rds")

##RQ: Do virulence traits (biofilm formation, hemolysin production) and thermal growth predict triazole resistance in Aspergillus fumigatus strains?
summary(AllAf_data_clean)
AllAf_data_clean$TriazoleClass <- factor(AllAf_data_clean$TriazoleClass)
AllAf_data_clean$TriazoleClass <- relevel(AllAf_data_clean$TriazoleClass, ref = "susceptible")
table(AllAf_data_clean$TriazoleClass) #check counts

##Predictor variables are TempGrowth41_48h, Standardized Biofilm formation, and Hemolysin Production after 48h
AllAf_data_clean <- AllAf_data_clean |>
  rename(TempGrowth41_48h = `TempGrowth41.48h.`) |>
  rename(ZOC_CD_48h = `ZOC.CD.48h.`)

##Fitting logistic GLM
glm_resistance <- glm(
  TriazoleClass ~ TempGrowth41_48h +
    StandardizedBiofilm_OD +
    ZOC_CD_48h,
  data = AllAf_data_clean,
  family = binomial)

summary(glm_resistance)

##tells me that a one-unit increase in Z0C_CD_48h decreases odds of resistance by 98% and is statistically significant. No evidence that thermotolerance or biofilm formation predicts resistance.
exp(coef(glm_resistance))

par(mfrow=c(2,2))
plot(glm_resistance)

##These 4 plots show us different things:
##Residuals vs Fitted: Most points clustered around 0 with some residuals as high as ~10-15 indicates poor prediction among some strains
##Q-Q Residuals: Pretty strong deviation at the upper tail likely the same strains that have large residuals
##Scale-Location: Increasing spread of residuals with fitted values. Upward trend in red smoothing line. This shows variance changes across fitted probabilities.
##Residuals vs Leverage: Some large residuals showing some influence among those strains.

##Add interaction between virulence and thermal growth
glm_resistance_int <- glm(TriazoleClass ~ TempGrowth41_48h * StandardizedBiofilm_OD * ZOC_CD_48h,
    family = binomial,
    data = AllAf_data_clean)

summary(glm_resistance_int)

par(mfrow=c(2,2))
plot(glm_resistance_int)
##Similar to what we see in previous model

##Overall, there seems to be no interaction effect between thermal growth and virulence-related factors in predictin resistance. These predictors overall may not strongly predict resistance (apart from hemolysin production). I believe this is probably more due to biological noise rather than the actual model used. Another thought I had is that resistance is rare among this dataset so they could account for the extreme residuals??