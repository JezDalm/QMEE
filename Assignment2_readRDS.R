library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
library(stats)
library(stringr)
library(ggplot2)
library(lme4)
library(emmeans)
library(lmerTest)
library(pbkrtest)

##Load RDS file
AllAf_data_clean <- read_rds("AllAf_data_clean.rds")

##1. Now let's make a histogram to look at the distributions of Biofilm formation among each cross
print(ggplot(AllAf_data_clean, aes(x=StandardizedBiofilm_OD))
      + geom_histogram(binwidth = .02)
      + labs(title = "All Cross Biofilms")
) #shows us all crosses combined (not what we want)

##Let's organize this by specific cross to see how the distributions are different
subset_progeny <- AllAf_data_clean|> 
  filter(Generation == "Progeny")|> 
  select(Strain, Generation, StandardizedBiofilm_OD) #filter by progeny and keep columns of interest

by_cross_biofilm <- subset_progeny |> 
  split(subset_progeny$Strain) #split by cross

for (i in seq_along(by_cross_biofilm)) { #creates sequence of numbers 1-5
  current_cross <- by_cross_biofilm[[i]] #if i is 3, takes by_cross_biofilm dataframe
  cross_name <- names(by_cross_biofilm)[i] #Extract by_cross_biofilm and its name
  
biofilm_plot <- ggplot(current_cross, aes(x = StandardizedBiofilm_OD)) + 
  geom_histogram(fill = "steelblue", color = "white") + 
  labs(title = paste("Biofilm Histogram for", cross_name)) + 
  theme_minimal()
  
  print(biofilm_plot) 
}

##2. I want to compare sub-MIC means between Strains that are resistant and strains that are susceptible while factoring in MIC values of itraconazole and voriconazole
progeny_Itra <- AllAf_data_clean|> 
  select(Strain, Generation, TriazoleClass, contains("Itraconazole_")) #selects for columns that only have itraconazole data

progeny_Vori <- AllAf_data_clean|> 
  select(Strain, Generation, TriazoleClass, contains("Voriconazole_")) #selects for columns that only have voriconazole data

print(progeny_Itra)
print(progeny_Vori)

progeny_Itra_rename <- progeny_Itra |>
  rename_with(
    ~ paste0("", str_extract(.x, "[0-9.]+")),
    starts_with("Itraconazole_OD_")
  ) #rename column headers for different drug concentrations to numeric values

progeny_Itra_long <- progeny_Itra_rename |>
  pivot_longer(cols = matches("^[0-9.]+$"), names_to = "Concentration", values_to = "OD") |>
  mutate(Concentration = as.numeric(Concentration)) |>
  filter(Concentration <= 0.12) #pivot longer and only include concentrations < 0.25

##Same for voriconazole
progeny_Vori_rename <- progeny_Vori |>
    rename_with(
    ~ paste0("", str_extract(.x, "[0-9.]+")),
    starts_with("Voriconazole_OD_")
  )

progeny_Vori_long <- progeny_Vori_rename |>
  pivot_longer(cols = matches("^[0-9.]+$"), names_to = "Concentration", values_to = "OD") |>
  mutate(Concentration = as.numeric(Concentration)) |>
  filter(Concentration <= 0.12)

##Fitting a model for itraconazole and seeing which one is better
model_itr <- lmer(
  OD ~ TriazoleClass + Concentration + Itraconazole_MIC + (1 | Strain),
  data = progeny_Itra_long, 
  REML = FALSE
)

Model_itr_2 <- lmer(
  OD ~ TriazoleClass + Concentration + Itraconazole_MIC +
    (1 + Concentration | Strain),
  data = progeny_Itra_long, 
  REML = FALSE
)

anova(model_itr, Model_itr_2) #second model is better since adding a random slope for concentration significantly improves the model (Chisq = 28.797, p = 5.58 × 10^-7)

##Get EMMs for plotting
emm_conc_itra <- emmeans(
  Model_itr_2,
  ~ TriazoleClass | Concentration,
  at = list(
    Concentration = sort(unique(progeny_Itra_long$Concentration))
  )
)

emm_df_itra <- as.data.frame(emm_conc_itra) #convert to df
print(emm_df_itra)

emm_df_itra <- emm_df_itra %>%
  mutate(Concentration_f = factor(Concentration)) #convert Concentration to factor for plot

#Plot model
ggplot(emm_df_itra,
       aes(x = Concentration_f,
           y = emmean,
           color = TriazoleClass,
           fill  = TriazoleClass,
           group = TriazoleClass)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = lower.CL, ymax = upper.CL),
    alpha = 0.2,
    color = NA
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression("Itraconazole Concentration ("*mu*"g/mL)"),
    y = "Estimated marginal mean growth (OD)"
  ) +
  theme(
    legend.title = element_blank()
  )