library(ggplot2)
library(tidyverse)
library(ggbeeswarm)
library(Hmisc)

pdf(width=10)

## JD: I really like that you are using a previously made .rds. Next time also mention that this script depends on your clearning script.

##Load RDS file
AllAf_data_clean <- read_rds("AllAf_data_clean.rds")
str(AllAf_data_clean)

##From Assignment 2, it was recommended I also present the raw data for the OD values I have at each drug concentration across the three drugs I tested. For the sake of this assignment, I will only select one cross to analyze.
progeny_OD <- AllAf_data_clean|> 
  select(Strain, Generation, TriazoleClass, contains("_OD_")) |>
  slice(326:437) #selects for columns that only have antifungal OD data and the specific rows that I want to analyze for one cross

progeny_Itra <- progeny_OD|> 
  select(Strain, Generation, TriazoleClass, contains("Itraconazole_")) #selects for columns that only have itraconazole data

progeny_Vori <- progeny_OD|> 
  select(Strain, Generation, TriazoleClass, contains("Voriconazole_")) #selects for columns that only have voriconazole data

progeny_AMB <- progeny_OD|> 
  select(Strain, Generation, TriazoleClass, contains("AmphotericinB_")) #selects for columns that only have AMB data

##I tried to make a function for above so I don't have to do it seperate each time for each drug but was struggling to figure it out.
# progeny_combined_long <- progeny_OD |>
#   pivot_longer(
#     cols = contains(c("Itraconazole_", "Voriconazole_", "AmphotericinB_")),
#     names_to = c("Drug", "Concentration"),
#     names_sep = "_",
#     values_to = "OD"
#   ) |>
#   mutate(Concentration = as.numeric(Concentration))

##Making a function to change column names for all drugs and pivot longer for easy readability by ggplot. Also change Concentration to a factor
drug_list <- list(
  Itraconazole = progeny_Itra,
  Voriconazole = progeny_Vori,
  AmphotericinB = progeny_AMB
)

progeny_combined_long <- drug_list |> 
  map_df(~ {
    .x |> 
      rename_with(~ str_extract(.x, "[0-9.]+"), matches("OD_")) |> 
      pivot_longer(cols = matches("^[0-9.]+$"), names_to = "Concentration", values_to = "OD",
      names_transform = list(Concentration = as.factor)
      )}, .id = "Drug")

##reorder based on numeric values
progeny_combined_long$Concentration <- fct_inseq(progeny_combined_long$Concentration)

##Beeswarm plot of my raw antifungal OD data with geom_smooth. I commented geom_boxplot() and geom_smooth() out because I think it looks best without them but feel free to provide any feedback about this. I mainly want to use this plot to demonstrate the quantitative variation when using OD to measure drug susceptibilities rather than qualitative MIC values.
ggplot(progeny_combined_long, aes(x = Concentration, y = OD, color = Drug, group = interaction(Drug, Concentration))) +
  #geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2), # sets 2 SD
               geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75), color = "black") +
  #geom_smooth(aes(group = Drug), method = "loess", se = FALSE, size = 2) +
  geom_quasirandom(alpha = 0.4, dodge.width = 0.75, bandwidth = 0.1) +
  coord_cartesian(ylim = c(-.1, 0.6)) + 
  xlab("Drug Concentration in µg/ml") +
  ylab("Mean OD at 530nm") +
  geom_point(aes(fill = Drug, shape = Strain),
             subset(progeny_combined_long, Strain %in% c("I162", "A3-4")),
             colour = "black",
             position = position_dodge(width = .7),
             size = 3) +
  scale_color_manual(values = c("mediumorchid1", "aquamarine", "#E7B800"),
                     breaks = c("Voriconazole", "AmphotericinB", "Itraconazole")) +
  scale_shape_manual(name = "Parents", values = c("I162" = 17, "A3-4" = 15)) +
  guides(fill = "none") +
  theme_bw(base_size = 13)

##Now let's make reaction norm plots for our temperature data within the same cross and compare triazole resistant versus susceptible strains
progeny_temp <- AllAf_data_clean|> 
  select(Strain, Generation, TriazoleClass, contains("Temp")) |>
  slice(326:437) #selects for columns that only have Temperature OD data and the specific rows that I want to analyze for one cross

##Change duplicate Strain names to unique so I can plot individual lines for each strain.
progeny_temp <- progeny_temp |>
  mutate(Strain = make.unique(as.character(Strain))) #Creates unique names

##pivot longer so that original columns are split into "Temperature" and "Time" and change those columns to a factor.
progeny_temp_long <- progeny_temp |>
  select(Strain, Generation, TriazoleClass, contains("Temp")) |>
  pivot_longer(
    cols = contains("Temp"),
    # Define the two new columns to create
    names_to = c("Temperature", "Time"), 
    names_pattern = "TempGrowth(\\d+).(\\d+)h", 
    values_to = "OD"
  ) |>
  mutate(
    Temperature = as.factor(Temperature),
    Time = as.factor(Time)
  )

## JD: I think the plot would be helped by transparency, and also by dodging. 
## I added the transparency for now; it's not completely obvious to me how to dodge with both lines and points, but I guess not too hard.
## group = Strain would be more beautifully moved to the plot-level aes
## I also added some code at the top that makes both plots wider (and plots them as pdfs, probably you can find a nicer way with rstudio). 
### The second plot definitely looks nicer when it's wide.
##This plot is not ideal and not really I want to present this data. Any feedback on how I could improve this to show if resistant versus susceptible strains have reaction norm differences to different temperatures? This might not be possible because of just how the data is.
reaction_norm <- ggplot(progeny_temp_long, aes(x = Temperature, y = OD, color = TriazoleClass)) +
  geom_line(aes(group = Strain), alpha = 0.4, linewidth = 1) +
  geom_point(aes(group = Strain), size = 2) +
  labs(x = "Temperature (°C)", y = "Mean Optical Density (OD)", color= "Triazole Phenotype") +
  theme_bw(base_size = 14) +
  scale_x_discrete() +
  scale_y_continuous(limits = c(0, 0.6)) + #manually removed an extreme point per day for one of the strains. Removing limits shows the outlier and makes it hard to tell the pattern of the different groups
  facet_wrap(~Time) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  guides(linetype = "none") +
  scale_color_manual(values = c("susceptible" = "#0072B2",
                                "resistant" = "#D81B60"))

plot(reaction_norm)

## JD: The move of just removing a point seems a bit weird. Also, it would be good to indicate that the facet variable is time, maybe with a paste command.

##I attempted to transform it but didn't really look better
#reaction_norm + scale_y_log10()

## Grade 2.1/3

## PS: For dodging, you might need to use something like group = interaction(TriazoleClass, Strain). Ask us on Teams if you have issues.
