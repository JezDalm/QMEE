
library(tidyverse) ## BMB: tidyverse loads all of the following packages
## library(dplyr)
## library(tidyr)
## library(readr)
## library(stringr)
## library(ggplot2)
## library(stats) ## BMB don't need to load this, it's automatic

AllAf_data <- read.csv("Af_CrossData.csv")
str(AllAf_data)

##First we want to convert MIC columns to characters just like Assignment 1
AllAf_data$Voriconazole_MIC <- as.character(AllAf_data$Voriconazole_MIC)
AllAf_data$AmphotericinB_MIC <- as.character(AllAf_data$AmphotericinB_MIC)

str(AllAf_data) #double check to see if columns were converted to characters

## BMB: could use mutate(across(...)); why not factors???

##Since this data includes all cross data, we want a column to assign each Strain to a generation and designate them as either 'Progeny' or 'Parent'. From this, we can also see if the count of parents adds up to 10 (5 biparental crosses) and if the progeny names match those of the parents.
print(AllAf_data
      |> count(Strain, name="count")
      |> arrange(desc(count))
) #gives us the frequency of different strain names so we know approximately how many parents and progeny strains there are

## BMB for only two possibilities you can use ifelse() or if_else()

AllAf_data_type <- AllAf_data |> 
  mutate(
    Generation = case_when(
      str_detect(Strain, "x") ~ "Progeny",
      TRUE ~ "Parent"
    )
  ) #creates new column called "Generation" where character "x" determines if the Strain is a "Progeny" and if condition not met, assign "Parent"

print(AllAf_data_type
      |> count(Strain, Generation, name = "count")
      |> arrange(desc(count))
      )

saveRDS(AllAf_data_type, "AllAf_data_clean.rds")
