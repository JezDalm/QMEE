library(dplyr)
library(tidyr)
library(readr)
library(stats)
library(corrplot)

##Read Dataset
Af_data <- read.csv("Af_J57xP20.csv")

##Can see some important info about lots of quantitative data here
summary(Af_data)
str(Af_data)

##Do this to see the frequency counts of categorical variables
summary(Af_data |> mutate_if(is.character, as.factor))

##Here we notice that some MIC columns are not being treated as characters. We want to change them to characters so that their frequencies are read when using summary()
Af_data$Voriconazole_MIC <- as.character(Af_data$Voriconazole_MIC)
Af_data$AmphotericinB_MIC <- as.character(Af_data$AmphotericinB_MIC)
str(Af_data) #double check to see if columns were converted to characters

##Now this should work
(Af_data |> mutate(across(where(is.character), as.factor)) 
  |> summary()
)

##Let's perform a pairwise correlation between the different quantitative traits to see if there are any relationships. Here we specify we only want correlations between numerical values
Af_correlations <- cor(Af_data[sapply(Af_data, is.numeric) #apply to only numeric columns
  ]
)

##Now, let's plot our correlation matrix
corrplot(Af_correlations,
         method = "color", #shades in tile according to color gradient
         outline = TRUE, #adds black tile outline
         type = "upper",
         tl.col = "black", #text color
         tl.cex = 0.6) #text font size