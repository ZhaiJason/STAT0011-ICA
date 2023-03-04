# library(tidyverse)

SP500 <- read.csv("data/SP500.csv")
NK225 <- read.csv("data/NK225.csv")

NK225.close <- NK225$Price
SP500.close <- SP500$Price

SP500.close <- as.numeric(gsub(",", "", SP500.close))

# The below two lines will be numeric 0's
NK225.lr <- na.omit(log(NK225.close / stats::lag(NK225.close)))
SP500.lr <- na.omit(log(SP500.close / stats::lag(SP500.close)))

# This two rows will provide meaningful log returns
NK225.lr_2 <- na.omit(log(NK225.close / dplyr::lag(NK225.close)))
SP500.lr_2 <- na.omit(log(SP500.close / dplyr::lag(SP500.close)))
