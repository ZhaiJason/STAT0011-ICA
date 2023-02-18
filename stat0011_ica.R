# Data =========================================================================
SP500.temp <- read.csv("data/SP500.csv")
SP500.price <- as.numeric(gsub(",", "", SP500.temp$Price))

CAC40.temp <- read.csv("data/CAC40.csv")
CAC40.price <- CAC40.temp$Close

# Log return
SP500 <- log(tail(SP500.price, -1)) - log(head(SP500.price, -1))
CAC40 <- log(tail(CAC40.price, -1)) - log(head(CAC40.price, -1))

rm(SP500.temp, CAC40.temp, SP500.price, CAC40.price)

# Plotting =====================================================================

# Boundaries for plot
ylim <- range(SP500, CAC40)

# Price trend plot
plot(SP500, type = "l", ylim = ylim, col = "blue",
     main = "Price History", ylab = "Price")
lines(CAC40, type = "l", col = "red")

# Scatter plot
plot(SP500, CAC40, pch = 20, col = rgb(0,0,0,0.2))

# Histograms
hist(SP500, breaks = 20, freq = FALSE)
hist(CAC40, breaks = 20, freq = FALSE)

# PIT ==========================================================================
library(MASS)
library(KScorrect)
library(fGarch)

# Standard t
SP500.params.t <- fitdistr(SP500, "t")$estimate
CAC40.params.t <- fitdistr(CAC40, "t")$estimate

hist(SP500, breaks = 20, freq = FALSE)
x <- seq(-0.4, 0.4, 0.001)
lines(x, dstd(x, SP500.params.t[1], SP500.params.t[2]*1.6, SP500.params.t[3]), type = "l")

hist(CAC40, breaks = 20, freq = FALSE)
x <- seq(-0.4, 0.4, 0.001)
lines(x, dstd(x, CAC40.params.t[1], CAC40.params.t[2]*1.2, CAC40.params.t[3]), type = "l")

SP500.u.t <- pstd(SP500, SP500.params.t[1], SP500.params.t[2]*1.6, SP500.params.t[3])
CAC40.u.t <- pstd(CAC40, CAC40.params.t[1], CAC40.params.t[2]*1.2, CAC40.params.t[3])

hist(SP500.u.t)
hist(CAC40.u.t)

LcKS(SP500.u.t, cdf = "punif")$p.value
LcKS(CAC40.u.t, cdf = "punif")$p.value

plot(SP500.u.t, CAC40.u.t)
model = VineCopula::BiCopSelect(SP500.u.t, CAC40.u.t, familyset = NA, selectioncrit = "AIC", indeptest = TRUE, level = 0.05)
model

# Packages =====================================================================
# Implementation of the copula analysis.
# install.packages("VineCopula")
# library(VineCopula)

# Implementation of the Anderson-Darling Goodness-of-Fit test.
# install.packages("goftest")
# library(goftest)

# Implementation of the (Lilliefors-Corrected) Kolmogorov-Smirnov Goodness-of-Fit test.
# install.packages("KScorrect")
# library(KScorrect)

# Implementation of the time series analysis.
# install.packages("fGarch")
# library(fGarch)
