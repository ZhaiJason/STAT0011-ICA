# Library ======================================================================
# library(MASS) # MLE distribution fit
# library(goftest)
# library(VineCopula) # For selecting copulas and MC simulation
# library(KScorrect) # LcKS test
# library(ADGofTest)
# library(fBasics)
library(fGarch)

# Function =====================================================================
plot.ac <- function(x) {
    par(mfrow = c(2, 2))
    acf(x)
    pacf(x)
    acf(x^2)
    par(mfrow = c(1, 1))
}

# Data =========================================================================
SP500.temp <- read.csv("data/SP500.csv")
SP500.price <- as.numeric(gsub(",", "", SP500.temp$Price))
SP500.price <- rev(SP500.price)

NK225.temp <- read.csv("data/NK225.csv")
NK225.price <- NK225.temp$Price
NK225.price <- rev(NK225.price)

# Log return
SP500 <- log(tail(SP500.price, -1)) - log(head(SP500.price, -1))
NK225 <- log(tail(NK225.price, -1)) - log(head(NK225.price, -1))

rm(SP500.temp, SP500.price, NK225.temp, NK225.price)

# Require GARCH
plot(SP500, type = "l", col = "red")
lines(NK225, col = "blue")

# Check for Autocorrelation ====================================================

Box.test(SP500, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Significant
plot.ac(SP500) # Cut-off at order 1

Box.test(NK225, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(NK225)

# TS Modelling =================================================================

# SP500.model <- garchFit(formula = ~arma(1,0) + garch(1,1), data = SP500, trace = FALSE, cond.dist = "norm")
SP500.model <- garchFit(formula = ~arma(1,0) + garch(1,1), data = SP500, trace = FALSE, cond.dist = "std")
SP500.res <- residuals(SP500.model, standardize = TRUE) # WHY DO WE NEED STANDARIZZED

Box.test(SP500.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
Box.test(SP500.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(SP500.res)
SP500.model@fit$ics

# PIT ==========================================================================

fBasics::jarqueberaTest(SP500)
fBasics::jarqueberaTest(NK225)



hist(SP500.res, breaks = 20, freq = FALSE)
x <- seq(-4, 4, 0.001)
lines(x, dstd(x, SP500.params[1], SP500.params[2]*1.2, SP500.params[3]), type = "l")


# MLE estimate, suggest using bbmle package
NK225.params <- MASS::fitdistr(NK225, "t")$estimate

hist(NK225, breaks = 20, freq = FALSE)
x <- seq(-0.4, 0.4, 0.001)
lines(x, dstd(x, NK225.params[1], NK225.params[2]*1.2, NK225.params[3]), type = "l")

rm(x)

# SP500.u <- pstd(SP500.res, SP500.params[1], SP500.params[2]*1.2, SP500.params[3])
SP500.u <- pstd(SP500.res, 0, 1, SP500.model@fit$coef["shape"])
NK225.u <- pstd(NK225, NK225.params[1], NK225.params[2]*1.2, NK225.params[3])

hist(SP500.u)
hist(NK225.u)

KScorrect::LcKS(SP500.u, cdf = "punif")$p.value
KScorrect::LcKS(NK225.u, cdf = "punif")$p.value
ADGofTest::ad.test(SP500.u, null = "punif")$p.value
ADGofTest::ad.test(NK225.u, null = "punif")$p.value

plot(SP500.u, NK225.u, pch = 20) # Preferrable combination, grumble copula maybe

# Construct Copula =============================================================

model <- VineCopula::BiCopSelect(SP500.u, NK225.u, familyset = NA, selectioncrit = "AIC", indeptest = TRUE, level = 0.05)
model

n <- 2000
sim.u <- VineCopula::BiCopSim(n, family = 2, model$par, model$par2)
plot(sim.u, pch = 20)

# SP500.sim <- qstd(sim.u[, 1], SP500.params[1], SP500.params[2]*1.2, SP500.params[3])
SP500.sim <- qstd(sim.u[, 1], 0, 1, SP500.model@fit$coef["shape"])
NK225.sim <- qstd(sim.u[, 2], NK225.params[1], NK225.params[2]*1.2, NK225.params[3])

# Re-introduce Autocorrelation & GARCH Effects =================================

SP500.coef <- SP500.model@fit$coef

# Introduce GARCH
mu <- SP500.coef["mu"]
ar1 <- SP500.coef["ar1"]
omega <- SP500.coef["omega"]
alpha1 <- SP500.coef["alpha1"]
beta1 <- SP500.coef["beta1"]
sigma2s <- numeric(n)
y <- numeric(n)
sigma2s[1] <- 0
y[1] <- mu + SP500.sim[1]

for (i in 2:n) {
    sigma2s[i] <- omega + alpha1 * sigma2s[i - 1] * SP500.sim[i - 1]^2 + beta1 * sigma2s[i - 1]
    y[i] <- mu + ar1 * y[i - 1] + sqrt(sigma2s[i]) * SP500.sim[i]
}

plot(y, type = "l")

# par(mfrow = c(1, 2))
plot(SP500, NK225, xlim = c(-0.5, 0.5), ylim = c(-0.3, 0.3))
# plot(y, NK225.sim, xlim = c(-0.5, 0.5), ylim = c(-0.3, 0.3))
points(y, NK225.sim, col = "blue")

# Compute VaR ==================================================================

portsim <- matrix(0, nrow = n, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

portsim <- log(1 + ((exp(SP500.sim) - 1) + (exp(NK225.sim) - 1)) * (1 / 2))
varsim <- quantile(portsim, c(0.01,0.05))
varsim
