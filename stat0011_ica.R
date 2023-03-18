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
    par(mfrow = c(1, 3))
    acf(x)
    pacf(x)
    acf(x^2)
    par(mfrow = c(1, 1))
}

# Data =========================================================================

# Read-in csv data of stock prices
SP500.temp <- read.csv("SP500.csv", encoding = "UTF-8")
SP500.price <- as.numeric(gsub(",", "", SP500.temp$Price))
SP500.price <- rev(SP500.price)
NK225.temp <- read.csv("NK225.csv", encoding = "UTF-8")
NK225.price <- NK225.temp$Price
NK225.price <- rev(NK225.price)

# Calculate log-return
SP500 <- log(tail(SP500.price, -1)) - log(head(SP500.price, -1))
NK225 <- log(tail(NK225.price, -1)) - log(head(NK225.price, -1))

# Remove temporary variables
rm(SP500.temp, SP500.price, NK225.temp, NK225.price)

# Check for Autocorrelation ====================================================

Box.test(SP500, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Significant
Box.test(NK225, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed

plot.ac(SP500) # Cut-off at order 1
plot.ac(NK225)

# TS Modelling =================================================================

# SP500.model <- garchFit(formula = ~arma(1,0) + garch(1,1), data = SP500, trace = FALSE, cond.dist = "norm")
SP500.model <- garchFit(formula = ~arma(1, 0) + garch(1, 1), data = SP500, trace = FALSE, cond.dist = "std")
SP500.res <- residuals(SP500.model, standardize = TRUE) # WHY DO WE NEED STANDARIZZED
SP500.model@fit$ics
SP500.coef <- SP500.model@fit$coef

Box.test(SP500.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
Box.test(SP500.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(SP500.res)

# Test GARCH for NK225
# NK225.model <- garchFit(formula = ~garch(1, 1), data = NK225, trace = FALSE, cond.dist = "std")
# NK225.res <- residuals(NK225.model, standardize = TRUE)
#
# Box.test(NK225.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
# Box.test(NK225.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
# plot.ac(NK225.res)
#
# plot(NK225, type = "l")
# plot(NK225.res, type = "l")

# PIT ==========================================================================

# Check if data from normal distribution
fBasics::jarqueberaTest(SP500.res, description = "Test of Normality") # Significant
fBasics::jarqueberaTest(NK225, description = NA) # Significant

par(mfrow = c(1, 2))

# Approximate suitable distribution for SP500 model residuals' PIT
SP500.params <- MASS::fitdistr(SP500.res, "t")$estimate

hist(SP500.res, breaks = 20, freq = FALSE)
x <- seq(-4, 4, 0.001)
lines(x, fGarch::dstd(x, SP500.params[1], SP500.params[2]*1.2, SP500.params[3])) # Parameters 0 and 1 to match the conditional distribution in garchFit
# lines(x, fGarch::dstd(x, 0, 1, SP500.coef["shape"]), type = "l") # Parameters 0 and 1 to match the conditional distribution in garchFit

# Perform MLE estimate on NK225
NK225.params <- MASS::fitdistr(NK225, "t")$estimate

hist(NK225, breaks = 20, freq = FALSE)
x <- seq(-0.4, 0.4, 0.001)
lines(x, fGarch::dstd(x, NK225.params[1], NK225.params[2]*1.2, NK225.params[3])) # Parameters manually adjusted to enhance distribution fit

par(mfrow = c(1, 1))

# Remove temporary variable
rm(x)

SP500.u <- pstd(SP500.res, SP500.params[1], SP500.params[2]*1.2, SP500.params[3])[-1]
# SP500.u <- pstd(SP500.res, 0, 1, SP500.coef["shape"])[-1]
NK225.u <- pstd(NK225, NK225.params[1], NK225.params[2]*1.2, NK225.params[3])[-1]

hist(SP500.u, breaks = 20)
hist(NK225.u, breaks = 20)

KScorrect::LcKS(SP500.u, cdf = "punif")$p.value
KScorrect::LcKS(NK225.u, cdf = "punif")$p.value
ADGofTest::ad.test(SP500.u, null = "punif")$p.value
ADGofTest::ad.test(NK225.u, null = "punif")$p.value

plot(SP500.u, NK225.u, pch = 20) # Preferrable combination, grumble copula maybe

# Construct Copula =============================================================

model <- VineCopula::BiCopSelect(SP500.u, NK225.u, familyset = NA, selectioncrit = "AIC", indeptest = TRUE, level = 0.05)
model

plot(SP500.u, NK225.u, pch = 20)

set.seed(10)

n <- 2000
sim.u <- VineCopula::BiCopSim(n, family = 2, model$par, model$par2)
plot(sim.u, pch = 20)

# SP500.sim <- qstd(sim.u[, 1], SP500.params[1], SP500.params[2]*1.2, SP500.params[3])
SP500.sim <- qstd(sim.u[, 1], SP500.params[1], SP500.params[2]*1.2, SP500.params[3])
# SP500.sim <- qstd(sim.u[, 1], 0, 1, SP500.model@fit$coef["shape"])
NK225.sim <- qstd(sim.u[, 2], NK225.params[1], NK225.params[2]*1.2, NK225.params[3])

# Re-introduce Autocorrelation & GARCH Effects =================================

plot(SP500.sim, NK225.sim)

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

SP500.res <- y
NK225.res <- NK225.sim

plot(SP500, NK225, pch = 20, col = rgb(0, 0, 0, 0.1), main = "Distribution of Log-returns")
points(SP500.res, NK225.res, pch = 20, col = rgb(1, 0, 0, 0.1))
legend("bottomright", pch = 20, col = rgb(c(0, 1), 0, 0, 0.5), legend = c("Original Data", "Simulated Data"), cex = 0.75, bty = "n")

# Compute VaR ==================================================================

portsim <- matrix(0, nrow = n, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

portsim <- log(1 + ((exp(SP500.sim) - 1) + (exp(NK225.sim) - 1)) * (1 / 2))
varsim <- quantile(portsim, c(0.01, 0.05))
varsim
