plot.ac <- function(x) {
    par(mfrow = c(1, 3))
    acf(x)
    pacf(x)
    acf(x^2)
    par(mfrow = c(1, 1))
}

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

Box.test(NK225^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Significant

library(fGarch)

# SP500.model <- garchFit(formula = ~arma(1,0) + garch(1,1), data = SP500, trace = FALSE, cond.dist = "norm")
SP500.model <- garchFit(formula = ~arma(1, 0) + garch(1, 1), data = SP500, trace = FALSE, cond.dist = "std")
SP500.res <- residuals(SP500.model, standardize = TRUE) # WHY DO WE NEED STANDARIZZED
SP500.model@fit$ics
SP500.coef <- SP500.model@fit$coef

Box.test(SP500.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
Box.test(SP500.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(SP500.res)

# NK225 GARCH(1,1) model
NK225.model <- garchFit(formula = ~garch(1, 1), data = NK225, trace = FALSE, cond.dist = "std")
NK225.res <- residuals(NK225.model, standardize = TRUE) # WHY DO WE NEED STANDARIZZED
NK225.model@fit$ics
NK225.coef <- NK225.model@fit$coef

Box.test(NK225.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
Box.test(NK225.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(NK225.res)

# PIT ==========================================================================
fBasics::jarqueberaTest(SP500.res, description = "Test of Normality") # Significant
fBasics::jarqueberaTest(NK225.res, description = NA) # Significant

par(mfrow = c(1, 2))

hist(SP500.res, breaks = 20, freq = FALSE, ylim = c(0, 0.5))
x <- seq(-4, 4, 0.001)
lines(x, fGarch::dstd(x, 0, 1, SP500.coef["shape"]), type = "l", col = "red")
legend("topleft", lty = 1, col = "red", legend = "Fitted Distribution", cex = 0.75, bty = "n")

hist(NK225.res, breaks = 20, freq = FALSE, ylim = c(0, 0.5))
x <- seq(-4, 4, 0.001)
lines(x, fGarch::dstd(x, 0, 1, NK225.coef["shape"]), type = "l", col = "red")
legend("topleft", lty = 1, col = "red", legend = "Fitted Distribution", cex = 0.75, bty = "n")

par(mfrow = c(1, 1))

SP500.u <- fGarch::pstd(SP500.res, 0, 1, SP500.coef["shape"])
NK225.u <- fGarch::pstd(NK225.res, 0, 1, NK225.coef["shape"])

par(mfrow = c(1, 2))
hist(SP500.u, breaks = 20)
hist(NK225.u, breaks = 20)
par(mfrow = c(1, 1))

KScorrect::LcKS(SP500.u, cdf = "punif")$p.value
KScorrect::LcKS(NK225.u, cdf = "punif")$p.value
ADGofTest::ad.test(SP500.u, null = "punif")$p.value
ADGofTest::ad.test(NK225.u, null = "punif")$p.value

plot(SP500.u, NK225.u, pch = 20)

model <- VineCopula::BiCopSelect(SP500.u, NK225.u, familyset = NA, selectioncrit = "AIC", indeptest = TRUE, level = 0.05)
model

n <- 2000
sim.u <- VineCopula::BiCopSim(n, family = 2, model$par, model$par2)

plot(sim.u, pch = 20)

# Inverse PIT ==================================================================
SP500.sim <- qstd(sim.u[, 1], 0, 1, SP500.coef["shape"])
NK225.sim <- qstd(sim.u[, 2], 0, 1, NK225.coef["shape"])

# Re-introduce GARCH for SP500 =================================================
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

SP500.fin <- y
rm(y, mu, ar1, omega, alpha1, beta1, sigma2s)

# Re-introduce GARCH for NK225
mu <- NK225.coef["mu"]
# ar1 <- NK225.coef["ar1"]
omega <- NK225.coef["omega"]
alpha1 <- NK225.coef["alpha1"]
beta1 <- NK225.coef["beta1"]
sigma2s <- numeric(n)
y <- numeric(n)
sigma2s[1] <- 0
y[1] <- mu + NK225.sim[1]

for (i in 2:n) {
    sigma2s[i] <- omega + alpha1 * sigma2s[i - 1] * NK225.sim[i - 1]^2 + beta1 * sigma2s[i - 1]
    y[i] <- mu + sqrt(sigma2s[i]) * NK225.sim[i]
}

NK225.fin <- y

rm(y, mu, omega, alpha1, beta1, sigma2s)

plot(SP500, NK225, pch = 20, col = rgb(0, 0, 0, 0.1), main = "Distribution of Log-returns")
points(SP500.fin, NK225.fin, pch = 20, col = rgb(1, 0, 0, 0.1))
legend("bottomright", pch = 20, col = rgb(c(0, 1), 0, 0, 0.5), legend = c("Original Data", "Simulated Data"), cex = 0.75, bty = "n")

portsim <- matrix(0, nrow = n, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

portsim <- log(1 + ((exp(SP500.res) - 1) + (exp(NK225.res) - 1)) * (1 / 2))
varsim <- quantile(portsim, c(0.01, 0.05))
varsim
