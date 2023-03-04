# Functions ====================================================================
plot.ac <- function(x) {
    par(mfrow = c(2, 2))
    acf(x)
    pacf(x)
    acf(x^2)
    par(mfrow = c(1, 1))
}

# Load Data ====================================================================

CAC40 <- read.csv("data/CAC40.csv")
NK225 <- read.csv("data/NK225.csv")

CAC40.price <- CAC40$Adj.Close
NK225.price <- NK225$Price

CAC40.lr <- na.omit(log(tail(CAC40.price, -1)) - log(head(CAC40.price, -1)))
NK225.lr <- na.omit(log(tail(NK225.price, -1)) - log(head(NK225.price, -1)))

rm(CAC40, NK225, CAC40.price, NK225.price)

hist(CAC40.lr, breaks = 20, freq = FALSE)
hist(NK225.lr, breaks = 20, freq = FALSE)
plot(CAC40.lr, NK225.lr)

# TS Adjustment ================================================================

Box.test(CAC40.lr, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Significant
Box.test(NK225.lr, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(CAC40.lr)

CAC40.model <- fGarch::garchFit(formula = ~arma(1, 0) + garch(1, 1), data = CAC40.lr, trace = FALSE, cond.dist = "norm")
CAC40.res <- residuals(CAC40.model, standardize = TRUE)
Box.test(CAC40.res, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
Box.test(CAC40.res^2, lag = 10, type = c("Ljung-Box"), fitdf = 1) # Passed
plot.ac(CAC40.res)

# PIT ==========================================================================

CAC40.params <- MASS::fitdistr(CAC40.res, "t")$estimate
NK225.params <- MASS::fitdistr(NK225.lr, "t")$estimate

hist(CAC40.res, breaks = 20, freq = FALSE)
x <- seq(-4, 4, 0.001)
lines(x, dstd(x, CAC40.params[1], CAC40.params[2]*1.1, CAC40.params[3]), type = "l")

hist(NK225.lr, breaks = 20, freq = FALSE)
x <- seq(-0.4, 0.4, 0.001)
lines(x, dstd(x, NK225.params[1], NK225.params[2]*1.2, NK225.params[3]), type = "l")

CAC40.u <- pstd(CAC40.res, CAC40.params[1], CAC40.params[2]*1.1, CAC40.params[3])
NK225.u <- pstd(NK225.lr, NK225.params[1], NK225.params[2]*1.2, NK225.params[3])

plot(CAC40.u, NK225.u)
