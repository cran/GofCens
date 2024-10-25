## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("GofCens")

## -----------------------------------------------------------------------------
set.seed(123)
survt <- round(rlnorm(300, 2, 1), 2)
censt <- round(rexp(300, 1 / 20), 2)

## -----------------------------------------------------------------------------
times <- pmin(survt, censt)
delta <- as.numeric(survt <= censt)

## -----------------------------------------------------------------------------
table(delta)

## -----------------------------------------------------------------------------
probPlot(times, delta, distr = "lognormal", prnt = FALSE, cex.lab = 1.3)
probPlot(times, delta, distr = "weibull", ggp = TRUE, prnt = FALSE)

## -----------------------------------------------------------------------------
probPlot(Surv(times, delta) ~ 1, distr = "lognormal", m = matrix(1:4, nrow = 1),
         params0 = list(location = 2, scale = 1.5), ggp = TRUE)

## -----------------------------------------------------------------------------
 cumhazPlot(times, delta, font.lab = 4, cex.lab = 1.3)

## -----------------------------------------------------------------------------
cumhazPlot(times, delta, distr = c("exponential", "beta", "lognormal"), 
           betaLimits = c(0, 100), ggp = TRUE, prnt = FALSE)


## -----------------------------------------------------------------------------
kmPlot(times, delta, ggp = TRUE, prnt = FALSE)

## -----------------------------------------------------------------------------
data("nba")
cumhazPlot(Surv(survtime, cens) ~ 1, nba, font.lab = 4, cex.lab = 1.3, 
           prnt = FALSE, lwd = 3, colour = "blue")

## -----------------------------------------------------------------------------
probPlot(Surv(survtime, cens) ~ 1, nba, distr = "logistic", ggp = TRUE, 
         degs = 2)
probPlot(Surv(survtime, cens) ~ 1, nba, distr = "normal", ggp = TRUE,
         prnt = FALSE)

