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
KScens(times, delta, distr = "lognormal")

## -----------------------------------------------------------------------------
set.seed(123)
summary(KScens(times, delta, distr = "weibull"), outp = "table")

## -----------------------------------------------------------------------------
KScens(times, delta, distr = "lognormal", boot = FALSE)

## -----------------------------------------------------------------------------
set.seed(123)
CvMcens(times, delta, distr = "lognormal")

## -----------------------------------------------------------------------------
set.seed(123)
CvMcens(times, delta, distr = "weibull")

## -----------------------------------------------------------------------------
set.seed(123)
summary(ADcens(times, delta, distr = "lognormal", 
             params0 = list(location = 2, scale = 1.5)), outp = "table")

## -----------------------------------------------------------------------------
set.seed(123)
gofcens(times, delta, distr = "lognormal")

## -----------------------------------------------------------------------------
set.seed(123)
summary(gofcens(times, delta, distr = "weibull"), outp = "table")

## -----------------------------------------------------------------------------
set.seed(123)
chisqcens(times, delta, M = 8, distr = "lognormal")

## -----------------------------------------------------------------------------
set.seed(123)
summary(chisqcens(times, delta, M = 8, distr = "weibull"), outp = "table")

## -----------------------------------------------------------------------------
data("nba")
set.seed(123)
nbasamp <- nba[sample(nrow(nba), 500), ]
set.seed(123)
gofcens(Surv(survtime, cens) ~ 1, nbasamp, distr = "logistic")

## -----------------------------------------------------------------------------
set.seed(123)
gofcens(Surv(survtime, cens) ~ 1, nbasamp, distr = "normal")

