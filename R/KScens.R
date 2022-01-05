KScens <-
function(times, cens = rep(1, length(times)),
                   distr = c("exponential", "gumbel", "weibull", "normal",
                             "lognormal", "logistic", "loglogistic", "beta"),
                   betaLimits = c(0, 1), igumb = c(10, 10), degs = 4,
                   params = list(shape = NULL, shape2 = NULL,
                                 location = NULL, scale = NULL)) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ", ceiling(max(times) + 1), ")'.")
    stop(msg)
  }
  n <- length(times)
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  alpha <- params$shape
  gamma <- params$shape2
  mu <- params$location
  beta <- params$scale
  if (distr == "exponential") {
    if (is.null(beta)) {
      muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                         dist = "exponential")))
      beta <- 1 / exp(-muu)
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pexp(x, 1/ beta)
    }
  }
  if (distr == "gumbel") {
    if (is.null(mu) || is.null(beta)) {
      param <- try(suppressMessages(fitdistcens(dd, "gumbel",
                                                 start = list(alpha = igumb[1],
                                                              scale = igumb[2]))),
                   silent = TRUE)
      if (attr(param, "class") == "try-error") {
        stop("Function failed to estimate the parameters. Try with other initial values.")
      }
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pgumbel(x, mu, beta)
    }
  }
  if (distr == "weibull") {
    if (is.null(alpha) || is.null(beta)) {
      param <- fitdistcens(dd, "weibull")
      alpha <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pweibull(x, alpha, beta)
    }
  }
  if (distr == "normal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "norm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pnorm(x, mu, beta)
    }
  }
  if (distr == "lognormal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "lnorm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - plnorm(x, mu, beta)
    }
  }
  if (distr == "logistic") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "logis")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - plogis(x, mu, beta)
    }
  }
  if (distr == "loglogistic") {
    if (is.null(alpha) || is.null(beta)) {
      param <- unname(survreg(Surv(times, cens) ~ 1,
                              dist = "loglogistic")$icoef)
      alpha <- 1 / exp(param[2])
      beta <- exp(param[1])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pllogis(x, alpha, scale = beta)
    }
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    if (is.null(alpha) || is.null(gamma)) {
      param <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      alpha <- unname(param$estimate[1])
      gamma <- unname(param$estimate[2])
    }
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pbeta((x - aBeta) / (bBeta - aBeta), alpha, gamma)
    }
  }
  sumSurvT <- survfit(Surv(times, cens) ~ 1, stype = 2, ctype = 2)
  survT <- unique(data.frame(times = sumSurvT$time, surv = sumSurvT$surv))
  stimes <- survT$time
  m <- length(stimes)
  svbefor <- c(1, survT$surv[-m])
  aux2 <- numeric(m)
  for (i in 1:m) {
    if (sumSurvT$n.censor[i] > 0) {
      aux2[i] <- with(sumSurvT, sum(1 / (n.risk[i] - n.event[i] -
                                         (0:(n.censor[i] - 1)))))
    }
  }
  alfatj <- exp(-c(0, cumsum(aux2))[-m])
  Atj <- sqrt(c(1, alfatj[-m])) *
         log(SofT0(c(0, stimes[-m]), alpha, gamma, mu, beta) /
             SofT0(stimes, alpha, gamma, mu, beta))
  Atj[is.nan(Atj)] <- 0
  Avec <- cumsum(Atj)
  Btj <- sqrt(c(1, alfatj[-m])) * log(svbefor / survT$surv)
  Btj[is.nan(Btj)] <- 0
  Bvec <- cumsum(Btj)
  Yl <- sqrt(n) / 2 * (svbefor + SofT0(stimes, alpha, gamma, mu, beta)) *
                      (Avec - c(0, Bvec[-m])) * ifelse(Bvec > 0, 1, 0)
  Y <-  sqrt(n) / 2 * (survT$surv + SofT0(stimes, alpha, gamma, mu, beta)) *
                      (Avec - Bvec) * ifelse(Bvec > 0, 1, 0)
  Ym <- sqrt(n) / 2 * (survT$surv[m] + SofT0(stimes[m], alpha, gamma, mu, beta)) *
                      (Avec[m] - Bvec[m])
  A <- max(abs(c(Yl, Y, Ym)))
  R <- 1 - 0.5 * (survT$surv[m] + SofT0(stimes[m], alpha, gamma, mu, beta))
  kr <- A / sqrt(R - R^2)
  sr <- sqrt((1 - R) / R)
  id <- 1:1000
  pvalue <- 2 * pnorm(-kr) - 2 * sum((-1)^id * exp(-2 * id^2 * A^2) *
                                     (pnorm(2 * id * A * sr + kr) -
                                      pnorm(2 * id * A * sr - kr)))
  output <- list(Test = round(c("p-value" = pvalue, A = A, "F(ym)" = R,
                                ym = stimes[m]), degs),
                 Distribution = distr,
                 Parameters = round(c(shape = alpha, shape2 = gamma,
                                      location = mu, scale = beta), degs))
  return(output)
}
