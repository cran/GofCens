gofcens <-
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
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  stimes <- sort(unique(times[cens == 1]))
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
    y0 <- c(pexp(stimes, 1 / beta), 1)
  }
  if (distr == "gumbel") {
    if (is.null(mu) || is.null(beta)) {
      param <- try(suppressMessages(fitdistcens(dd, "gumbel",
                                                start = list(alpha = igumb[1],
                                                             scale = igumb[2]))),
                   silent = TRUE)
      if (attr(param, "class") == "try-error") {
        stop("Function failed to estimate the parameters.\n
              Try with other initial values.")
      }
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    y0 <- c(pgumbel(stimes, mu, beta), 1)
  }
  if (distr == "weibull") {
    if (is.null(alpha) || is.null(beta)) {
      param <- fitdistcens(dd, "weibull")
      alpha <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    y0 <- c(pweibull(stimes, alpha, beta), 1)
  }
  if (distr == "normal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "norm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    y0 <- c(pnorm(stimes, mu, beta), 1)
  }
  if (distr == "lognormal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "lnorm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    y0 <- c(plnorm(stimes, mu, beta), 1)
  }
  if (distr == "logistic") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "logis")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    y0 <- c(plogis(stimes, mu, beta), 1)
  }
  if (distr == "loglogistic") {
    if (is.null(alpha) || is.null(beta)) {
      param <- unname(survreg(Surv(times, cens) ~ 1,
                              dist = "loglogistic")$icoef)
      alpha <- 1 / exp(param[2])
      beta <- exp(param[1])
    }
    y0 <- c(pllogis(stimes, alpha, scale = beta), 1)
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    if (is.null(alpha) || is.null(gamma)) {
      param <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      alpha <- unname(param$estimate[1])
      gamma <- unname(param$estimate[2])
    }
    y0 <- c(pbeta((stimes - aBeta) / (bBeta - aBeta), alpha, gamma), 1)
  }
  KM <- summary(survfit(Surv(times, cens) ~ 1))$surv
  nc <- length(KM)
  Fn <- c(1 - KM, NA)
  CvM <- nc * (sum(Fn[-(nc + 1)] * (y0[-1] - y0[-(nc + 1)]) *
                   (Fn[-(nc + 1)] - (y0[-1] + y0[-(nc + 1)]))) + 1 / 3)
  Fn <- Fn[-(nc + 1)]
  y0 <- y0[-(nc + 1)]
  AD <- nc * (-1 - log(y0[nc]) - log(1 - y0[nc])+
               sum(Fn[-nc]^2 *
                   (-log(1 - y0[-1]) + log(y0[-1]) + log(1 - y0[-nc]) - log(y0[-nc]))) -
               2 * sum(Fn[-nc] * (-log(1 - y0[-1]) + log(1 - y0[-nc]))))
  KS <- as.vector(KScens(times, cens, distr, betaLimits)$Test[2])
  output <- list("Test statistics" = round(c(KS = KS, CvM = CvM, AD = AD), degs),
                 Distribution = distr,
                 Parameters = round(c(shape = alpha, shape2 = gamma,
                                      location = mu, scale = beta), degs))
  return(output)
}
