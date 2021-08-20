chisqcens1 <-
function(times, cens = rep(1, length(times)), M,
                       distr = c("exponential", "gumbel", "weibull", "normal",
                                 "lognormal", "logistic", "loglogistic", "beta",
                                 "uniform"),
                       betaLimits=c(0, 1), igumb = c(10, 10), degs = 4,
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
  survKM <- survfit(Surv(times, cens) ~ 1)
  alpha <- params$shape
  gamma <- params$shape2
  mu <- params$location
  beta <- params$scale
  Morig <- M
  cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
  if (anyNA(cb)) {
    cb <- cb[!is.na(cb)]
  }
  M <- length(cb) - 1
  cb[M + 1] <- cb[M + 1] + 1
  cellsCut <- cut(times[cens == 1], cb, right = FALSE)
  obsfreq <- as.vector(table(cellsCut))
  if (distr == "exponential") {
    if (is.null(beta)) {
      muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                         dist = "exponential")))
      beta <- exp(-muu)
    }
    expProb <- pexp(cb[1:M + 1], beta) - pexp(cb[1:M], beta)
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
    expProb <- pgumbel(cb[1:M + 1], beta, mu) - pgumbel(cb[1:M], beta, mu)
  }
  if (distr == "weibull") {
    if (is.null(alpha) || is.null(beta)) {
      param <- fitdistcens(dd, "weibull")
      alpha <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    expProb <- pweibull(cb[1:M + 1], alpha, beta) - pweibull(cb[1:M], alpha, beta)
  }
  if (distr == "normal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "norm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    expProb <- pnorm(cb[1:M + 1], mu, beta) - pnorm(cb[1:M], mu, beta)
  }
  if (distr == "lognormal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "lnorm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    expProb <- plnorm(cb[1:M + 1], mu, beta) - plnorm(cb[1:M], mu, beta)
  }
  if (distr == "logistic") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "logis")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    expProb <- plogis(cb[1:M + 1], mu, beta) - plogis(cb[1:M], mu, beta)
  }
  if (distr == "loglogistic") {
    if (is.null(alpha) || is.null(beta)) {
      param <- unname(survreg(Surv(times, cens) ~ 1,
                              dist = "loglogistic")$icoef)
      alpha <- 1 / exp(param[2])
      beta <- exp(param[1])
    }
    expProb <- pllogis(cb[1:M + 1], alpha, beta) - pllogis(cb[1:M], alpha, beta)
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    if (is.null(alpha) || is.null(gamma)) {
      param <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      alpha <- unname(param$estimate[1])
      gamma <- unname(param$estimate[2])
    }
    expProb <- pbeta((cb[1:M + 1] - aBeta) / (bBeta - aBeta), alpha, gamma) -
               pbeta((cb[1:M] - aBeta) / (bBeta - aBeta), alpha, gamma)
  }
  if (distr == "uniform") {
    if (is.null(alpha) || is.null(gamma)) {
      param <- fitdistcens(dd, "unif")
      alpha <- unname(param$estimate[1])
      gamma <- unname(param$estimate[2])
    }
    expProb <- punif(cb[1:M + 1], alpha, gamma) - punif(cb[1:M], alpha, gamma)
  }
  if (is.element(0, expProb)) {
    stop("Some of the expected probabilities are 0.")
  }
  v <- (obsfreq - n * expProb) / sqrt(n * expProb)
  tn <- as.vector(t(v) %*% v)
  output <- list(Statistic = tn, Distribution = distr,
                 Parameters = round(c(shape = alpha, shape2 = gamma,
                                      location = mu, scale = beta), degs),
                 Cellnumber = c("Original" = Morig, "Final" = M))
  return(output)
}
