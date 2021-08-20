chisqcens2 <-
function(times, cens, M, distrData = c("weibull", "lognormal", "loglogistic"),
         distrCens = c("weibull", "lognormal", "loglogistic", "uniform"), 
         BS = 1000, degs = 4, 
         params = list(shape = NULL, location = NULL, scale = NULL)) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  distrData <- match.arg(distrData)
  distrCens <- match.arg(distrCens)
  n <- length(times)
  m <- max(times)
  est <- chisqcens1(times, cens, M, distrData, params)
  tn <- est$Statistic
  MF <- unname(est$Cellnumber[2])
  alpha <- mu <- beta <- NULL
  if ("shape" %in% names(est$Parameters)) {
    alpha <- unname(est$Parameters["shape"])
  }
  if ("location" %in% names(est$Parameters)) {
    mu <- unname(est$Parameters["location"])
  }
  if ("scale" %in% names(est$Parameters)) {
    beta <- unname(est$Parameters["scale"])
  }
  censC <- 1 - cens
  estCens <- chisqcens1(times, censC, M, distrCens)
  alphaCens <- gammaCens <- muCens <- betaCens <- NULL
  if ("shape" %in% names(estCens$Parameters)) {
    alphaCens <- unname(estCens$Parameters["shape"])
  }
  if ("shape2" %in% names(estCens$Parameters)) {
    gammaCens <- unname(estCens$Parameters["shape2"])
  }
  if ("location" %in% names(estCens$Parameters)) {
    muCens <- unname(estCens$Parameters["location"])
  }
  if ("scale" %in% names(estCens$Parameters)) {
    betaCens <- unname(estCens$Parameters["scale"])
  }
  t <- numeric(BS)
  if (distrData == "weibull") {
    if (distrCens == "weibull") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "weibull", alpha, -log(beta),
                                dist.cens = "weibull", alphaCens, -log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "weibull",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
    if (distrCens == "lognormal") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "weibull", alpha, -log(beta),
                                dist.cens = "lnorm", betaCens, muCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "weibull",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
    if (distrCens == "loglogistic") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "weibull", alpha, -log(beta),
                                dist.cens = "llogistic", 1 / alphaCens,
                                log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "weibull",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
    if (distrCens == "unif") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "weibull", alpha, -log(beta),
                                dist.cens = "unif", gammaCens, alphaCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "weibull",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
  }
  if (distrData == "lognormal") {
    if (distrCens == "weibull") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf, dist.ev = "lnorm", beta, mu,
                                dist.cens = "weibull", alphaCens, -log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "lognormal",
                           params = list(shape = alpha, location = mu))$Statistic
      }
    }
    if (distrCens == "lognormal") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf, dist.ev = "lnorm", beta, mu,
                                dist.cens = "lnorm", betaCens, muCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "lognormal",
                           params = list(shape = alpha, location = mu))$Statistic
      }
    }
    if (distrCens == "loglogistic") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf, dist.ev = "lnorm", beta, mu,
                                dist.cens = "llogistic", 1 / alphaCens,
                                log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "lognormal",
                           params = list(shape = alpha, location = mu))$Statistic
      }
    }
    if (distrCens == "unif") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf, dist.ev = "lnorm", beta, mu,
                                dist.cens = "unif", gammaCens, alphaCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "lognormal",
                           params = list(shape = alpha, location = mu))$Statistic
      }
    }
  }
  if (distrData == "loglogistic") {
    if (distrCens == "weibull") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "llogistic", 1 / alpha,
                                log(beta), dist.cens = "weibull", alphaCens,
                                -log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "loglogistic",
                           params = list(shape = alpha, location = mu))$Statistic
      }
    }
    if (distrCens == "lognormal") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "llogistic", 1 / alpha,
                                log(beta), dist.cens = "lnorm", betaCens, muCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "loglogistic",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
    if (distrCens == "loglogistic") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "llogistic", 1 / alpha,
                                log(beta), dist.cens = "llogistic",
                                1 / alphaCens, log(betaCens))
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "loglogistic",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
    if (distrCens == "unif") {
      for (i in 1:BS) {
        rand <- simple.surv.sim(n, Inf,  dist.ev = "llogistic", 1 / alpha,
                                log(beta), dist.cens = "unif", gammaCens,
                                alphaCens)
        t[i] <- chisqcens1(rand$stop, rand$status, M, distr = "loglogistic",
                           params = list(shape = alpha, scale = beta))$Statistic
      }
    }
  }
  pvalue <- 1 - ecdf(t)(tn)
  output <- list(test = c("Statistic" = tn), "p-value" = pvalue,
                 distrData = distrData, distrCens = distrCens,
                 Parameters = round(c(shape = alpha, location = mu,
                                      scale = beta), degs),
                 Cellnumber = c("Original" = M, "Final" = MF))
  return(output)
}
