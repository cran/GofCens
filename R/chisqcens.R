chisqcens <-
  function(times, cens = rep(1, length(times)), M,
           distr = c("exponential", "gumbel", "weibull", "normal",
                     "lognormal", "logistic", "loglogistic", "beta"),
           betaLimits=c(0, 1), igumb = c(10, 10), degs = 3, BS = 999,
           params = list(shape = NULL, shape2 = NULL,
                         location = NULL, scale = NULL),
           prnt = TRUE, outp = "list", tol = 1e-04) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Status indicator must be either 0 or 1!")
  }
  if (prnt && !outp %in% c("list", "table")) {
    stop("Invalid value of outp. Use 'table' or 'list'.")
  }
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ", ceiling(max(times) + 1), ")'.")
    stop(msg)
  }
  n <- length(times)
  rnd <- -log(tol, 10)
  times <- round(pmax(times, tol), rnd)
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  survKM <- survfit(Surv(times, cens) ~ 1)
  alpha <- params$shape
  gamma <- params$shape2
  mu <- params$location
  beta <- params$scale
  censKM <- survfit(Surv(times, 1 - cens) ~ 1)
  Morig <- M
  cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
  if (anyNA(cb)) {
    cb <- cb[!is.na(cb)]
  }
  Mout <- length(cb) - 1
  if (distr == "exponential") {
    if (is.null(beta)) {
      muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                         dist = "exponential")))
      beta <- 1 / exp(-muu)
    }
    expStat <- function(dat) {
      muu <- unname(coefficients(survreg(Surv(dat$times, dat$cens) ~ 1,
                                         dist = "exponential")))
      betahat <- 1 / exp(-muu)
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      expProb <- diff(pexp(cb[1:(Mred + 1)], 1 / betahat))
      v <- (obsfreq - n * expProb) / sqrt(n * expProb)
      tn <- t(v) %*% v
    }
    expRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rexp(n, mle),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), expStat, R = BS, sim = "parametric",
                ran.gen = expRnd, mle = 1 / beta)
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
    gumbStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens(dd, "gumbel", start = list(alpha = mu,
                                                       scale = beta))
      muhat <- unname(params$estimate[1])
      betahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      gumbProb <- diff(pgumbel(cb[1:(M + 1)], muhat, betahat))
      gumbProb[1] <- gumbProb[1] + pgumbel(cb[1], muhat, betahat)
      v <- (obsfreq - n * gumbProb) / sqrt(n * gumbProb)
      tn <- t(v) %*% v
    }
    gumbRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rgumbel(n, mle[1], mle[2]), tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), gumbStat, R = BS, sim = "parametric",
                ran.gen = gumbRnd, mle = c(mu, beta))
  }
  if (distr == "weibull") {
    if (is.null(alpha) || is.null(beta)) {
      param <- fitdistcens(dd, "weibull")
      alpha <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    weiStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens(dd, "weibull")
      alphahat <- unname(params$estimate[1])
      betahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      weiProb <- diff(pweibull(cb[1:(Mred + 1)], alphahat, betahat))
      v <- (obsfreq - n * weiProb) / sqrt(n * weiProb)
      tn <- t(v) %*% v
    }
    weiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rweibull(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), weiStat, R = BS, sim = "parametric",
                ran.gen = weiRnd, mle = c(alpha, beta))
  }
  if (distr == "normal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "norm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    normStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens(dd, "norm")
      muhat <- unname(params$estimate[1])
      betahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      normProb <- diff(pnorm(cb[1:(Mred + 1)], muhat, betahat))
      normProb[1] <- normProb[1] + pnorm(cb[1], muhat, betahat)
      v <- (obsfreq - n * normProb) / sqrt(n * normProb)
      tn <- t(v) %*% v
    }
    normRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rnorm(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), normStat, R = BS, sim = "parametric",
                ran.gen = normRnd, mle = c(mu, beta))
  }
  if (distr == "lognormal") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "lnorm")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    lnormStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens(dd, "lnorm")
      muhat <- unname(params$estimate[1])
      betahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      lnormProb <- diff(plnorm(cb[1:(Mred + 1)], muhat, betahat))
      v <- (obsfreq - n * lnormProb) / sqrt(n * lnormProb)
      tn <- t(v) %*% v
    }
    lnormRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rlnorm(n, mle[1], mle[2]), tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), lnormStat, R = BS, sim = "parametric",
                ran.gen = lnormRnd, mle = c(mu, beta))
  }
  if (distr == "logistic") {
    if (is.null(mu) || is.null(beta)) {
      param <- fitdistcens(dd, "logis")
      mu <- unname(param$estimate[1])
      beta <- unname(param$estimate[2])
    }
    logiStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens(dd, "logis")
      muhat <- unname(params$estimate[1])
      betahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      logiProb <- diff(plogis(cb[1:(Mred + 1)], muhat, betahat))
      logiProb[1] <- logiProb[1] + plogis(cb[1], muhat, betahat)
      v <- (obsfreq - n * logiProb) / sqrt(n * logiProb)
      tn <- t(v) %*% v
      }
    logiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rlogis(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), logiStat, R = BS, sim = "parametric",
                ran.gen = logiRnd, mle = c(mu, beta))
  }
  if (distr == "loglogistic") {
    if (is.null(alpha) || is.null(beta)) {
      param <- unname(survreg(Surv(times, cens) ~ 1,
                              dist = "loglogistic")$icoef)
      alpha <- 1 / exp(param[2])
      beta <- exp(param[1])
    }
    llogiStat <- function(dat) {
      params <- unname(survreg(Surv(dat$times, dat$cens) ~ 1,
                               dist = "loglogistic")$icoef)
      alphahat <- 1 / exp(params[2])
      betahat <- exp(params[1])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      llogiProb <- diff(pllogis(cb[1:(Mred + 1)], alphahat, scale = betahat))
      v <- (obsfreq - n * llogiProb) / sqrt(n * llogiProb)
      tn <- t(v) %*% v
    }
    llogiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rllogis(n, mle[1], scale = mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), llogiStat, R = BS, sim = "parametric",
                ran.gen = llogiRnd, mle = c(alpha, beta))
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    if (is.null(alpha) || is.null(gamma)) {
      param <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      alpha <- unname(param$estimate[1])
      gamma <- unname(param$estimate[2])
    }
    betaStat <- function(dat) {
      dd <- data.frame(left = as.vector(dat$times),
                       right = ifelse(dat$cens == 1, dat$times, NA))
      params <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      alphahat <- unname(params$estimate[1])
      gammahat <- unname(params$estimate[2])
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      betaProb <- diff(pbeta(cb[1:(M + 1)] * (bBeta - aBeta) + aBeta, alphahat,
                             gammahat))
      v <- (obsfreq - n * betaProb) / sqrt(n * betaProb)
      tn <- t(v) %*% v
    }
    betaRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rbeta(n, alpha, gamma) * (bBeta - aBeta) + aBeta,
                              tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    bts <- boot(data.frame(times, cens), betaStat, R = BS, sim = "parametric",
                ran.gen = betaRnd, mle = c(alpha, gamma))
  }
  tn <- bts$t0
  pval <- (sum(bts$t[, 1] > bts$t0[1]) + 1) / (bts$R + 1)
  output <- list(Test = round(c("Statistic" = tn, "p-value" = pval), degs),
                 Distribution = distr,
                 Parameters = round(c(shape = alpha, shape2 = gamma,
                                      location = mu, scale = beta), degs),
                 Cellnumber = c("Original" = Morig, "Final" = Mout))
  if (prnt) {
    if (outp == "table") {
      cat("Distribution: ", output$Distribution, "\n")
      cat("\nCvM Test Results:\n")
      header <- c("Metric", "Value")
      max_col_width <- max(nchar(header), nchar(names(output$Test)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      cat(sprintf("%-*s | %-*s\n", max_col_width, header[1],
                  max_col_width, header[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      for (i in 1:length(output$Test)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width, names(output$Test)[i],
                    max_col_width, unname(output$Test)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                  strrep("-", max_col_width)))
      cat("\nDistribution Parameters:\n")
      header1 <- c("Parameter", "Value")
      max_col_width1 <- max(nchar(header1), nchar(names(output$Parameters)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                  max_col_width1, header1[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      for (i in 1:length(output$Parameters)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width1, names(output$Parameters)[i],
                    max_col_width1, unname(output$Parameters)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      invisible(output)
    } else {
      cat("\nDistribution: ", output$Distribution, "\n")
      cat("\nChi-squared Test Results:\n")
      print(output$Test)
      cat("\nDistribution Parameters:\n")
      print(output$Parameters)
      cat("\nCell numbers:\n")
      print(output$Cellnumber)
      invisible(output)
    }
  } else {
    invisible(output)
  }
}
