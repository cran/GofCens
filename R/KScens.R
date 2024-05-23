KScens <- function(times, cens = rep(1, length(times)),
                   distr = c("exponential", "gumbel", "weibull", "normal",
                             "lognormal", "logistic", "loglogistic", "beta"),
                   betaLimits = c(0, 1), igumb = c(10, 10), degs = 3,
                   params0 = list(shape = NULL, shape2 = NULL,
                                  location = NULL, scale = NULL),
                   prnt = TRUE, outp ="list") {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  if (prnt && !outp %in% c("list", "table")) {
    stop("Invalid value of outp. Use 'table' or 'list'.")
  }
  if (!is.list(params0)) {
    stop("params0 must be a list!")
  }
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ",
                  ceiling(max(times) + 1), ")'.")
    stop(msg)
  }
  if (!all(sapply(params0, is.null))) {
    if (distr == "exponential" && is.null(params0$scale)) {
      stop("Argument 'params0' requires a value for the scale parameter.")
    }
    if (distr %in% c("weibull", "loglogistic") &&
        (is.null(params0$shape) || is.null(params0$scale))) {
      stop("Argument 'params0' requires values for the shape and scale parameters.")
    }
    if (distr %in% c("gumbel", "normal", "lognormal", "logistic") &&
        (is.null(params0$location) || is.null(params0$scale))) {
      stop("Argument 'params0' requires values for the location and scale parameters.")
    }
    if (distr == "beta" && (is.null(params0$shape) || is.null(params0$shape2))) {
      stop("Argument 'params0' requires values for both shape parameters.")
    }
  }
  n <- length(times)
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  alpha0 <- params0$shape
  gamma0 <- params0$shape2
  mu0 <- params0$location
  beta0 <- params0$scale
  alphaML <- gammaML <- muML <- betaML <- NULL
  if (distr == "exponential") {
    if (!is.null(beta0)) {
      hypo <- c(scale = beta0)
    }
    muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                       dist = "exponential")))
    betaML <- 1 / exp(-muu)
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pexp(x, 1/ beta)
    }
  }
  if (distr == "gumbel") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- try(suppressMessages(fitdistcens(dd, "gumbel",
                                                 start = list(alpha = igumb[1],
                                                              scale = igumb[2]))),
                    silent = TRUE)
    if (attr(paramsML, "class") == "try-error") {
      stop("Function failed to estimate the parameters.\n
            Try with other initial values.")
    }
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pgumbel(x, mu, beta)
    }
  }
  if (distr == "weibull") {
    if (!is.null(alpha0) && !is.null(beta0)) {
      hypo <- c(shape = alpha0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "weibull")
    alphaML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pweibull(x, alpha, beta)
    }
  }
  if (distr == "normal") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "norm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pnorm(x, mu, beta)
    }
  }
  if (distr == "lognormal") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "lnorm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - plnorm(x, mu, beta)
    }
  }
  if (distr == "logistic") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "logis")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - plogis(x, mu, beta)
    }
  }
  if (distr == "loglogistic") {
    if (!is.null(alpha0) && !is.null(beta0)) {
      hypo <- c(shape = alpha0, scale = beta0)
    }
    paramsML <- unname(survreg(Surv(times, cens) ~ 1,
                               dist = "loglogistic")$icoef)
    alphaML <- 1 / exp(paramsML[2])
    betaML <- exp(paramsML[1])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pllogis(x, alpha, scale = beta)
    }
  }
  if (distr == "beta") {
    if (!is.null(alpha0) && !is.null(gamma0)) {
      hypo <- c(shape = alpha0, shape2 = gamma0)
    }
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    paramsML <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
    alphaML <- unname(paramsML$estimate[1])
    gammaML <- unname(paramsML$estimate[2])
    SofT0 <- function(x, alpha, gamma, mu, beta) {
      1 - pbeta((x - aBeta) / (bBeta - aBeta), alpha, gamma)
    }
  }
  if (!all(sapply(params0, is.null))) {
    alpha <- alpha0
    gamma <- gamma0
    mu <- mu0
    beta <- beta0
  } else {
    alpha <- alphaML
    gamma <- gammaML
    mu <- muML
    beta <- betaML
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
  Ym <- sqrt(n) / 2 * (survT$surv[m] + SofT0(stimes[m], alpha, gamma, mu,
                        beta)) * (Avec[m] - Bvec[m])
  A <- max(abs(c(Yl, Y, Ym)))
  R <- 1 - 0.5 * (survT$surv[m] + SofT0(stimes[m], alpha, gamma, mu, beta))
  kr <- A / sqrt(R - R^2)
  sr <- sqrt((1 - R) / R)
  id <- 1:1000
  pvalue <- 2 * pnorm(-kr) - 2 * sum((-1)^id * exp(-2 * id^2 * A^2) *
                                       (pnorm(2 * id * A * sr + kr) -
                                          pnorm(2 * id * A * sr - kr)))
  if (all(sapply(params0, is.null))) {
    output <- list(Distribution = distr,
                   Test = round(c(A = A, "p-value" = pvalue, "F(ym)" = R,
                                  ym = stimes[m]), degs),
                   Estimates = round(c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML), degs))
  } else {
    output <- list(Distribution = distr,
                   Hypothesis = hypo,
                   Test = round(c(A = A, "p-value" = pvalue, "F(ym)" = R,
                                  ym = stimes[m]), degs),
                   Estimates = round(c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML), degs))
  }
  if (prnt){
    if (outp == "table") {
      cat("Distribution:", output$Distribution, "\n")
      if (!all(sapply(params0, is.null))) {
        cat("\nNull hypothesis:\n")
        header1 <- c("Parameter", "Value")
        max_col_width1 <- max(nchar(header1), nchar(names(output$Hypothesis)))
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
        cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                    max_col_width1, header1[2]))
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
        for (i in 1:length(output$Hypothesis)) {
          cat(sprintf("%-*s | %-*s\n", max_col_width1, names(output$Hypothesis)[i],
                      max_col_width1, unname(output$Hypothesis)[i]))
        }
        cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                    strrep("-", max_col_width1)))
      }
      cat("\nKS Test results:\n")
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
      cat("\nParameter estimates:\n")
      header1 <- c("Parameter", "Value")
      max_col_width1 <- max(nchar(header1), nchar(names(output$Estimates)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                  max_col_width1, header1[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      for (i in 1:length(output$Estimates)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width1, names(output$Estimates)[i],
                    max_col_width1, unname(output$Estimates)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      invisible(output)
    } else {
      cat("Distribution:", output$Distribution, "\n")
      if (!all(sapply(params0, is.null))) {
        cat("\nNull hypothesis:\n")
        print(output$Hypothesis)
      }
      cat("\nKS Test results:\n")
      print(output$Test)
      cat("\nParameter estimates:\n")
      print(output$Estimates)
      invisible(output)
    }
  } else{
    invisible(output)
  }
}
