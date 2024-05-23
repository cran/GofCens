probPlot <- function(times, cens = rep(1, length(times)),
                     distr = c("exponential", "gumbel", "weibull", "normal",
                               "lognormal", "logistic", "loglogistic", "beta"),
                     plots = c("PP", "QQ", "SP", "ER"),
                     colour = c("green4", "deepskyblue4", "yellow3",
                                "mediumvioletred"), mtitle = TRUE, ggp = FALSE,
                     m = NULL, betaLimits = c(0, 1), igumb = c(10, 10),
                     degs = 3, prnt = TRUE,
                     params0 = list(shape = NULL, shape2 = NULL,
                                    location = NULL, scale = NULL), ...) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  if (!is.logical(ggp) || !is.logical(prnt)) {
    stop("ggp and prnt must be logicals!")
  }
  if (!is.list(params0)) {
    stop("params0 must be a list!")
  }
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ", ceiling(max(times) + 1), ")'.")
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
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  survKM <- survfit(Surv(times, cens) ~ 1)
  tim <- summary(survKM)$time
  survTim <- summary(survKM)$surv
  uPointSurv <- survfit(Surv(tim, rep(1, length(tim))) ~ 1)
  uPoint <- 1 - uPointSurv$surv
  empiricF <- rbind(c(0, tim, Inf), c(0, uPoint, 1))
  uEstim <- rep(0, length(uPoint))
  alpha0 <- params0$shape
  gamma0 <- params0$shape2
  mu0 <- params0$location
  beta0 <- params0$scale
  alphaML <- gammaML <- muML <- betaML <- NULL
  if (distr == "exponential") {
    muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                       dist = "exponential")))
    betaML <- 1 / exp(-muu)
    if (is.null(beta0)) {
      rateExp <- exp(-muu)
      outp <- list(Distribution = "exponential", Estimates = betaML)
    } else {
      rateExp <- 1 / beta0
      hypo <- c(scale = beta0)
      outp <- list(Distribution = "exponential", Parameters = hypo,
                   Estimates = betaML)
    }
    theorPP <- pexp(tim, rateExp)
    theorQQ <- qexp(1 - survTim, rateExp)
  }
  if (distr == "gumbel") {
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
    if (is.null(mu0) || is.null(beta0)) {
      locGum <- muML
      scaleGum <- betaML
      outp <- list(Distribution = "Gumbel",
                   Estimates = c(location = muML, scale = betaML))
    } else {
      locGum <- mu0
      scaleGum <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "Gumbel", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML))
    }
    theorPP <- pgumbel(tim, locGum, scaleGum)
    theorQQ <- qgumbel(1 - survTim, locGum, scaleGum)
  }
  if (distr == "weibull") {
    paramsML <- fitdistcens(dd, "weibull")
    alphaML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    if (is.null(alpha0) || is.null(beta0)) {
      shapeWei <- alphaML
      scaleWei <- betaML
      outp <- list(Distribution = "Weibull",
                   Estimates = c(shape = alphaML, scale = betaML))
    } else {
      shapeWei <- alpha0
      scaleWei <- beta0
      hypo <- c(shape = alpha0, scale = beta0)
      outp <- list(Distribution = "Weibull", Parameters = hypo,
                   Estimates = c(shape = alphaML, scale = betaML))
    }
    theorPP <- pweibull(tim, shapeWei, scaleWei)
    theorQQ <- qweibull(1 - survTim, shapeWei, scaleWei)
  }
  if (distr == "normal") {
    paramsML <- fitdistcens(dd, "norm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    if (is.null(mu0) || is.null(beta0)) {
      locNorm <- muML
      scaleNorm <- betaML
      outp <- list(Distribution = "normal",
                   Estimates = c(location = muML, scale = betaML))
    } else {
      locNorm <- mu0
      scaleNorm <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "normal", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML))
    }
    theorPP <- pnorm(tim, locNorm, scaleNorm)
    theorQQ <- qnorm(1 - survTim, locNorm, scaleNorm)
  }
  if (distr == "lognormal") {
    paramsML <- fitdistcens(dd, "lnorm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    if (is.null(mu0) || is.null(beta0)) {
      locLnorm <- muML
      scaleLnorm <- betaML
      outp <- list(Distribution = "log-normal",
                   Estimates = c(location = muML, scale = betaML))
    } else {
      locLnorm <- mu0
      scaleLnorm <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "log-normal", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML))
    }
    theorPP <- plnorm(tim, locLnorm, scaleLnorm)
    theorQQ <- qlnorm(1 - survTim, locLnorm, scaleLnorm)
  }
  if (distr == "logistic") {
    paramsML <- fitdistcens(dd, "logis")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    if (is.null(mu0) || is.null(beta0)) {
      locLogis <- muML
      scaleLogis <- betaML
      outp <- list(Distribution = "logistic",
                   Estimates = c(location = muML, scale = betaML))
    } else {
      locLogis <- mu0
      scaleLogis <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "logistic", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML))
    }
    theorPP <- plogis(tim, locLogis, scaleLogis)
    theorQQ <- qlogis(1 - survTim, locLogis, scaleLogis)
  }
  if (distr == "loglogistic") {
    paramsML <- unname(survreg(Surv(times, cens) ~ 1,
                               dist = "loglogistic")$icoef)
    alphaML <- 1 / exp(paramsML[2])
    betaML <- exp(paramsML[1])
    if (is.null(alpha0) || is.null(beta0)) {
      shapeLoglog <- alphaML
      scaleLoglog <- betaML
      outp <- list(Distribution = "log-logistic",
                   Estimates = c(shape = alphaML, scale = betaML))
    } else {
      shapeLoglog <- alpha0
      scaleLoglog <- beta0
      hypo <- c(shape = alpha0, scale = beta0)
      outp <- list(Distribution = "log-logistic", Parameters = hypo,
                   Estimates = c(shape = alphaML, scale = betaML))
    }
    theorPP <- pllogis(tim, shapeLoglog, scale = scaleLoglog)
    theorQQ <- qllogis(1 - survTim, shapeLoglog, scale = scaleLoglog)
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    paramsML <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
    alphaML <- unname(paramsML$estimate[1])
    gammaML <- unname(paramsML$estimate[2])
    if (is.null(alpha0) || is.null(gamma0)) {
      shape1Beta <- alphaML
      shape2Beta <- gammaML
      outp <- list(Distribution = "beta",
                   Estimates = c(shape = alphaML, shape2 = gammaML),
                   interval.domain = betaLimits)
    } else {
      shape1Beta <- alpha0
      shape2Beta <- gamma0
      hypo <- c(shape = alpha0, shape2 = gamma0)
      outp <- list(Distribution = "beta", Parameters = hypo,
                   Estimates = c(shape = alphaML, shape2 = gammaML),
                   interval.domain = betaLimits)
    }
    theorPP <- pbeta((tim - aBeta)/(bBeta - aBeta), shape1Beta, shape2Beta)
    theorQQ <- qbeta((1 - survTim), shape1Beta, shape2Beta) * (bBeta - aBeta)
               + aBeta
}
  index <- 0
  for (i in 1:length(uPoint)) {
    while (theorQQ[i] > empiricF[1, index + 1]) {
      index <- index + 1
    }
    if (index != 0) {
      uEstim[i] <- empiricF[2, index]
    }
  }
  howmany <- length(plots)
  nCol <- length(colour)
  if (is.null(m)) {
    if (howmany == 1) {
      m <- matrix(1)
    } else {
      if (howmany == 2) {
        m <- matrix(1:2, nrow = 1)
      } else {
        if (howmany == 3) {
          m <- matrix(c(1, 2, 3, 0), byrow = TRUE, nrow = 2)
        } else {
          m <- matrix(c(1, 3, 2, 4), nrow = 2)
        }
      }
    }
  }
  if (!ggp) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(m)
    par(col = 1, las = 1, mar = c(4.5, 5, 2, 1), oma = c(0, 0, 1, 0), pch = 16,
        ...)
    if ("PP" %in% plots) {
      plot(1 - survTim, theorPP, col = colour[0 %% nCol + 1],
           xlab = "Empirical cumulative distribution",
           ylab = "Theoretical cumulative distribution",
           main = "P-P plot")
      abline(0, 1)
    }
    if ("QQ" %in% plots) {
      plot(tim, theorQQ, col = colour[1 %% nCol + 1],
           xlab = "Sample quantiles",
           ylab = "Theoretical quantiles",
           main = "Q-Q plot")
      abline(0, 1)
    }
    if ("SP" %in% plots) {
      plot(2 / pi * asin(sqrt(1 - survTim)), 2 / pi * asin(sqrt(theorPP)),
           col = colour[2 %% nCol + 1],
           xlab = expression(bold(2 / pi %*% arcsin(hat(F)[n](t)^{1/2}))),
           ylab = expression(bold(2 / pi %*% arcsin(hat(F)[0](t)^{1/2}))),
           main = "SP plot")
      abline(0, 1)
    }
    if ("ER" %in% plots) {
      plot(uPoint, uEstim, col = colour[3 %% nCol + 1],
           xlab = expression(bold(hat(F)[u](t))),
           ylab = expression(bold(hat(F)[u](paste(hat(F)[0]^{-1})(hat(F)(t))))),
           main = "ER plot")
      abline(0, 1)
    }
    if (mtitle) {
      title(paste("Probability plots for", outp$Distribution, "distribution"),
            outer = TRUE, cex.main = 1.5)
    }
  } else {
    ERx <- ERy <- PPx <- PPy <- QQx <- QQy <- SPx <- SPy <- NULL
    ggdat <- data.frame(PPx = 1 - survTim, PPy = theorPP,
                        QQx = tim, QQy = theorQQ,
                        SPx = 2 / pi * asin(sqrt(1 - survTim)),
                        SPy = 2 / pi * asin(sqrt(theorPP)),
                        ERx = uPoint, ERy = uEstim)
    PP <- ggplot(data = ggdat, aes(x = PPx, y = PPy)) +
      geom_point(colour = colour[0 %% nCol + 1]) +
      xlab("Empirical cumulative distribution") +
      ylab("Theoretical cumulative distribution") +
      geom_abline(intercept = 0) +
      annotate("text", label = "P-P plot", x = Inf, y = -Inf, hjust = 1,
               vjust = -1, size = 6, fontface = "bold")
    QQ <- ggplot(data = ggdat, aes(x = QQx, y = QQy)) +
      geom_point(colour = colour[1 %% nCol + 1]) +
      xlab("Sample quantiles") +
      ylab("Theoretical quantiles") +
      geom_abline(intercept  =  0) +
      annotate("text", label = "Q-Q plot", x = Inf, y = -Inf, hjust = 1,
               vjust = -1, size = 6, fontface = "bold")
    SP <- ggplot(data = ggdat, aes(x = SPx, y = SPy)) +
      geom_point(colour = colour[2 %% nCol + 1]) +
      xlab(expression(bold(2 / pi %*% arcsin(hat(F)[n](t)^{1/2})))) +
      ylab(expression(bold(2 / pi %*% arcsin(hat(F)[0](t)^{1/2})))) +
      geom_abline(intercept = 0) +
      annotate("text", label = "SP plot", x = Inf, y = -Inf, hjust = 1,
               vjust = -1, size = 6, fontface = "bold")
    ER <- ggplot(data = ggdat, aes(x = ERx, y = ERy)) +
      geom_point(colour = colour[3 %% nCol + 1]) +
      xlab(expression(hat(F)[u](t))) +
      ylab(expression(hat(F)[u](paste(hat(F)[0]^{-1})(hat(F)(t))))) +
      geom_abline(intercept = 0) +
      annotate("text", label = "ER plot", x = Inf, y = -Inf, hjust = 1,
               vjust = -1, size = 6, fontface = "bold")
    plolis <- list(PP, QQ, SP, ER)[c("PP", "QQ", "SP", "ER") %in% plots]
    if (mtitle) {
      grid.arrange(grobs = plolis, layout_matrix = m,
                   top = textGrob(paste("Probability plots for",
                                        outp$Distribution, "distribution"),
                   gp = gpar(fontsize = 15, font = 2)))
    } else {
      grid.arrange(grobs = plolis, layout_matrix = m)
    }
  }
  if (prnt) {
    cat("Distribution:", outp$Distribution, "\n")
    cat(" ", "\n")
    if (!all(sapply(params0, is.null))) {
      cat("Parameters used in probability plots:\n")
      for (dist in distr) {
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(outp$Parameters)) {
            cat("Location:", outp$Parameters[1], "\n")
          }
          if ("shape" %in% names(outp$Parameters)) {
            cat("   Shape:", outp$Parameters[1], "\n")
          }
          cat("   Scale:", outp$Parameters[2], "\n")
        } else if (dist == "exponential") {
          cat("   Scale:",  outp$Parameters, "\n")
        } else {
          cat("  Shape1:", outp$Parameters[1], "\n")
          cat("  Shape2:", outp$Parameters[2], "\n")
          cat("  Domain:", outp$interval.domain[1], "-",
                           outp$interval.domain[2], "\n")
        }
      }
      cat(" ", "\n")
    }
    cat("Parameter estimates:\n")
    for (dist in distr) {
      if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                      "loglogistic")) {
        if ("location" %in% names(outp$Estimates)) {
          cat("Location:", round(outp$Estimates[1], degs), "\n")
        }
        if ("shape" %in% names(outp$Estimates)) {
          cat("   Shape:", round(outp$Estimates[1], degs), "\n")
        }
        cat("   Scale:", round(outp$Estimates[2], degs), "\n")
      } else if (dist == "exponential") {
        cat("   Scale:",  round(outp$Estimates,degs), "\n")
      } else {
        cat("  Shape1:", round(outp$Estimates[1], degs), "\n")
        cat("  Shape2:", round(outp$Estimates[2], degs), "\n")
        cat("  Domain:", round(outp$interval.domain[1], degs), "-",
                         round(outp$interval.domain[2], degs), "\n")
      }
    }
  }
}
