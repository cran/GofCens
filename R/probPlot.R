probPlot <-
function(times, cens = rep(1, length(times)),
                     distr = c("exponential", "gumbel", "weibull", "normal",
                               "lognormal", "logistic", "loglogistic", "beta"),
                     plots = c("PP", "QQ", "SP", "ER"),
                     colour = c("green4", "deepskyblue4", "yellow3",
                                "mediumvioletred"), betaLimits = c(0, 1),
                     igumb = c(10, 10), mtitle = TRUE, ggp = FALSE, m = NULL,
                     prnt = TRUE, degs = 3,
                     params = list(shape = NULL, shape2 = NULL,
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
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ", ceiling(max(times) + 1), ")'.")
    stop(msg)
  }
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  survKM <- survfit(Surv(times, cens) ~ 1)
  tim <- summary(survKM)$time
  survTim <- summary(survKM)$surv
  uPointSurv <- survfit(Surv(tim, rep(1, length(tim))) ~ 1)
  uPoint <- 1 - uPointSurv$surv
  empiricF <- rbind(c(0, tim, Inf), c(0, uPoint, 1))
  uEstim <- rep(0, length(uPoint))
  if (distr == "exponential") {
    if (is.null(params$scale)) {
      muu <- unname(coefficients(survreg(Surv(times, cens) ~ 1,
                                         dist = "exponential")))
      rateExp <- exp(-muu)
    } else {
      rateExp <- 1 / params$scale
    }
    theorPP <- pexp(tim, rateExp)
    theorQQ <- qexp(1 - survTim, rateExp)
    outp <- list(Distribution = "Exponential", Parameters = 1 / rateExp)
  }
  if (distr == "gumbel") {
    if (is.null(params$location) || is.null(params$scale)) {
      fitGum <- try(suppressMessages(fitdistcens(dd, "gumbel",
                                                 start = list(alpha = igumb[1],
                                                              scale = igumb[2]))),
                    silent = TRUE)
      if (attr(fitGum, "class") == "try-error") {
        stop("Function failed to estimate the parameters. Try with other initial values.")
      }
      locGum <- unname(fitGum$estimate[1])
      scaleGum <- unname(fitGum$estimate[2])
    } else {
      locGum <- params$location
      scaleGum <- params$scale
    }
    theorPP <- pgumbel(tim, locGum, scaleGum)
    theorQQ <- qgumbel(1 - survTim, locGum, scaleGum)
    outp <- list(Distribution = "Gumbel", Parameters = c(location = locGum,
                                                         scale = scaleGum))
  }
  if (distr == "weibull") {
    if (is.null(params$shape) || is.null(params$scale)) {
      fitWei <- fitdistcens(dd, "weibull")
      shapeWei <- unname(fitWei$estimate[1])
      scaleWei <- unname(fitWei$estimate[2])
    } else {
      shapeWei <- params$shape
      scaleWei <- params$scale
    }
    theorPP <- pweibull(tim, shapeWei, scaleWei)
    theorQQ <- qweibull(1 - survTim, shapeWei, scaleWei)
    outp <- list(Distribution = "Weibull", Parameters = c(shape = shapeWei,
                                                          scale = scaleWei))
  }
  if (distr == "normal") {
    if (is.null(params$location) || is.null(params$scale)) {
      fitNorm <- fitdistcens(dd, "norm")
      locNorm <- unname(fitNorm$estimate[1])
      scaleNorm <- unname(fitNorm$estimate[2])
    } else {
      locNorm <- params$location
      scaleNorm <- params$scale
    }
    theorPP <- pnorm(tim, locNorm, scaleNorm)
    theorQQ <- qnorm(1 - survTim, locNorm, scaleNorm)
    outp <- list(Distribution = "Normal", Parameters = c(location = locNorm,
                                                         scale = scaleNorm))
  }
  if (distr == "lognormal") {
    if (is.null(params$location) || is.null(params$scale)) {
      fitLnorm <- fitdistcens(dd, "lnorm")
      locLnorm <- unname(fitLnorm$estimate[1])
      scaleLnorm <- unname(fitLnorm$estimate[2])
    } else {
      locLnorm <- params$location
      scaleLnorm <- params$scale
    }
    theorPP <- plnorm(tim, locLnorm, scaleLnorm)
    theorQQ <- qlnorm(1 - survTim, locLnorm, scaleLnorm)
    outp <- list(Distribution = "Log-normal", Parameters = c(location = locLnorm,
                                                             scale = scaleLnorm))
  }
  if (distr == "logistic") {
    if (is.null(params$location) || is.null(params$scale)) {
      fitLogis <- fitdistcens(dd, "logis")
      locLogis <- unname(fitLogis$estimate[1])
      scaleLogis <- unname(fitLogis$estimate[2])
    } else {
      locLogis <- params$location
      scaleLogis <- params$scale
    }
    theorPP <- plogis(tim, locLogis, scaleLogis)
    theorQQ <- qlogis(1 - survTim, locLogis, scaleLogis)
    outp <- list(Distribution = "Logistic", Parameters = c(location = locLogis,
                                                           scale = scaleLogis))
  }
  if (distr == "loglogistic") {
    if (is.null(params$shape) || is.null(params$scale)) {
      fitLoglog <- unname(survreg(Surv(times, cens) ~ 1,
                                  dist = "loglogistic")$icoef)
      shapeLoglog <- 1 / exp(fitLoglog[2])
      scaleLoglog <- exp(fitLoglog[1])
    } else {
      shapeLoglog <- params$shape
      scaleLoglog <- params$scale
    }
    theorPP <- pllogis(tim, shapeLoglog, scale = scaleLoglog)
    theorQQ <- qllogis(1 - survTim, shapeLoglog, scale = scaleLoglog)
    outp <- list(Distribution = "Log-logistic", Parameters = c(shape = shapeLoglog,
                                                               scale = scaleLoglog))
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    if (is.null(params$shape) || is.null(params$shape2)) {
      fitBeta <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      shape1Beta <- unname(fitBeta$estimate[1])
      shape2Beta <- unname(fitBeta$estimate[2])
    } else {
      shape1Beta <- params$shape
      shape2Beta <- params$shape2
    }
    theorPP <- pbeta((tim - aBeta)/(bBeta - aBeta), shape1Beta, shape2Beta)
    theorQQ <- qbeta((1 - survTim), shape1Beta, shape2Beta) * (bBeta - aBeta)
               + aBeta
    outp <- list(Distribution = "Beta", Parameters = c(shape1 = shape1Beta,
                                                       shape2 = shape2Beta),
                 interval.domain = betaLimits)
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
      title(paste("Probability plots for a", outp$Distribution, "distribution"),
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
                   top = textGrob(paste("Probability plots for a",
                                        outp$Distribution, "distribution"),
                   gp = gpar(fontsize = 15, font = 2)))
    } else {
      grid.arrange(grobs = plolis, layout_matrix = m)
    }
  }
  if (prnt) {
    cat("Parameter Estimates\n")
    cat(" ", "\n")
    for (dist in distr) {
      if (!is.null(outp$Parameters)) {
        cat(dist, ":\n", sep = "")
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(outp$Parameters)) {
            cat("Location:", round(outp$Parameters[1], degs), "\n")
          }
          cat("   Scale:", round(outp$Parameters[2], degs), "\n")
          if ("shape" %in% names(outp$Parameters)) {
            cat("   Shape:", round(outp$Parameters[1], degs), "\n")
          }
        } else if (dist == "exponential") {
          cat("Scale:",  round(outp$Parameters,degs), "\n")
        } else {
          cat("  Shape1:", round(outp$Parameters[1], degs), "\n")
          cat("  Shape2:", round(outp$Parameters[2], degs), "\n")
          cat("  Domain:", round(outp$interval.domain[1], degs), "-",
                           round(outp$interval.domain[2], degs), "\n")
        }
        cat("\n")
      }
    }
  }
}
