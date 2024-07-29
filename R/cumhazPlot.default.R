cumhazPlot.default <- function(times, cens = rep(1, length(times)),
                               distr = "all6", colour = 1, betaLimits = c(0, 1),
                               igumb = c(10, 10), ggp = FALSE, m = NULL,
                               prnt = TRUE, degs = 3, ...) {
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
  if (length(distr) == 1 && distr == "all6") {
    distributions <- c("weibull", "loglogistic", "lognormal", "gumbel",
                       "logistic", "normal")
  } else {
    distributions <- match.arg(distr, c("exponential", "weibull", "loglogistic",
                                        "lognormal", "gumbel", "logistic",
                                        "normal", "beta"), several.ok = TRUE)
  }
  if ("beta" %in% distributions && any(times < betaLimits[1] |
                                       times > betaLimits[2])) {
    warning("Beta distribution is ignored because of out-of-bounds values.
  Try with 'betaLimits = c(", pmax(0, min(times) - 1), ", ",
            ceiling(max(times) + 1), ")'.",
            immediate. = TRUE)
    cat("\n")
    distributions <- distributions[distributions != "beta"]
  }
  if (length(distributions) == 0) {
    stop("No distribution to be fitted!")
  }
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  survNA <- survfit(Surv(times, cens) ~ 1, stype = 2, ctype = 1)
  Haz <- -log(summary(survNA)$surv)
  tim <- summary(survNA)$time
  genlis <- vector("list", length(distributions))
  names(genlis) <- distributions
  for (v in c("params", "xscale", "yscale", "xlabs", "ylabs", "regline",
              "titl")) {
    assign(v, genlis)
  }
  xla = "Time"
  xlogla = "Log(Time)"
  if ("exponential" %in% distributions) {
    tryCatch({
      fitExpo <- fitdistcens(dd, "exp")
      rateExpo <- fitExpo$estimate[1]
      params$exponential <- 1 / rateExpo
      names(params$exponential) <- "scale"
      xscale$exponential <- tim
      yscale$exponential <- Haz
      xlabs$exponential <- xla
      ylabs$exponential <- expression(bolditalic(hat(Lambda)(Time)))
      regline$exponential <- rateExpo * tim
      titl$exponential <- "Exponential"
    }, error = function(e) e)
  }
  if ("gumbel" %in% distributions) {
    tryCatch({
      fitGumb <- fitdistcens(dd, "gumbel", start = list(alpha = 0, scale = 3))
      locGumb <- fitGumb$estimate[1]
      scaleGumb <- fitGumb$estimate[2]
      params$gumbel <- c(locGumb, scaleGumb)
      names(params$gumbel) <- c("location", "scale")
      xscale$gumbel <- tim
      yscale$gumbel <- -log(-log(1 - exp(-Haz)))
      xlabs$gumbel <- xla
      ylabs$gumbel <- expression(bolditalic(-log(-log(1 - exp(-hat(Lambda)(Time))))))
      regline$gumbel <- (tim - locGumb) / scaleGumb
      titl$gumbel <- "Gumbel"
    }, error = function(e) e)
  }
  if ("weibull" %in% distributions) {
    tryCatch({
      fitWei <- fitdistcens(dd, "weibull")
      shapeWei <- fitWei$estimate[1]
      scaleWei <- fitWei$estimate[2]
      params$weibull <- c(shapeWei, scaleWei)
      xscale$weibull <- log(tim)
      yscale$weibull <- log(Haz)
      xlabs$weibull <- xlogla
      ylabs$weibull <- expression(bolditalic(log(hat(Lambda)(Time))))
      regline$weibull <- shapeWei * (-log(scaleWei) + log(tim))
      titl$weibull <- "Weibull"
    }, error = function(e) e)
  }
  if ("normal" %in% distributions) {
    tryCatch({
      fitNorm <- fitdistcens(dd, "norm")
      locNorm <- fitNorm$estimate[1]
      scaleNorm <- fitNorm$estimate[2]
      params$normal <- fitNorm$estimate
      names(params$normal) <- c("location", "scale")
      xscale$normal <- tim
      yscale$normal <- qnorm(1 - exp(-Haz))
      xlabs$normal <- xla
      ylabs$normal <- expression(bolditalic(paste(Phi^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$normal <- (tim - locNorm) / scaleNorm
      titl$normal <- "Normal"
    }, error = function(e) e)
  }
  if ("logistic" %in% distributions) {
    tryCatch({
      fitLog <- fitdistcens(dd, "logis")
      locLogis <- fitLog$estimate[1]
      scaleLogis <- fitLog$estimate[2]
      params$logistic <- fitLog$estimate
      xscale$logistic <- tim
      yscale$logistic <- log(exp(Haz) - 1)
      xlabs$logistic <- xla
      ylabs$logistic <- expression(bolditalic(log(exp(hat(Lambda)(Time)) - 1)))
      regline$logistic <- (tim - locLogis) / scaleLogis
      titl$logistic <- "Logistic"
    }, error = function(e) e)
  }
  if ("lognormal" %in% distributions) {
    tryCatch({
      fitLnorm <- fitdistcens(dd, "lnorm")
      locLnorm <- fitLnorm$estimate[1]
      scaleLnorm <- fitLnorm$estimate[2]
      params$lognormal <- c(locLnorm, scaleLnorm)
      names(params$lognormal) <- c("location", "scale")
      xscale$lognormal <- log(tim)
      yscale$lognormal <- qnorm(1 - exp(-Haz))
      xlabs$lognormal <- xlogla
      ylabs$lognormal <- expression(bolditalic(paste(Phi^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$lognormal <- (log(tim) - locLnorm)/scaleLnorm
      titl$lognormal <- "Lognormal"
    }, error = function(e) e)
  }
  if ("loglogistic" %in% distributions) {
    tryCatch({
      fitLoglog <- fitdistcens(dd, "llogis")
      shapeLoglogis <- fitLoglog$estimate[1]
      scaleLoglogis <- fitLoglog$estimate[2]
      params$loglogistic <- c(shapeLoglogis, scaleLoglogis)
      xscale$loglogistic <- log(tim)
      yscale$loglogistic <- log(exp(Haz) - 1)
      xlabs$loglogistic <- xlogla
      ylabs$loglogistic <- expression(bolditalic(log(exp(hat(Lambda)(Time) - 1))))
      regline$loglogistic <- shapeLoglogis * (log(tim) - log(scaleLoglogis))
      titl$loglogistic <- "Log-logistic"
    }, error = function(e) e)
  }
  if ("beta" %in% distributions) {
    tryCatch({
      aBeta <- betaLimits[1]
      bBeta <- betaLimits[2]
      fitBeta <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      shape1Beta <- fitBeta$estimate[1]
      shape2Beta <- fitBeta$estimate[2]
      params$beta <- list(parameters = c(shape1Beta, shape2Beta),
                          domain = betaLimits)
      xscale$beta <- tim
      yscale$beta <- qbeta(1 - exp(-Haz), shape1Beta, shape2Beta)
      xlabs$beta <- xla
      ylabs$beta <- expression(bolditalic(paste(F^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$beta <- (tim - aBeta) / (bBeta - aBeta)
      titl$beta <- "Beta"
    }, error = function(e) e)
  }
  output <- list(times = times, cens = cens, distributions = distributions,
                 params = params, xscale = xscale, yscale = yscale,
                 xlabs = xlabs, ylabs = ylabs, regline = regline, titl = titl,
                 colour = colour, ggp = ggp, m = m, betaLimits = betaLimits,
                 igumb = igumb, degs = degs, prnt = prnt)
  class(output) <- "cumhazPlot"
  plot(output)
  print(output)
}
