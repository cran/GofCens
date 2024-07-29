kmPlot.default <- function(times, cens = rep(1, length(times)),
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
    stop("Status indicator must be either 0 or 1!")
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
  survKM <- survfit(Surv(times, cens) ~ 1)
  mxtime <- max(survKM$time)
  sqtime <- seq(0, mxtime, length.out = 1001)
  genlis <- vector("list", length(distributions))
  names(genlis) <- distributions
  for (v in c("params", "srvline", "titl")) {
    assign(v, genlis)
  }
  if ("exponential" %in% distributions) {
    tryCatch({
      fitExpo <- fitdistcens(dd, "exp")
      rateExpo <- fitExpo$estimate[1]
      params$exponential <- 1 / rateExpo
      names(params$exponential) <- "scale"
      srvline$exponential <- 1 - pexp(sqtime, rateExpo)
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
      srvline$gumbel <- 1 - pgumbel(sqtime, locGumb, scaleGumb)
      titl$gumbel <- "Gumbel"
    }, error = function(e) {warning("Problem fitting Gumbel distribution, try other initial values.", immediate. = TRUE)})
  }
  if ("weibull" %in% distributions) {
    tryCatch({
      fitWei <- fitdistcens(dd, "weibull")
      shapeWei <- fitWei$estimate[1]
      scaleWei <- fitWei$estimate[2]
      params$weibull <- c(shapeWei, scaleWei)
      srvline$weibull <- 1 - pweibull(sqtime, shapeWei, scaleWei)
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
      srvline$normal <- 1 - pnorm(sqtime, locNorm, scaleNorm)
      titl$normal <- "Normal"
    }, error = function(e) e)
  }
  if ("logistic" %in% distributions) {
    tryCatch({
      fitLog <- fitdistcens(dd, "logis")
      locLogis <- fitLog$estimate[1]
      scaleLogis <- fitLog$estimate[2]
      params$logistic <- fitLog$estimate
      srvline$logistic <- 1 - plogis(sqtime, locLogis, scaleLogis)
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
      srvline$lognormal <- 1 - plnorm(sqtime, locLnorm, scaleLnorm)
      titl$lognormal <- "Lognormal"
    }, error = function(e) e)
  }
  if ("loglogistic" %in% distributions) {
    tryCatch({
      fitLoglog <- fitdistcens(dd, "llogis")
      shapeLoglogis <- fitLoglog$estimate[1]
      scaleLoglogis <- fitLoglog$estimate[2]
      params$loglogistic <- c(shapeLoglogis, scaleLoglogis)
      srvline$loglogistic <- 1 - pllogis(sqtime, shapeLoglogis, scale = scaleLoglogis)
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
      srvline$beta <- 1 - pbeta((sqtime - aBeta)/(bBeta - aBeta), shape1Beta,
                                shape2Beta)
      titl$beta <- "Beta"
    }, error = function(e) e)
  }
  output <- list(times = times, cens = cens, distributions = distributions,
                 params = params, survKM = survKM, sqtime = sqtime,
                 srvline = srvline, titl = titl, colour = colour, ggp = ggp,
                 m = m, prnt = prnt, degs = degs)
  class(output) <- "kmPlot"
  plot(output)
  print(output)
}
