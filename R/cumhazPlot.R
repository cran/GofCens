cumhazPlot <-
function(times, cens = rep(1, length(times)),
                       distr = "all6", colour = 1, betaLimits = c(0, 1),
                       igumb = c(10, 10), ggplo = FALSE, m = NULL,
                       prnt = TRUE, decdig = 7, ...) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  if (!is.logical(ggplo) || !is.logical(prnt)) {
    stop("ggplo and prnt must be logicals!")
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
    warning("Beta distributions is ignored because of out-of-bounds values.",
            immediate. = TRUE)
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
      params$exponential <- rateExpo
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
      params$gumbel <- c("location" = locGumb, scaleGumb)
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
      params$beta <- list(param = c(shape1Beta, shape2Beta), domain = betaLimits)
      xscale$beta <- tim
      yscale$beta <- qbeta(1 - exp(-Haz), shape1Beta, shape2Beta)
      xlabs$beta <- xla
      ylabs$beta <- expression(bolditalic(paste(F^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$beta <- (tim - aBeta) / (bBeta - aBeta)
      titl$beta <- "Beta"
    }, error = function(e) e)
  }
  old <- options(digits = 2)
  on.exit(options(old))
  options(digits = decdig)
  if (is.null(m)) {
    nplots = length(distributions)
    nro <- ifelse(nplots %in% 1:3, 1, ifelse(nplots %in% 4:6, 2, 3))
    m <- matrix(1:nplots, byrow = TRUE, nrow = nro)
  }
  if (!ggplo) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(m)
    par(pch = 16, las = 1, mar = c(4, 4.5, 2, 1), font.lab = 3, ...)
    for (i in distributions) {
      if (!is.null(xscale[[i]])) {
        plot(xscale[[i]], yscale[[i]], col = colour, xlab = xlabs[[i]],
             ylab = ylabs[[i]], main = titl[[i]])
        lines(xscale[[i]], regline[[i]])
      }
    }
  } else {
    plolis <- vector("list", length(distributions))
    names(plolis) <- distributions
    for (i in distributions) {
      if (!is.null(xscale[[i]])) {
        plolis[[i]] <- local({
          i <- i
          p1 <- ggplot(mapping = aes(x = xscale[[i]], y = yscale[[i]])) +
            geom_point(colour = colour) +
            geom_segment(aes(x = xscale[[i]][1], y = regline[[i]][1],
                             xend = rev(xscale[[i]])[1], yend = rev(regline[[i]])[1])) +
            labs(title = titl[[i]], size = 6, fontface = "bold") +
            xlab(xlabs[[i]]) +
            ylab(ylabs[[i]])
        })
      }
    }
    grid.arrange(grobs = plolis, layout_matrix = m)
  }
  if (prnt) {
    cat("\nParameter estimates", fill = TRUE)
    cat(rep("=", nchar("Parameter estimates")), sep = "", fill = TRUE)
    return(params)
  }
}
