kmPlot <-
  function(times, cens = rep(1, length(times)),
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
      warning("Beta distributions is ignored because of out-of-bounds values.
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
    if (is.null(m)) {
      nplots = length(distributions)
      nro <- ifelse(nplots %in% 1:3, 1, ifelse(nplots %in% 4:6, 2, 3))
      m <- matrix(1:nplots, byrow = TRUE, nrow = nro)
    }
    if (!ggp) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      layout(m)
      par(pch = 16, las = 1, mar = c(4, 4.5, 2, 1), font.lab = 4, ...)
      for (i in distributions) {
        plot(survKM, col = colour, xlab = "Time",
             ylab = expression(bolditalic(hat(S)(t))), main = titl[[i]])
        lines(sqtime, srvline[[i]])
      }
    } else {
      plolis <- vector("list", length(distributions))
      names(plolis) <- distributions
      for (i in distributions) {
        if (!is.null(srvline[[i]])) {
          plolis[[i]] <- local({
            i <- i
            tmpdat <- data.frame(x = sqtime, y = srvline[[i]])
            p1 <- ggsurvplot(survKM, data = data.frame(times, cens),
                             ggtheme = theme_minimal(),
                             xlab = expression(bolditalic(Time)),
                             ylab = expression(bolditalic(hat(S)(t))),
                             censor = FALSE, legend = "none",
                             title = titl[[i]],
                             font.main = c(14, "bold", "black"),
                             palette = colour)$plot +
              geom_point(aes(tmpdat$x, tmpdat$y), size = 1, data = tmpdat) +
              geom_line(aes( tmpdat$x, tmpdat$y), data = tmpdat)
          })
        }
      }
      suppressWarnings(grid.arrange(grobs = plolis, layout_matrix = m))
    }
    if (prnt) {
      cat("Parameter estimates\n")
      cat(" ", "\n")
      for (dist in distributions) {
        if (!is.null(params[[dist]])) {
          cat(dist, "\n", sep = "")
          if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                          "loglogistic")) {
            if ("location" %in% names(params[[dist]])) {
              cat("Location:", round(params[[dist]][1], degs), "\n")
            }
            if ("shape" %in% names(params[[dist]])) {
              cat("   Shape:", round(params[[dist]][1], degs), "\n")
            }
            cat("   Scale:", round( params[[dist]][2], degs), "\n")
          } else if(dist == "exponential") {
            cat("   Scale:", round(params[[dist]], digits = degs), "\n")
          } else if (dist == "beta") {
            cat("   Shape1:", round(params[[dist]]$parameters[1], degs), "\n")
            cat("   Shape2:", round(params[[dist]]$parameters[2], degs), "\n")
            cat("   Domain:", round(params[[dist]]$domain[1],degs), "-",
                             round(params[[dist]]$domain[2],degs), "\n")
          }
          cat("\n")
        }
      }
   }
}
