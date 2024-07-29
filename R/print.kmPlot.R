print.kmPlot <- function(x, ...) {
  if (x$prnt) {
    cat("Parameter estimates\n")
    cat(" ", "\n")
    for (dist in x$distributions) {
      if (!is.null(x$params[[dist]])) {
        cat(dist, "\n", sep = "")
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(x$params[[dist]])) {
            cat("Location:", round(x$params[[dist]][1], x$degs), "\n")
          }
          if ("shape" %in% names(x$params[[dist]])) {
            cat("   Shape:", round(x$params[[dist]][1], x$degs), "\n")
          }
          cat("   Scale:", round(x$params[[dist]][2], x$degs), "\n")
        } else {
          if (dist == "exponential") {
            cat("   Scale:", round(x$params[[dist]], x$degs), "\n")
          } else {
            cat("  Shape1:", round(x$params[[dist]]$parameters[1], x$degs), "\n")
            cat("  Shape2:", round(x$params[[dist]]$parameters[2], x$degs), "\n")
            cat("  Domain:", round(x$params[[dist]]$domain[1], x$degs), "-",
                round(x$params[[dist]]$domain[2],x$degs), "\n")
          }
        }
        cat("\n")
      }
    }
  }
  if ("beta" %in% x$distributions) {
    x$params$beta <- x$params$beta$parameters
  }
  invisible(lapply(x$params, round, x$degs))
}
