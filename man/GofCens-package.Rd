\name{GofCens-package}
\alias{GofCens-package}
\alias{GofCens}
\docType{package}
\title{Goodness-of-Fit Methods for Complete and Right-Censored Data.}
\description{This package provides both graphical tools and goodness-of-fit tests for analyzing complete and right-censored data. It includes:
    \enumerate{
      \item Kolmogorov-Smirnov, Cramér-von Mises, and Anderson-Darling tests,
      which utilize the empirical distribution function for complete data and are extended to handle right-censored data.

      \item Generalized chi-squared-type test, which is based on the squared
      differences between observed and expected counts using random
      cells with right-censored data.

      \item Graphical tools, such as probability and cumulative hazard plots, to help guide decisions about the most appropriate parametric model for the data.
    }
}

\details{
  The \code{GofCens} package can be used to assess the goodness of fit for the following eight distributions. The list below displays the parameterizations of their survival functions.

  \enumerate{
    \item Exponential Distribution [Exp\eqn{(\beta)}] \deqn{S(t)=e^{-\frac{t}{\beta}}}
    \item Weibull Distribution [Wei(\eqn{\alpha,\,\beta})]  \deqn{S(t)=e^{-(\frac{t}{\beta})^\alpha}}
    \item Gumbel Distribution [Gum(\eqn{\mu,\,\beta})]  \deqn{S(t)=1 - e^{-e^{-\frac{t-\mu}{\beta}}}}
    \item Log-Logistic Distribution [LLogis(\eqn{\alpha, \beta})]  \deqn{S(t)=\frac{1}{1 + \left(\frac{t}{\beta}\right)^\alpha}}
    \item Logistic Distribution [Logis(\eqn{\mu,\beta})]  \deqn{S(t)=\frac{e^{-\frac{t -\mu}{\beta}}}{1 + e^{-\frac{t - \mu}{\beta}}}}
    \item Log-Normal Distribution [LN(\eqn{\mu,\beta})]  \deqn{S(t)=\int_{\frac{\log t - \mu}{\beta}}^\infty \!\frac{1}{\sqrt{2 \pi}}}
    \item Normal Distribution [N(\eqn{\mu,\beta})]  \deqn{S(t)=\int_t^\infty \! \frac{1}{\beta\sqrt{2\pi}}e^{-\frac{(x - \mu)^2}{2 \beta^2}} dx}
    \item 4-Param. Beta Distribution [Beta(\eqn{\alpha, \gamma, a, b})] \deqn{S(t)=1 - \frac{B_{(\alpha, \gamma, a, b)}(t)}{B(\alpha, \gamma)}}
  }

  The parameters of the theoretical distribution can be set manually using the \code{params0} argument in each function. 
  In this case, the correspondences are as follows: \eqn{\alpha} represents the \code{shape}, \eqn{\gamma} the \code{shape2}, 
  \eqn{\mu} the \code{location}, and \eqn{\beta} the \code{scale} parameter.


  \tabular{ll}{
    Package: \tab GofCens\cr
    Type: \tab Package\cr
    Version: \tab 1.2.1\cr
    Date: \tab 2024-11-9\cr
    License: \tab GPL (>= 2)\cr}
}
\author{Klaus Langohr, Mireia Besalú, Matilde Francisco, Arnau Garcia, Guadalupe Gómez

  Maintainer: Klaus Langohr <klaus.langohr@upc.edu>}
