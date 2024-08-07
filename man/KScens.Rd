\name{KScens}
\alias{KScens}
\alias{KScens.default}
\alias{KScens.formula}
\title{Kolmogorov-Smirnov test for complete and right-censored data}
\description{
  Function \code{KScens} computes the Kolmogorov-Smirnov statistic and p-value for complete 
  and right-censored data against eight possible distributions.
}

\usage{
\method{KScens}{default}(times, cens = rep(1, length(times)), 
       distr = c("exponential", "gumbel", "weibull", "normal", 
                 "lognormal", "logistic", "loglogistic", "beta"), 
       betaLimits = c(0, 1), igumb = c(10, 10), degs = 3, 
       params0 = list(shape = NULL, shape2 = NULL, location = NULL, 
                      scale = NULL), ...)
\method{KScens}{formula}(formula, data, ...)                     
}

\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not 
              provided, all times are assumed to be exact.}
  \item{distr}{A string specifying the name of the distribution to be studied. 
               The possible distributions are the exponential (\code{"exponential"}), 
               the Weibull (\code{"weibull"}), the Gumbel (\code{"gumbel"}), 
               the normal (\code{"normal"}), the lognormal (\code{"lognormal"}),
               the logistic (\code{"logistic"}), the loglogistic (\code{"loglogistic"}),
               and the beta (\code{"beta"}) distribution.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the 
                    Beta distribution. This argument is only required, if the beta
                    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation 
               of the Gumbel distribution parameters.}
  \item{degs}{Integer indicating the number of decimal places of the numeric 
              results of the output.}
  \item{params0}{List specifying the parameters of the theoretical distribution. 
                 By default, parameters are set to \code{NULL} and estimated with 
                 the maximum likelihood method. This argument is only considered, 
                 if all parameters of the studied distribution are specified.}
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{\dots}{Additional arguments.}
}

\details{Fleming et al. (1980) proposed a modified Kolmogorov-Smirnov test to use 
         with right-censored data. This function reproduces this test for a 
         given survival data and a theorical distribution. The p-value is computed 
         following the results of Koziol and Byar (1975) and the output of the 
         function follows the notation of Fleming et al. (1980). 
         
         In presence of ties, different authors provide slightly different 
         definitions of \eqn{\widehat{F}_n(t)}, with which other values of 
         the test statistic might be obtained.
         
         An alternative with complete data is the function 
         \code{\link[stats]{ks.test}} of the \pkg{stats} package. 
         
         The parameter estimation is acomplished with the \code{fitdistcens} 
         function of the \pkg{fitdistrplus} package.}

\value{
\code{KScens} returns an object of \link[base]{class} \code{"KScens"}.

An object of class \code{"KScens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Null hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{A}}{Value of the modified Kolmogorov-Smirnov statistic.}
  \item{\code{p-value}}{P-value.}
  \item{\code{F(y_m)}}{Estimation of the image of the last recorded time.}
  \item{\code{y_m}}{Last recorded time.}
  \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
                           of the distribution under study.}
}

\references{
T. R. Fleming et al. \emph{Modified Kolmogorov-Smirnov test procedure with application 
to arbitrarily right-censored data}. In: Biometrics 36 (1980), 607-625.

J.A. Koziol and P. Byar. \emph{Percentage Points of the Asymptotic Distributions 
of One and Two Sample K-S statistics for Truncated or Censored Data}. In: 
Technometrics 17 (4) (1975), 507-510.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}

\seealso{Function \link[stats]{ks.test} (Package \pkg{stats}) for complete data and 
         \link{gofcens} for statistics and p-value of Kolmogorov-Smirnov, Cramér von-Mises 
         and Anderson-Darling together for right-censored data.
}

\examples{
# Complete data
set.seed(123)
KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull")

# Censored data
KScens(Surv(time, status) ~ 1, colon, distr = "norm")

data(nba)
KScens(Surv(survtime, cens) ~ 1, nba, "logis", degs = 2)
KScens(Surv(survtime, cens) ~ 1, nba, "beta", betaLimits = c(0, 80))
}
