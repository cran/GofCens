\name{KScens}
\alias{KScens}
\title{Kolmogorov-Smirnov test for complete and right-censored data}
\description{
\code{KScens} computes the Kolmogorov-Smirnov statistic and p-value for complete 
and right-censored data against eight possible distributions.
}
\usage{KScens(times, cens = rep(1, length(times)), 
       distr = c("exponential", "gumbel", "weibull", "normal", 
                 "lognormal", "logistic", "loglogistic", "beta"), 
       betaLimits = c(0, 1), igumb = c(10, 10), degs = 4, 
       params = list(shape = NULL, shape2 = NULL, location = NULL, 
                     scale = NULL))
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
  \item{params}{List specifying the parameters of the theoretical distribution. 
                By default, parameters are set to \code{NULL} and estimated with 
                the maximum likelihood method. This argument is only considered, 
                if all parameters of the studied distribution are specified.}
}
\details{Fleming et al. (1980) proposed a modified Kolmogorov-Smirnov test to use 
         with right-censored data. This function reproduces this test for a 
         given survival data and a theorical distribution. The p-value is computed 
         following the results of Koziol and Byar (1975) and the output of the 
         function follows the notation of Fleming et al. (1980). 
         
         In presence of ties, different authors provide slightly different 
         definitions of \eqn{\widehat{F}_n(t)}, with which other values of 
         the test statistic might be obtained.
         
         When dealing with complete data, we recommend to use the function 
         \code{\link[stats]{ks.test}} of the \pkg{stats} package. 
         
         The parameter estimation is acomplished with the \code{fitdistcens} 
         function of the \pkg{fitdistrplus} package.}
\value{
A list containing the following components:
\item{\code{p-value}}{Estimated p-value.}
\item{\code{A}}{Value of the modified Kolmogorov-Smirnov statistic.}
\item{\code{F(y_m)}}{Estimation of the image of the last recorded time.}
\item{\code{y_m}}{Last recorded time.}
\item{\code{Distribution}}{Null distribution.}
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
\author{K. Langohr, M. Besalú, G. Gómez.}
\seealso{\link[stats]{ks.test} (Package \code{stats}) for complete data and 
         \link{gofcens} for Crámer von-Mises and Anderson-Darling statistics for 
         right-censored data.}
\examples{
# Complete data
set.seed(123)
KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull")

# Censored data
library(survival)
KScens(aml$time, aml$status, distr = "norm")

data(nba)
KScens(nba$survtime, nba$cens, "logis", degs = 2)
KScens(nba$survtime, nba$cens, "beta", betaLimits = c(0, 70))
}
