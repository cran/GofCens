\name{gofcens}
\alias{gofcens}
\title{Kolmogorov-Smirnov, Crámer-von Mises, and Anderson-Darling statistics for 
       complete and right-censored data}
\description{
\code{gofcens} computes the Kolmogorov-Smirnov, Crámer-von Mises, and 
Anderson-Darling statistics for complete and right-censored data against 
eight possible distributions.
}
\usage{
gofcens(times, cens = rep(1, length(times)), 
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
\details{Fleming et al. (1980) proposed a modified Kolmogorov-Smirnov test to be 
         used with right-censored data. Koziol and Green (1976) proposed a Crámer-von 
         Mises statistic for randomly censored data. This function reproduces 
         this test for a given survival data and a theorical distribution. In 
         presence of ties, different authors provide slightly different definitions 
         of the product-limit estimator, what might provide different values of 
         the test statistic.  
         
         When dealing with complete data, we recommend the use of functions 
         \code{\link[stats]{ks.test}} of the \pkg{stats} package and 
         \code{\link[goftest]{cvm.test}} and \code{\link[goftest]{ad.test}} of 
         the \pkg{goftest} package. 
}
\value{
A list containing the following components
 \item{\code{Tests statistics}}{Values of the Kolmogovor-Smirnov, Crámer-von Mises,
       and Anderson-Darling test statistics }
 \item{\code{Distribution}}{Null distribution}
 \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
                         of the distribution under study.}
}
\references{
T. R. Fleming et al. \emph{Modified Kolmogorov-Smirnov test procedure with application 
to arbitrarily right-censored data}. In: Biometrics 36 (1980), 607-625.

J. A. Koziol and S. B. Green. \emph{A Crámer-von Mises statistic for randomly 
censored data}. In: Biometrika, 63 (3) (1976), 465-474.

A. N. Pettitt and M. A. Stephens. \emph{Modified Crámer-von Mises statistics 
for censored data.} In: Biometrika, 63 (2) (1976), 291-298.
}
\author{K. Langohr, M. Besalú, G. Gómez.}
\seealso{\link[stats]{ks.test} (Package \code{stats}), \link[goftest]{cvm.test} 
         (Package \code{goftest}), and \link[goftest]{ad.test} 
         (Package \code{goftest}) for complete data, and \link{KScens} for the 
         Kolmogorov-Smirnov test for right-censored data, which returns the p-value.
}
\examples{
# Complete data
set.seed(123)
gofcens(times = rweibull(1000, 12, scale = 4), distr = "weibull")

# Censored data
library(survival)
gofcens(aml$time, aml$status, distr = "norm")

data(nba)
gofcens(nba$survtime, nba$cens, "lognorm")
gofcens(nba$survtime, nba$cens, "norm")
}
