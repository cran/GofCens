\name{chisqcens2}
\alias{chisqcens2}
\title{General chi-squared test for right-censored data.}
\description{
\code{chisqcens2} computes the general chi-squared statistic for right-censored 
data introduced by Kim (1993). The p-value is computed using bootstrapping. 
}
\usage{
chisqcens2(times, cens, M, 
           distrData = c("weibull", "lognormal", "loglogistic"), 
           distrCens = c("weibull", "lognormal", "loglogistic", "uniform"), 
           BS = 1000, degs = 4, 
           params = list(shape = NULL, location = NULL, scale = NULL))
}
\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time).}
  \item{M}{Number indicating the number of cells that will be considered.}
  \item{distrData}{A string specifying the name of the distribution to be studied. 
                   The possible distributions are the Weibull (\code{"weibull"}), 
                   the lognormal (\code{"lognormal"}), and the 
                   loglogistic (\code{"loglogistic"}) distribution.}
  \item{distrCens}{A string specifying the name of the distribution of the 
                   censoring times. The possible distributions are the Weibull
                   (\code{"weibull"}), the lognormal (\code{"lognormal"}), the
                   loglogistic (\code{"loglogistic"}), and the uniform 
                   (\code{"uniform"}) distribution.}
  \item{BS}{Number of randomly censored samples under the null distribution. 
            Default value: \code{BS = 1000}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric 
              results of the output.}
  \item{params}{List specifying the parameters of the theoretical distribution. 
                By default, parameters are set to \code{NULL} and estimated with 
                the maximum likelihood method. This argument is only considered, 
                if all parameters of the studied distribution are specified.}
}
\details{The function implements the test introduced by Kim (1993). 

         Different from the function \code{chisqcens1}, this function does not only
         provide the value of the test statistic, but also a p-value. For this
         purpose, random censored samples are generated under the null 
         distribution and a given distribution of the censoring times, and the
         value of the test statistic is obtained for each sample. The empirical
         distribution of the test statistics obtained is used to determine the
         p-value.

         The random censored samples are generated with the \code{simple.surv.sim} 
         function of the \pkg{survsim} package.}
\value{
A list containing the following components
 \item{\code{Statistic}}{Value of the test statistic.}
 \item{\code{pvalue}}{p-value computed using bootstrapping.}
 \item{\code{distrData}}{Null distribution.}
 \item{\code{distrCens}}{Distribution of the censoring times.}
 \item{\code{Parameters}}{The values of the parameters of the null distribution. 
                          If the user has set the parameters manually, these will 
                          be the returned parameters, otherwise the maximum 
                          likelihood estimates are returned.}
 \item{\code{CellNumber}}{Vector with two values: the original cell number 
                          introduced by the user and the final cell number used.}
}
\references{
J. H. Kim. \emph{Chi-Square Goodness-of-Fit Tests for Randomly Censored Data.} 
In: The Annals of Statistics, 21 (3) (1993), 1621-1639.
}
\author{K. Langohr, M. Besalú, G. Gómez.}
\section{Warning }{If the amount of data is large, the execution time of the 
                   function can be elevated. The parameter \code{BS} can 
                   limit the number of random censored samples generated and 
                   reduce the execution time. Also, notice that this function 
                   only works with right-censored data.} 
\seealso{
\link{chisqcens1} for the computation of the test statistic of different distributions.
}
\examples{
\dontrun{
set.seed(123)
library(survsim)
n <- 50
datos <- simple.surv.sim(n, Inf,  dist.ev = "lnorm", 3, 2, 
                         dist.cens = "lnorm", 3, 2)
chisqcens2(datos$stop, datos$status, M = 8, "lognormal", "lognormal")
chisqcens2(datos$stop, datos$status, M = 8, "lognormal", "unif")

data(nba)
chisqcens2(nba$survtime, nba$cens, 10, "loglog", BS = 100)
}
}
