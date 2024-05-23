\name{chisqcens}
\alias{chisqcens}
\title{General chi-squared statistics for right-censored data.}
\description{
  Function \code{chisqcens} computes the general chi-squared test statistic for 
  right-censored data introduced by Kim (1993) and the respective p-value 
  using bootstrapping.}

\usage{
chisqcens(times, cens = rep(1, length(times)), M,
          distr = c("exponential", "gumbel", "weibull", "normal",
                    "lognormal", "logistic", "loglogistic", "beta"),
          betaLimits=c(0, 1), igumb = c(10, 10), degs = 3, BS = 999,
          params0 = list(shape = NULL, shape2 = NULL,
                         location = NULL, scale = NULL),
          prnt = TRUE, outp = "list", tol = 1e-04)
}

\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not 
              provided, all times are assumed to be exact.}
  \item{M}{Number indicating the number of cells that will be considered.}
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
  \item{BS}{Number of bootstrap samples.}
  \item{params0}{List specifying the parameters of the theoretical distribution. 
                 By default, parameters are set to \code{NULL} and estimated with 
                 the maximum likelihood method. This argument is only considered, 
                 if all parameters of the studied distribution are specified.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are                        
              \code{list} and \code{table}.}
  \item{prnt}{Logical to indicate if the estimations of the Anderson-Darling statistic and 
              p-value should be printed. Default is \code{TRUE}.}
  \item{tol}{Precision of survival times.}
}

\details{The function implements the test introduced by Kim (1993) and returns 
         the value of the test statistic. 
         
         The cell boundaries of the test are obtained via the quantiles, which 
         are based on the Kaplan-Meier estimate of the distribution function. 
         In the presence of right-censored data, it is possible that not all 
         quantiles are estimated, and in this case, the value of \code{M} 
         provided by the user is reduced. 

         The parameter estimation is acomplished with the \code{fitdistcens} 
         function of the \pkg{fitdistrplus} package.
         
         The precision of the survival times is important mainly in 
         the data generation step of the bootstrap samples.
}

\value{
If \code{prnt = TRUE}, a list containing the following components:
 \item{\code{Distribution}}{Null distribution.}
  \item{\code{Null hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
 \item{\code{Statistic}}{Value of the test statistic.}
 \item{\code{p-value}}{Estimated p-value.}
 \item{\code{Parameters }}{The values of the parameters of the null distribution. 
                           If the user has set the parameters manually, these will 
                           be the returned parameters, otherwise the maximum 
                           likelihood estimates are returned.}
 \item{\code{CellNumber}}{Vector with two values: the original cell number 
                          introduced by the user and the final cell number used.}
The list is also returned invisibly.
}

\references{
  J. H. Kim. \emph{Chi-Square Goodness-of-Fit Tests for Randomly Censored Data.} 
  In: The Annals of Statistics, 21 (3) (1993), 1621-1639.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}

\examples{
# Complete data
set.seed(123)
chisqcens(time = rgumbel(100, 12, scale = 4), M = 8, distr = "gumbel",
          BS = 99)

\dontrun{
# Censored data
library(survival)
colonsamp <- colon[sample(nrow(colon), 100), ]
chisqcens(colonsamp$time, colonsamp$status, M = 6, distr = "normal")
}
}
