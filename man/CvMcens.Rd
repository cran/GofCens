\name{CvMcens}
\alias{CvMcens}
\title{Cramér-von Mises test for complete and right-censored data}
\description{
\code{CvMcens} computes the Cramér-von Mises statistic and p-value for complete 
and right-censored data against eight possible distributions.
}
\usage{
CvMcens(times, cens = rep(1, length(times)),
        distr = c("exponential", "gumbel", "weibull", "normal",
                  "lognormal", "logistic", "loglogistic", "beta"),
        betaLimits = c(0, 1), igumb = c(10, 10), degs = 3, BS = 999,
        params = list(shape = NULL, shape2 = NULL,
                      location = NULL, scale = NULL),
        prnt = TRUE, outp = "list", tol = 1e-04)
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
  \item{BS}{Number of bootstrap samples.}
  \item{params}{List specifying the parameters of the theoretical distribution. 
                By default, parameters are set to \code{NULL} and estimated with 
                the maximum likelihood method. This argument is only considered, 
                if all parameters of the studied distribution are specified.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are                        
  \code{list} and \code{table}.}
  \item{prnt}{Logical to indicate if the estimations of the Anderson-Darling statistic and p-value should be printed. Default is \code{TRUE}.}
  \item{tol}{Precision of survival times.}
}
\details{
         Koziol and Green (1976) proposed a Cramér-von 
         Mises statistic for randomly censored data. This function reproduces 
         this test for a given survival data and a theorical distribution. In 
         presence of ties, different authors provide slightly different definitions 
         of the product-limit estimator, what might provide different values of 
         the test statistic.  
         When dealing with complete data, we recommend to use the function 
         \code{cvm.test} of the \pkg{goftest} package.
         
         The parameter estimation is acomplished with the \code{fitdistcens} 
         function of the \pkg{fitdistrplus} package.
         
         The precision of the survival times is important mainly in 
         the data generation step of the bootstrap samples.
}
\value{
If \code{prnt = TRUE}, a list containing the following components:
  \item{\code{Distribution }}{Null distribution.}
  \item{\code{CvM}}{Value of Cramér-von Mises statistic.}
  \item{\code{p-value}}{Estimated p-value.}
  \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
                           of the distribution under study.}
The list is also returned invisibly.
}

\references{
J. A. Koziol and S. B. Green. \emph{A Cramér-von Mises statistic for randomly 
censored data}. In: Biometrika, 63 (3) (1976), 465-474.

A. N. Pettitt and M. A. Stephens. \emph{Modified Cramér-von Mises statistics 
for censored data.} In: Biometrika, 63 (2) (1976), 291-298.
}
\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}

\section{Warning }{If the amount of data is large, the execution time of the 
                   function can be elevated. The parameter \code{BS} can 
                   limit the number of random censored samples generated and 
                   reduce the execution time.} 

\seealso{Function \code{cvm.test} (Package \pkg{goftest}) for complete data and 
         \link{gofcens} for statistics and p-value of Kolmogorov-Smirnov, 
         Cramér von-Mises and Anderson-Darling together for right-censored data.
}

\examples{
# Complete data
set.seed(123)
CvMcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
       BS = 199)

\dontrun{
# Censored data
library(survival)
colonsamp <- colon[sample(nrow(colon), 100), ]
CvMcens(colonsamp$time, colonsamp$status, distr = "normal")
}
}
