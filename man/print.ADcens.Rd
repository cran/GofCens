\name{print.ADcens}
\alias{print.ADcens}
\title{Printing method for \code{ADcens} object.}
\description{Printing method for \code{ADcens} object.}
\usage{
\method{print}{ADcens}(x, prnt = TRUE, outp = c("list", "table"), ...)
}
\arguments{
  \item{x}{An object of class \code{ADcens}.}
  \item{prnt}{Logical to indicate if the estimations of the Anderson-Darling statistic and p-value should be printed. Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Null hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{AD}}{Value of the Anderson-Darling statistic.}
  \item{\code{p-value}}{Estimated p-value.}
  \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
    of the distribution under study.}
The list is also returned invisibly.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}
\examples{
# List output
set.seed(123)
ADcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
       BS = 149)
  
# Table output
set.seed(123)
print(ADcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
             BS = 99), outp = "table")
}
