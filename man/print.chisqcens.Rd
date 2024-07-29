\name{print.chisqcens}
\alias{print.chisqcens}
\title{Printing method for \code{chisqcens} object.}
\description{Printing method for \code{chisqcens} object.}
\usage{
\method{print}{chisqcens}(x, prnt = TRUE, outp = c("list", "table"), ...)
}
\arguments{
  \item{x}{An object of class \code{chisqcens}.}
  \item{prnt}{Logical to indicate if the estimations of the chi-squared statistic and p-value should be printed. Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
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

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}
\examples{
# List output
set.seed(123)
chisqcens(times = rweibull(100, 12, scale = 4), M = 8, distr = "weibull",
          BS = 149)

# Table output
set.seed(123)
print(chisqcens(times = rweibull(100, 12, scale = 4), M = 8, distr = "weibull",
                BS = 99), outp = "table")
}
