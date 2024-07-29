\name{print.KScens}
\alias{print.KScens}
\title{Printing method for \code{KScens} object.}
\description{Printing method for \code{KScens} object.}
\usage{
\method{print}{KScens}(x, prnt = TRUE, outp = c("list", "table"), ...)
}
\arguments{
  \item{x}{An object of class \code{KScens}.}
  \item{prnt}{Logical to indicate if the estimations of the Kolmogorov-Smirnov statistic and p-value should be printed. Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Null hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{A}}{Value of the modified Kolmogorov-Smirnov statistic.}
  \item{\code{p-value}}{P-value.}
  \item{\code{F(y_m)}}{Estimation of the image of the last recorded time.}
  \item{\code{y_m}}{Last recorded time.}
  \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
                           of the distribution under study.}
The list is also returned invisibly.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}
\examples{
# List output
set.seed(123)
KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull")
  
# Table output
set.seed(123)
print(KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull"), 
      outp = "table")
}
