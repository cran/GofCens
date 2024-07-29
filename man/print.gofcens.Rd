\name{print.gofcens}
\alias{print.gofcens}
\title{Printing method for \code{gofcens} object.}
\description{Printing method for \code{gofcens} object.}
\usage{
\method{print}{gofcens}(x, prnt = TRUE, outp = c("list", "table"), ...)
}
\arguments{
  \item{x}{An object of class \code{gofcens}.}
  \item{prnt}{Logical to indicate if the values of the Kolmogovor-Smirnov, Cramér-von Mises,
              and Anderson-Darling test statistics along with the p-values should be printed. 
              Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Null hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Tests statistics}}{Values of the Kolmogovor-Smirnov, Cramér-von Mises,
        and Anderson-Darling test statistics.}
  \item{\code{p-value}}{P-values associated with each test.}
  \item{\code{Parameters}}{List with the maximum likelihood estimates of the parameters 
                           of the distribution under study.}
The list is also returned invisibly.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}
\examples{
\dontrun{
# List output
set.seed(123)
gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
        BS = 149)
  
# Table output
set.seed(123)
print(gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
              BS = 149), outp = "table")
}
}
