\name{nba}
\alias{nba}
\docType{data}
\title{Survival times of former NBA players.}
\description{Survival times of former NBA players after their NBA career.}
\usage{data("nba")}
\format{
  A data frame with 3501 observations on the following 3 variables.
  \describe{
    \item{\code{id}}{Player ID}
    \item{\code{survtime}}{Time (in years) from end of NBA career until either
    death or April 15, 2014.}
    \item{\code{cens}}{Death indicator (1, exact survival time; 0, right-censored 
    survival time).}
  }
}
\details{
The survival times of former NBA players were analyzed by Martínez et al. (2019). 
}
\source{
J. A. Martínez, K. Langohr, J. Felipo, and M. Casals. \emph{Mortality of NBA 
players: Risk factors and comparison with the general US population.} In: Applied 
Sciences, 9 (3) (2019).
}
\examples{
data(nba)
cumhazPlot(nba$survtime, nba$cens)
}
\keyword{datasets}