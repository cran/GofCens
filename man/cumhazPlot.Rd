\name{cumhazPlot}
\alias{cumhazPlot}
\title{Cumulative hazard plots to check the goodness of fit of parametric models}
\description{
  Function \code{cumhazPlot} uses the cumulative hazard plot to check if a certain distribution 
  is an appropiate choice for the data.
}

\usage{
cumhazPlot(times, cens = rep(1, length(times)), distr = "all6", colour = 1, 
           betaLimits = c(0, 1), igumb = c(10, 10), ggp = FALSE, m = NULL, 
           prnt = TRUE, degs = 3, ...)
}

\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not 
              provided, all times are assumed to be exact.}
  \item{distr}{A string specifying the names of the distributions to be studied. 
               The possible distributions are the exponential (\code{"exponential"}), 
               the Weibull (\code{"weibull"}), the Gumbel (\code{"gumbel"}), 
               the normal (\code{"normal"}), the lognormal (\code{"lognormal"}),
               the logistic (\code{"logistic"}), the loglogistic (\code{"loglogistic"}),
               and the beta (\code{"beta"}) distribution. By default, \code{distr} 
               is set to \code{"all6"}, which means that the cumulative hazard 
               plots are drawn for the Weibull, loglogistic, lognormal, Gumbel, 
               logistic, and normal distributions.} 
  \item{colour}{Colour of the points. Default colour: black.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the 
                    Beta distribution. This argument is only required, if the beta
                    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation 
               of the Gumbel distribution parameters.}
  \item{ggp}{Logical to use or not the \pkg{ggplot2} package to draw the plots. 
             Default is \code{FALSE}.}
  \item{m}{Optional layout for the plots to be displayed.}
  \item{prnt}{Logical to indicate if the maximum likelihood estimates of the 
              parameters of all distributions considered should be printed. 
              Default is \code{TRUE}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric 
              results of the output.}
  \item{\dots}{Optional arguments for function \code{par}, if \code{ggplo = FALSE}.}
}

\details{The cumulative hazard plot is based on transforming the cumulative 
         hazard function \eqn{\Lambda} in such a way that it becomes linear in \eqn{t} 
         or \eqn{\log(t)}. This transformation is specific for each distribution. 
         The function uses the data to compute the Nelson-Aalen estimator of the 
         cumulative hazard function, \eqn{\widehat{\Lambda}}, and the 
         maximum likelihood estimators of the parameters of the theoretical 
         distribution under study. If the distribution fits the data, the plot is 
         expected to be a straight line.
         
         The parameter estimation is acomplished with the \code{fitdistcens} 
         function of the \pkg{fitdistrplus} package.}

\value{If \code{prnt = TRUE}:
  \item{Parameter estimates}{A list with the maximum likelihood estimates of the 
                             parameters of all distributions considered.}}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez.}

\examples{
# Complete data and default distributions
set.seed(123)
x <- rlogis(1000, 50, 5)
cumhazPlot(x, lwd = 2)

# Censored data comparing three distributions
data(nba)
cumhazPlot(nba$survtime, nba$cens, distr = c("expo", "normal", "gumbel"))
}
