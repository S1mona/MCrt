% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00output.R
\name{MCcont}
\alias{MCcont}
\title{Monte Carlo Simulations}
\usage{
MCcont(x, y, ymin, ymax, test, alpha, betabarj, j, MC, method, par)
}
\arguments{
\item{x}{matrix of exogenous variables}

\item{y}{vector of endogenous variables}

\item{ymin}{lower bound of endogenous variable}

\item{ymax}{upper bound of endogenous variable}

\item{test}{specifies test either t - or robust t - test ("t" or "r" respectively)}

\item{alpha}{significance level of the test}

\item{j}{index of the coefficient of interest}

\item{MC}{numeric value of Monte Carlo simulations}

\item{method}{defines what method is used to generate data: either division to the "Intervals" or k-nearest neighbors "KNN"}

\item{par}{in Intervals case defines in how many intervals range of Y is divided; in KNN case defines for how many nearest neighbours of each row to include in the generation process itself excluded.}

\item{betabrj}{supposedly true value of coefficient of interest under H_0}
}
\value{
empirical p-value
}
\description{
Calculates empirical p-value of t - or robust t - test
}
\examples{
library(datasets)
 Z <- as.matrix(mtcars)
 X <- Z[,6]
 X <- cbind(rep(1, length(Z[,1])), X)
 Y <- Z[,1]
 #Calling method Intervals
 MCcont(X,Y, 10,40,"t",.05,0,2,10000,"Intervals",1) # for t-test, the range of Y is devided in 1 Interval
 MCcont(X,Y, 10,40,"r",.05,0,2,10000,"Intervals",3) # for robust t-test, the range of Y is devided in 3 Intervals
 ##Calling method KNN
 MCcont(X,Y, 10,40,"t",.05,0,2,10000,"KNN",3) # for t-test, with 3 nearest neighbors initialy assign to the generation process
 MCcont(X,Y, 10,40,"r",.05,0,2,10000,"KNN",6) # for robust t-test, with 6 nearest neighbors initialy assign to the generation process t-test, with 6 nearest neighbors initialy assign to the generation process
}
\author{
RT
}

