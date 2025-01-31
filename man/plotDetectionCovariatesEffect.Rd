% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{plotDetectionCovariatesEffect}
\alias{plotDetectionCovariatesEffect}
\title{Estimates of covariates on detection probability}
\usage{
plotDetectionCovariatesEffect(modelResults)
}
\arguments{
\item{modelResults}{Output of the function \code{runModel}}
}
\value{
A plot with of the estimates. The interval represents the 
95\% credible interval.
}
\description{
Plots the estimates of the covariates on the detection probability.
}
\examples{

plotDetectionCovariatesEffect(sampleResults)

}
