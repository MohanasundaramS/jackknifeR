% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jackknife.lm.R
\name{jackknife.lm}
\alias{jackknife.lm}
\title{Delete-d Jackknife Estimate for Linear Regression}
\usage{
jackknife.lm(formula, d = 1, data, conf = 0.95, numCores = detectCores())
}
\arguments{
\item{formula}{Simple or multiple  linear regression formula with dependent and independent variables}

\item{d}{Number of observations to be deleted from data to make jackknife samples. The default is 1 (for delete-1 jackknife).}

\item{data}{Data frame with dependent and independent independent variables specified in the formula}

\item{conf}{Confidence level, a positive number < 1. The default is 0.95.}

\item{numCores}{Number of processors to be used}
}
\value{
A list containing a summary data frame of jackknife estimates
with bias, standard error. t-statistics, and confidence intervals,
linear regression model of original data and a data frame with
coefficient estimates of jackknife samples.
}
\description{
This function creates jackknife samples from the data by
sequentially removing \emph{d} observations from the data,
fits models linear regression model using the jackknife samples
as specified in the formula and estimates the jackknife coefficients
bias standard error, standard error and confidence intervals.
}
\examples{
## library(jackknifeR)
jk <- jackknife.lm(mpg ~ wt + hp, d = 2, data = mtcars, numCores = 2)
summary(jk)
}
\references{
Quenouille, M. H. (1956). Notes on Bias in Estimation.
\emph{Biometrika}, \emph{43}(3/4), 353-360.
\doi{10.2307/2332914}

Tukey, J. W. (1958). Bias and Confidence in Not-quite Large Samples.
\emph{Annals of Mathematical Statistics}, \emph{29}(2), 614-623.
\doi{10.1214/aoms/1177706647}

Shi, X. (1988). A note on the delete-d jackknife variance estimators.
\emph{Statistics & Probability Letters}, \emph{6}(5), 341-347.
\doi{10.1016/0167-7152(88)90011-9}
}
\seealso{
\code{\link[=lm]{lm()}} which is used for linear regression.
}
