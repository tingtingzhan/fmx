% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moment.R
\name{moment2fmx}
\alias{moment2fmx}
\title{Creates \linkS4class{fmx} Object with Given Component-Wise Moments}
\usage{
moment2fmx(distname, w, ...)
}
\arguments{
\item{distname}{\link[base]{character} scalar}

\item{w}{\link[base]{numeric} \link[base]{vector}}

\item{...}{\link[base]{numeric} scalars,
some or all of \code{mean}, \code{sd}, \code{skewness} and \code{kurtosis}
(length will be recycled), see \link[param2moment]{moment2param}}
}
\value{
Function \code{\link[=moment2fmx]{moment2fmx()}} returns a \linkS4class{fmx} object.
}
\description{
Creates \linkS4class{fmx} Object with Given Component-Wise Moments
}
\keyword{internal}
