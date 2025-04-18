% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3.R
\name{vcov.fmx}
\alias{vcov.fmx}
\title{Variance-Covariance of \linkS4class{fmx} Object}
\usage{
\method{vcov}{fmx}(object, internal = FALSE, ...)
}
\arguments{
\item{object}{\linkS4class{fmx} object}

\item{internal}{\link[base]{logical} scalar, either for the user-friendly parameters (\code{FALSE}, default)
(e.g., \verb{mean,sd} for normal mixture, and \verb{A,B,g,h} for Tukey \eqn{g}-and-\eqn{h} mixture), or
for the internal/unconstrained parameters (\code{TRUE}).}

\item{...}{place holder for S3 naming convention}
}
\value{
Function \code{\link[=vcov.fmx]{vcov.fmx()}} returns a \link[base]{matrix}.
}
\description{
..
}
\details{
Function \code{\link[=vcov.fmx]{vcov.fmx()}} returns
the approximate asymptotic variance-covariance \link[base]{matrix} of the user-friendly parameters via delta-method (\code{parm = 'user'}),
or the asymptotic variance-covariance matrix of the internal/unconstrained parameters (\code{parm = 'internal'}).
When the distribution has constraints on one or more parameters,
function \code{\link[=vcov.fmx]{vcov.fmx()}} does not return the variance/covariance involving the constrained parameters.
}
\keyword{internal}
