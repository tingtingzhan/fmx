
#' @title Moment of Each Component in an \linkS4class{fmx} Object
#' 
#' @description
#' To find moments of each component in an \linkS4class{fmx} object.
#' 
#' @param object an \linkS4class{fmx} object
#' 
#' @details
#' Function [moment_fmx()] calculates the \link[param2moment:moment-class]{moment}s 
#' and distribution characteristics of each mixture component of 
#' an S4 \linkS4class{fmx} object.
#' 
#' @returns 
#' Function [moment_fmx()] returns a \link[param2moment:moment-class]{moment} object.
#' 
#' @examples
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' moment_fmx(d2)
#' 
#' @importFrom param2moment moment_GH moment_norm moment_sn moment_st
#' @export
moment_fmx <- function(object) {
  pars <- object@pars 
  par_nm <- colnames(pars)
  x <- lapply(seq_len(dim(pars)[2L]), FUN = function(i) pars[,i])
  names(x) <- par_nm
  do.call(what = paste0('moment_', object@distname), args = x)
}




#' @title Creates \linkS4class{fmx} Object with Given Component-Wise Moments
#' 
#' @param distname \link[base]{character} scalar
#' 
#' @param w \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... \link[base]{numeric} scalars, 
#' some or all of `mean`, `sd`, `skewness` and `kurtosis`
#' (length will be recycled), see \link[param2moment]{moment2param}
#' 
#' @returns 
#' Function [moment2fmx()] returns a \linkS4class{fmx} object.
#' 
#' @examples
#' m = c(-1.5, 1.5)
#' s = c(.9, 1.1)
#' sk = c(.2, -.3)
#' kt = c(.5, .75)
#' w = c(2, 3)
#' (d1 = moment2fmx(distname='GH', w=w, mean=m, sd=s, skewness=sk, kurtosis=kt))
#' moment_fmx(d1)
#' (d2 = moment2fmx(distname='st', w=w, mean=m, sd=s, skewness=sk, kurtosis=kt))
#' moment_fmx(d2)
#' library(ggplot2)
#' ggplot() + 
#'  geom_function(aes(color = 'GH'), fun = dfmx, args = list(dist=d1), n = 1001) + 
#'  geom_function(aes(color = 'st'), fun = dfmx, args = list(dist=d1), n = 1001) +
#'  xlim(-5, 6)
#' # two curves looks really close, but actually not identical
#' x = rfmx(n = 1e3L, dist = d1)
#' range(dfmx(x, dist = d1) - dfmx(x, dist = d2))
#' 
#' @importFrom param2moment moment2param
#' @export
moment2fmx <- function(distname, w, ...) {
  tmp <- moment2param(distname = distname, ...)
  pars <- do.call(rbind, args = tmp)
  w1 <- cbind(w, pars)[, 1L] # length recycle
  w <- w1/sum(w1)
  new(Class = 'fmx', distname = distname, pars = pars, w = w)
}

