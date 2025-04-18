
#' @title Moment of Each Component in an \linkS4class{fmx} Object
#' 
#' @description
#' To find moments of each component in an \linkS4class{fmx} object.
#' 
#' @param object an \linkS4class{fmx} object
#' 
#' @details
#' Function [moment_fmx()] calculates the \linkS4class{moment}s 
#' and distribution characteristics of each mixture component of 
#' an S4 \linkS4class{fmx} object.
#' 
#' @returns 
#' Function [moment_fmx()] returns a \linkS4class{moment} object.
#' 
#' @keywords internal
#' @importFrom param2moment moment_GH moment_norm moment_sn moment_st
#' @importClassesFrom param2moment moment
#' @export
moment_fmx <- function(object) {
  pars <- object@pars 
  par_nm <- colnames(pars)
  x <- lapply(seq_len(dim(pars)[2L]), FUN = \(i) pars[,i])
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
#' @keywords internal
#' @importFrom param2moment moment2param
#' @export
moment2fmx <- function(distname, w, ...) {
  pars <- distname |> 
    moment2param(...) |>
    do.call(what = rbind)
  w1 <- cbind(w, pars)[, 1L] # length recycle
  w <- w1/sum(w1)
  new(Class = 'fmx', distname = distname, pars = pars, w = w)
}

