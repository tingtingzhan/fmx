

#' @title Create \linkS4class{fmx} Object for Finite Mixture Distribution
#' 
#' @description 
#' To create \linkS4class{fmx} object for finite mixture distribution.
#' 
#' @param distname \link[base]{character} scalar
#' 
#' @param w (optional) \link[base]{numeric} \link[base]{vector}.  
#' Does not need to sum up to 1; `w/sum(w)` will be used internally.
#' 
#' @param ... mixture distribution parameters.
#' See function \link[TukeyGH77]{dGH} for the names and default values of Tukey \eqn{g}-&-\eqn{h} distribution parameters, 
#' or \link[stats]{dnorm} for the names and default values of normal distribution parameters.
#' 
#' @returns 
#' Function [fmx()] returns an \linkS4class{fmx} object.
#' 
#' @importFrom TukeyGH77 dGH
#' @export
fmx <- function(distname, w = 1, ...) {
  if (!is.character(distname) || (length(distname) != 1L) || anyNA(distname) || !nzchar(distname)) stop('distname must be len-1 char')
  anm <- distArgs(distname)
  ddist <- paste0('d', distname)
  farg <- formals(ddist)[anm] # if `farg` has empty element, `do.call(cbind, farg)` will err
  
  K <- length(w)
  if (!is.numeric(w) || !K || anyNA(w) || any(w <= 0)) stop('illegal mixing proportions `w`')
  
  arg <- list(...)[anm]
  names(arg) <- anm
  if (!any(la <- lengths(arg, use.names = FALSE))) {
    message('Using default arguments of ', sQuote(ddist))
    if (K != 1L) stop('must specify at least one non-equal parameter for K>2 mixture')
    return(new(Class = 'fmx', pars = do.call(cbind, args = farg), w = 1, distname = distname))
  }
  
  for (id in which(la == 0L)) arg[id] <- farg[id]
  if (!all(vapply(arg, FUN = is.numeric, FUN.VALUE = NA))) stop('distribution parameters must be numeric')
  if (anyNA(arg, recursive = TRUE)) stop('do not allow NA in `arg`')
  pars <- do.call(cbind, args = arg) # vector recycling
  if (dim(pars)[1L] != K) stop('parameter (formal args and `w`) lengths not match')
  if (is.unsorted(loc <- pars[,1L], strictly = FALSE)) {
    message('Re-ordered by location parameter')
    pars <- pars[order(loc, decreasing = FALSE), , drop = FALSE]
  }
  new(Class = 'fmx', pars = pars, w = unname(w/sum(w)), distname = distname)
}
