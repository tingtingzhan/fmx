
#' @title \linkS4class{fmx} Class: Finite Mixture Parametrization
#' 
#' @description 
#' An S4 object to specify 
#' the parameters and type of distribution 
#' of a one-dimensional finite mixture distribution.
#'
#' @slot distname \link[base]{character} scalar, 
#' name of parametric distribution of the mixture components.
#' Currently, normal (`'norm'`) and Tukey \eqn{g}-&-\eqn{h} (`'GH'`) distributions are supported.
#' 
#' @slot pars \link[base]{double} \link[base]{matrix}, 
#' all distribution parameters in the mixture. 
#' Each row corresponds to one component. Each column includes the same parameters of all components.
#' The order of rows corresponds to the (non-strictly) increasing order of the component location parameters.
#' The columns match the formal arguments of the corresponding distribution, 
#' e.g., `'mean'` and `'sd'` for \link[stats]{Normal} mixture, 
#' or `'A'`, `'B'`, `'g'` and `'h'` for Tukey \eqn{g}-&-\eqn{h} mixture.
#' 
#' @slot w \link[base]{numeric} \link[base]{vector} of mixing proportions that must sum to 1
#' 
#' @slot data (optional) \link[base]{numeric} \link[base]{vector}, the one-dimensional observations
#' 
#' @slot data.name (optional) \link[base]{character} scalar, a human-friendly name of the observations
#'
#' @slot vcov_internal (optional) variance-covariance \link[base]{matrix} of the internal (i.e., unconstrained) estimates
#' 
#' @slot vcov (optional) variance-covariance \link[base]{matrix} of the mixture distribution (i.e., constrained) estimates
#' 
#' @slot dist.ks (optional) \link[base]{numeric} scalar, Kolmogorov-Smirnov distance, via \link[stats]{ks.test}
#' 
#' @slot dist.cvm (optional) \link[base]{numeric} scalars, Cramer von Mises distance, via \link[goftest]{cvm.test}
#' 
#' @slot dist.kl (optional) \link[base]{numeric} scalars, Kullback-Leibler distance
#' 
#' @export
setClass(Class = 'fmx', slots = c(
  distname = 'character',
  pars = 'matrix',
  w = 'numeric',
  ### all below: optional
  data = 'numeric', 
  data.name = 'character',
  vcov_internal = 'matrix',
  vcov = 'matrix',
  ### all below: diagnostics
  dist.ks = 'numeric',
  dist.cvm = 'numeric',
  dist.kl = 'numeric'
), prototype = prototype(
  w = 1 # for 1-component
), validity = function(object) {
  pars <- object@pars
  if (anyNA(pars)) stop('do not allow NA in `fmx` distribution parameter')
  if (is.unsorted(pars[,1L], strictly = FALSE)) { 
    # ?base::is.unsorted.  Note here I use `strictly = FALSE`
    # since log(d) may have large negative value, such that exp(log(d)) may be numerically 0.
    stop('location parameter must be sorted: ', sQuote(pars[,1L])) 
  }
  K <- dim(pars)[1L]
  w <- object@w
  if (K == 1L) {
    if (!isTRUE(all.equal.numeric(w, 1))) stop('slot `w` should be `1` for 1-component') # may not be ?base::identical
  } else {
    if (length(w) != K) stop('slot `w` should be length-K')
    if (!isTRUE(all.equal.numeric(sum(w), 1))) stop('slot `w` must sum up to be `1`') # may not be ?base::identical
  }
}) 


#' @importFrom stats ks.test
#' @importFrom goftest cvm.test
# @seealso `dgof::cvm.test`
#' @importFrom LaplacesDemon KLD
setMethod(f = initialize, signature = 'fmx', definition = function(.Object, ...) {
  
  x <- callNextMethod(.Object, ...)
  
  data <- x@data
  if (!length(data)) return(x)
  
  if (anyNA(data)) stop('Observations in \'fmx\' must be free of NA_real_')

  # overwrite existing
  x@dist.ks <- ks.test(x = data, y = pfmx, dist = x)$statistic |> unname() |> suppressWarnings()
  
  x@dist.cvm <- cvm.test(
    x = data, 
    #null = pfmx, dist = x, # error, do not know why..
    null = function(q, dist = x) {
      pfmx(q, dist = dist) |>
        pmax.int(.Machine$double.eps) |>
        pmin.int(1 - .Machine$double.eps)
    }, 
    nullname = ''
  )$statistic |>
    unname() |> 
    suppressWarnings()
  
  kl.px <- dfmx(data, dist = x, log = FALSE)
  kl.py <- if (x@distname %in% c(distType('continuous'), distType('nonNegContinuous'))) {
    approxdens(data)(data)
  } else (tabulate(data, nbins = max(data)) / length(data))[data]
  x@dist.kl <- KLD(px = kl.px, py = kl.py)$sum.KLD.py.px
  
  return(x)
  
})






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







#' @title Subset of Components in \linkS4class{fmx} Object
#' 
#' @description 
#' 
#' Taking subset of components in \linkS4class{fmx} object
#' 
#' @param x \linkS4class{fmx} object
#' 
#' @param i \link[base]{integer} or \link[base]{logical} \link[base]{vector}, 
#' the *row* indices of *components* to be chosen, see \link[base]{[}
#' 
#' @returns 
#' 
#' An \linkS4class{fmx} object consisting of a subset of components.
#' Information about the observations (e.g. slots `@@data` and `@@data.name`),
# as well as other estimation related slots (e.g., `@@init`) 
#' will be lost.
#' 
#' @keywords internal
#' @export
`[.fmx` <- function(x, i) {
  if (length(x@data)) message('Subset the estimates and drop `@data` etc.')
  pars <- x@pars[i, , drop = FALSE]
  w <- x@w[i]
  w <- unname(w / sum(w)) # adjust mixing proportions
  o <- order(pars[, 1L])
  new(Class = 'fmx', pars = pars[o, , drop = FALSE], w = w[o], distname = x@distname)
}


#' @title Show \linkS4class{fmx} Object
#' 
#' @description
#' Print the parameters of an \linkS4class{fmx} object and plot its density curves.
#' 
#' @param object an \linkS4class{fmx} object
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{fmx} object 
#' does not have a returned value.
#' 
#' @keywords internal
#' @export
setMethod(f = show, signature = signature(object = 'fmx'), definition = function(object) {
  print.fmx(object)
})








