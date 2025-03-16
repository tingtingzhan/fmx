
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
#' e.g., `'mean'` and `'sd'` for \link[stats:dnorm]{normal} mixture, 
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
#' @slot Kolmogorov,CramerVonMises,KullbackLeibler (optional) \link[base]{numeric} scalars
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
  Kolmogorov = 'numeric',
  CramerVonMises = 'numeric',
  KullbackLeibler = 'numeric'
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
  
  # all below: optional
  if (length(object@data)) {
    if (anyNA(object@data)) stop('Observations in \'fmx\' must be free of NA_real_')
  }
}) 



