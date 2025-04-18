

#' @title Turn Various Objects to \linkS4class{fmx} Class
#' 
#' @description 
#' 
#' Turn various objects created in other R packages 
#' to \linkS4class{fmx} class.
#' 
#' @param x an R object
#' 
#' @param ... additional parameters, see **Arguments** in individual S3 dispatches
#' 
#' @details 
#' Various mixture distribution estimates obtained from other R packages
#' are converted to \linkS4class{fmx} class, 
#' so that we could take advantage of all methods defined for \linkS4class{fmx} objects.
#' 
#' @returns
#' S3 generic function [as.fmx()] returns an \linkS4class{fmx} object.
#' 
#' @export
as.fmx <- function(x, ...) UseMethod('as.fmx')


#' @export
as.fmx.fmx <- function(x, ...) x




#' @title Convert `mixEM` Objects to \linkS4class{fmx} Class
#' 
#' @description
#' To convert `mixEM` objects (from package \CRANpkg{mixtools}) 
#' to \linkS4class{fmx} class.
#' 
#' Currently only the returned value of 
#' \link[mixtools]{normalmixEM} and \link[mixtools]{gammamixEM} are supported
#' 
#' @param x `mixEM` object
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @note 
#' \link[mixtools]{plot.mixEM} not plot \link[mixtools]{gammamixEM} returns, as of 2022-09-19.
#' 
#' @returns 
#' Function [as.fmx.mixEM()] returns an \linkS4class{fmx} object.
#' 
#' @keywords internal
#' @method as.fmx mixEM
#' @export as.fmx.mixEM
#' @export
as.fmx.mixEM <- function(x, data = x[['x']], ...) {
  if (!length(data)) stop('wont happen')
  x <- sort.mixEM(x, decreasing = FALSE)
  
  switch(x[['ft']], normalmixEM = {
    pars <- cbind(mean = x[['mu']], sd = x[['sigma']])
    distname <- 'norm'
  }, gammamixEM = {
    pars <- t.default(x[['gamma.pars']])
    colnames(pars) <- c('shape', 'scale') # names of parameters of ?stats::dgamma
    # read \link[mixtools]{gammamixEM} carefully: 'beta' is actually `scale`
    distname <- 'gamma'
  }, stop(x[['ft']], ' not supported yet'))
  new(Class = 'fmx', 
      pars = pars, w = x[['lambda']], distname = distname,
      data = data)
}









#' @title Convert `Skew.normal` Object to \linkS4class{fmx}
#' 
#' @description 
#' To convert `Skew.normal` object (from package \CRANpkg{mixsmsn}) 
#' to \linkS4class{fmx} class.
#' 
#' @param x `'Skew.normal'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Skew.normal'`.
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Skew.normal()] returns an \linkS4class{fmx} object.
#' 
#' @keywords internal
#' @method as.fmx Skew.normal
#' @export as.fmx.Skew.normal
#' @export
as.fmx.Skew.normal <- function(x, data, ...) {
  x <- sort.Skew.normal(x, decreasing = FALSE)
  
  ret <- new(Class = 'fmx', pars = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']]
  ), 
  w = x[['pii']],
  distname = 'sn')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
}





#' @title Convert `Normal` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description
#' To convert `Normal` object (from package \CRANpkg{mixsmsn}) 
#' to \linkS4class{fmx} class.
#' 
#' @param x `'Normal'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Normal'`
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Normal()] returns an \linkS4class{fmx} object.
#' 
#' @keywords internal
#' @method as.fmx Normal
#' @export as.fmx.Normal
#' @export
as.fmx.Normal <- function(x, data, ...) {
  x <- sort.Normal(x, decreasing = FALSE)
  
  ret <- new(Class = 'fmx', pars = cbind(
    mean = x[['mu']],
    sd = sqrt(x[['sigma2']])
  ), 
  w = x[['pii']],
  distname = 'norm')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
  
}






#' @title Convert `Skew.t` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description
#' To convert `Skew.t` object (from package \CRANpkg{mixsmsn}) 
#' to \linkS4class{fmx} class.
#' 
#' @param x `'Skew.t'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Skew.t'` 
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Skew.t()] returns an \linkS4class{fmx} object.
#' 
#' @keywords internal
#' @method as.fmx Skew.t
#' @export as.fmx.Skew.t
#' @export
as.fmx.Skew.t <- function(x, data, ...) {
  x <- sort.Skew.t(x, decreasing = FALSE)
  
  K <- length(x[['mu']]) # number of components
  if (length(x[['nu']]) != 1L) stop('\\pkg{mixsmsn} update to enable multiple `nu`? Modify ?npar.fmx') 
  
  ret <- new(Class = 'fmx', pars = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']],
    nu = x[['nu']]
  ), 
  w = x[['pii']],
  distname = 'st')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
}






#' @title Convert `t` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description
#' To convert `t` object (from package \CRANpkg{mixsmsn}) 
#' to \linkS4class{fmx} class.
#' 
#' @param x `'t'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 't'`
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.t()] has not been completed yet
#' 
#' @keywords internal
#' @method as.fmx t
#' @export
as.fmx.t <- function(x, data, ...) {
  if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
  x <- sort.t(x, decreasing = FALSE)
  stop('need to programe scale_and_shift_t')
  #new(Class = 'fmx', pars = cbind(
  #  mean = x[['mu']],
  #  sd = sqrt(x[['sigma2']])
  #), 
  #w = x[['pii']],
  #distname = 'norm', 
  #data = data)
}
