



#' @title Distribution Type
#' 
#' @description ..
#' 
#' @param type \link[base]{character} scalar
#' 
#' @returns 
#' Function [distType()] returns a \link[base]{character} \link[base]{vector}. 
#' 
#' @keywords internal
#' @export
distType <- function(type = c('discrete', 'nonNegContinuous', 'continuous')) {
  switch(match.arg(type), discrete = c(
    # integer distribution on [0, Inf)
    'binom', 'geom', 'hyper', 'nbinom', 'pois', 'wilcox', 'signrank', # \pkg{stats}
    'genpois1', # ?VGAM::dgenpois1
    NULL
  ), nonNegContinuous = c(
    # *p*ositive continuous distribution on (0, Inf) or [0, Inf)
    'beta', 'chisq', 'exp', 'f', 'gamma', 'lnorm', 'weibull', # \pkg{stats}
    'invgauss', # ?statmod::dinvgauss
    'pareto', # ?sads::dpareto
    'invweibull', # ?actuar::dinvweibull
    'llogis', # function name clash! ?STAR::dllogis (error when x = 0) and ?actuar::dllogis
    'normp', # ?normalp::dnormp not supported, as ?normalp::qnormp has `pr` as first arg
    NULL
  ), continuous = c(
    # continuous distribution on (-Inf, Inf)
    'cauchy', 'norm', 't', 'unif', 'logis', # \pkg{stats}
    'RevGumbel', 'RevWeibull', # ?DescTools::dRevGumbel, ?DescTools::dRevWeibull
    'gumbel', 'lgamma', # ?ordinal::dgumbel, ?ordinal::dlgamma.  Name clash ?VGAM::dgumbel, ?VGAM::dlgamma
    # 'lnorm3', # needs rethink
    'sn', 'st', # skew-normal, skew-t from \CRANpkg{sn}
    'GH', # my \link[TukeyGH77]{dGH}
    # MCMCpack::dinvgamma # do not have p* and q* function
    NULL
  ))
}







#' @title Name(s) of Formal Argument(s) of Distribution
#' 
#' @description 
#' To obtain the name(s) of distribution parameter(s).
#' 
#' @param distname \link[base]{character} scalar, name of distribution
#' 
#' @returns 
#' Function [distArgs()] returns a \link[base]{character} \link[base]{vector}.
#' 
#' @seealso 
#' \link[methods]{formalArgs}
#' 
#' @keywords internal
#' @export
distArgs <- function(distname) {
  switch(distname, 
         beta = c('shape1', 'shape2', 'ncp'), # ?stats::dbeta
         binom = c('size', 'prob'), # ?stats::dbinom
         cauchy = c('location', 'scale'), # ?stats:::dcauchy
         chisq = c('df', 'ncp'), # ?stats:::dchisq
         gamma = c('shape', 'scale'), # ?stats::dgamma
         # read \link[mixtools]{gammamixEM} carefully: 'beta' is actually `scale`
         nbinom = c('size', 'prob'), # ?stats::dnbinom
         norm = c('mean', 'sd'), # ?stats::dnorm
         GH = c('A', 'B', 'g', 'h'), # TukeyGH77::dGH
         sn = c('xi', 'omega', 'alpha'), # ?sn::dsn # I do not understand 'tau'
         st = c('xi', 'omega', 'alpha', 'nu'), # ?sn::dst
         stop(sQuote(distname), ' not supported yet'))
}


#' @title Distribution Parameters that needs to have a log-transformation
#' 
#' @description ..
#' 
#' @param distname \link[base]{character} scalar, name of distribution
#' 
#' @returns 
#' Function [dist_logtrans()] returns an \link[base]{integer} scalar
#' 
#' @keywords internal
#' @export
dist_logtrans <- function(distname) {
  switch(distname, norm = {
    2L # `sdlog`
  }, GH = {
    c(2L, 4L) # `Blog` & `hlog`
  }, stop('distribution', sQuote(distname), 'not supported yet'))
}



