


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




#' @title S3 \link[base]{print} of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param x an \linkS4class{fmx} object
#' 
#' @param ... additional parameters, not currently in use
#' 
#' @returns 
#' Function [print.fmx()] returns the input \linkS4class{fmx} object invisibly.
#' 
#' @keywords internal
#' @export print.fmx
#' @export
print.fmx <- function(x, ...) {
  pars <- x@pars
  pars[] <- sprintf(fmt = '%.2f', pars)
  if (length(id_constr <- fmx_constraint(x))) pars[id_constr] <- '.'
  K <- dim(pars)[1L]
  rownames(pars) <- paste0(seq_len(K), '-comp.')
  obj <- if (K == 1L) pars else cbind(pars, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  names(dimnames(obj)) <- c(x@distname, 'Parameters')
  
  if (length(x@data)) {
    # theoretically does not make sense to talk about confidence interval, without a sample
    ci <- confint.fmx(x, level = .95, internal = FALSE)
    id_constr <- fmx_constraint(x)
    if (length(ci) && !anyNA(ci)) {
      ci0 <- sprintf(fmt = '(%.2f~%.2f)', ci[,1L], ci[,2L])
      if (length(id_constr)) {
        obj[id_constr] <- '.'
        obj[-id_constr] <- paste(obj[-id_constr], ci0)
      } else obj[] <- paste(obj, ci0)
    }
  }
  
  print.default(obj, quote = FALSE)
  return(invisible(x))
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
#' the row indices of *components* to be chosen, see \link[base]{[}
#' 
#' @details 
#' 
#' Using definitions as S3 method dispatch \code{`[.fmx`} won't work 
#' for \linkS4class{fmx} objects.
#' 
#' @returns 
#' 
#' An \linkS4class{fmx} object consisting of a subset of components.
#' information about the observations (e.g. slots `@@data` and `@@data.name`),
# as well as other estimation related slots (e.g., `@@init`) 
#' will be lost.
#' 
#' @examples 
#' 
#' (d = fmx('norm', mean = c(1, 4, 7), w = c(1, 1, 1)))
#' d[1:2]
#' 
#' @keywords internal
#' @export
setMethod(`[`, signature(x = 'fmx', i = 'ANY'), definition = function(x, i) {
  if (length(x@data)) message('Subset the estimates and drop `@data` etc.')
  pars <- x@pars[i, , drop = FALSE]
  w <- x@w[i]
  w <- unname(w / sum(w)) # adjust mixing proportions
  o <- order(pars[, 1L])
  new(Class = 'fmx', pars = pars[o, , drop = FALSE], w = w[o], distname = x@distname)
})
















#' @title Confidence Interval of \linkS4class{fmx} Object
#' 
#' @description ...
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param level confidence level, default \eqn{95\%}.
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' Function [confint.fmx()] returns the Wald-type confidence intervals based on the user-friendly parameters (`parm = 'user'`),
#'  or the internal/unconstrained parameters (`parm = 'internal'`).
#' When the distribution has constraints on one or more parameters, 
#' function [confint.fmx()] does not return the confident intervals of for the constrained parameters.
#'  
#' @returns 
#' Function [confint.fmx()] returns a \link[base]{matrix}
#' 
#' @keywords internal
#' @importFrom stats confint qnorm
#' @export confint.fmx
#' @export
confint.fmx <- function(object, ..., level = .95) {
  # essentially ?stats::confint.default
  cf <- coef.fmx(object, ...)
  if (!length(vv <- vcov.fmx(object, ...))) return(invisible())
  ses <- sqrt(diag(vv))
  p1 <- (1 - level) / 2
  p <- c(p1, 1 - p1)
  ret <- cf + ses %*% t.default(qnorm(p))
  dimnames(ret) <- list(names(cf), sprintf('%.1f%%', 1e2*p))
  return(ret)
}



#' @title Variance-Covariance of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param internal \link[base]{logical} scalar, either for the user-friendly parameters (`FALSE`, default)
#' (e.g., `mean,sd` for normal mixture, and `A,B,g,h` for Tukey \eqn{g}-and-\eqn{h} mixture), or
#' for the internal/unconstrained parameters (`TRUE`).
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' Function [vcov.fmx()] returns 
#' the approximate asymptotic variance-covariance \link[base]{matrix} of the user-friendly parameters via delta-method (`parm = 'user'`), 
#' or the asymptotic variance-covariance matrix of the internal/unconstrained parameters (`parm = 'internal'`). 
#' When the distribution has constraints on one or more parameters, 
#' function [vcov.fmx()] does not return the variance/covariance involving the constrained parameters.
#' 
#' @returns 
#' 
#' Function [vcov.fmx()] returns a \link[base]{matrix}.
#' 
#' @keywords internal
#' @importFrom stats vcov
#' @export vcov.fmx
#' @export
vcov.fmx <- function(object, internal = FALSE, ...) {
  
  if (!length(object@data)) return(invisible())
  
  if (!internal && length(vv <- object@vcov)) return(vv) # 'fitdist' objects
  
  int_vv <- object@vcov_internal
  if (internal) return(int_vv)
  
  if (!length(int_vv)) return(int_vv) # wont be able to computer user-vcov if internal-vcov is wrong
  
  distname <- object@distname
  pars <- object@pars
  K <- dim(pars)[1L]
  int_nm <- dimnames(int_vv)[[1L]]
  
  int_p <- fmx2dbl(object) # internal parameters
  anm <- distArgs(distname)
  n_anm <- length(anm)
  user_nm <- c(t.default(outer(c(anm, if (K > 1L) 'w'), 1:K, FUN = paste0)))
  jacob <- array(0, dim = c(length(int_p), length(user_nm)), dimnames = list(names(int_p), user_nm))
  # location parameters A_1 -> A_1
  jacob[1L, 1:K] <- 1
  if (K > 1L) {
    # location parameters A_2, .., A_k -> d_2, .., d_k
    for (k in 2:K) jacob[k, k:K] <- exp(int_p[k]) 
    # mixture parameters w_1, .., w_k -> pi_2, .., pi_k
    id_pi <- (n_anm*K+1L):((n_anm+1L)*K-1L)
    e_pi <- exp(int_p[id_pi])
    sum_pi <- sum(1 + e_pi)^2 # 1 = exp(pi_1) = exp(0)
    jacob[id_pi, n_anm*K+1L] <- - e_pi / sum_pi
    jacob[id_pi, id_pi+1L] <- - tcrossprod(e_pi) / sum_pi
  }
  
  switch(distname, norm = {
    id_exp <- (K+1L):(2*K) # 'sd'
    id_identity <- id_constr <- NULL
  }, GH = {
    id_exp <- c((K+1L):(2*K), (3*K+1L):(4*K)) # 'B' and 'h'
    id_identity <- (2*K+1L):(3*K) # 'g'
    id_constr <- fmx_constraint(object)
  })
  
  if (length(id_exp)) jacob[cbind(id_exp, id_exp)] <- exp(int_p[id_exp])
  if (length(id_identity)) jacob[cbind(id_identity, id_identity)] <- 1
  jacob_free <- if (length(id_constr)) jacob[int_nm, -id_constr] else jacob
  
  return(t.default(jacob_free) %*% int_vv %*% jacob_free)
  
}





# ?stats::residuals 
# @export
# residuals.fmx <- function(object, ...) stop('useful?')





#' @title Parameter Estimates of \linkS4class{fmx} object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param internal \link[base]{logical} scalar, either for the user-friendly parameters (`FALSE`, default)
#' (e.g., `mean,sd` for normal mixture, and `A,B,g,h` for Tukey \eqn{g}-and-\eqn{h} mixture), or
#' for the internal/unconstrained parameters (`TRUE`).
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' Function [coef.fmx()] returns the estimates of the user-friendly parameters (`parm = 'user'`), 
#' or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' When the distribution has constraints on one or more parameters, 
#' function [coef.fmx()] does not return the estimates (which is constant \code{0}) of the constrained parameters.
#' 
#' @returns 
#' 
#' Function [coef.fmx()] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @keywords internal
#' @importFrom stats coef
#' @export coef.fmx
#' @export
coef.fmx <- function(object, internal = FALSE, ...) {
  anm <- distArgs(object@distname)
  K <- dim(object@pars)[1L]
  cf0 <- if (internal) fmx2dbl(object) else if (K == 1L) {
    setNames(c(object@pars), nm = c(t.default(outer(anm, 1:K, FUN = paste0))))
  } else {
    setNames(c(object@pars, object@w), nm = c(t.default(outer(c(anm, 'w'), 1:K, FUN = paste0))))
  }
  if (!length(id_constr <- fmx_constraint(object))) return(cf0)
  return(cf0[-id_constr])
}






#' @title Log-Likelihood of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param data \link[base]{double} \link[base]{vector}, actual observations
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' Function [logLik.fmx()] returns a \link[stats]{logLik} object indicating the log-likelihood.
#' An additional attribute `attr(,'logl')` indicates the point-wise log-likelihood, 
#' to be use in Vuong's closeness test.
#' 
#' @returns 
#' 
#' Function [logLik.fmx()] returns a \link[stats]{logLik} object with 
#' an additional attribute `attr(,'logl')`.
#' 
#' @keywords internal
#' @importFrom stats logLik
#' @export logLik.fmx
#' @export
logLik.fmx <- function(object, data = object@data, ...) {
  
  # for developer to batch-calculate AIC/BIC quickly
  #if (length(objF <- attr(object, which = 'objF', exact = TRUE))) {
  #  if (inherits(objF[[1L]], what = 'logLik')) return(objF[[1L]])
  #}
  # ?step_fmx no longer uses logLik
  
  if (!length(data)) return(invisible())
  
  logd <- dfmx(x = data, dist = object, log = TRUE, ...)
  if (!all(is.finite(logd))) {
    #print(logd)
    #object <<- object
    #stop('malformed fit (?.dGH has been well debug-ged)')
    # very likely to be `B = 0`
    # do not stop.  settle with -Inf log-likelihood
  }
  ret <- sum(logd)
  attr(ret, which = 'logl') <- logd # additional attributes; needed in Vuong's test
  attr(ret, which = 'nobs') <- length(data)
  attr(ret, which = 'df') <- npar.fmx(object)
  class(ret) <- 'logLik'
  return(ret)
}






