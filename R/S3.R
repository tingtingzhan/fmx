



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
  
  if (length(x@data)) {
    
    cat('\n')
    
    x@data |> 
      length() |>
      sprintf(fmt = 'Number of Observations: %d\n') |> 
      cat()
    
    cat('\n')
    
    x@logLik |>
      print()
    
    cat('\n')
    
    x@dist.ks |>
      sprintf(fmt = 'Kolmogorov-Smirnov Distance: %.5f\n') |> 
      cat()
    
    x@dist.cvm |>
      sprintf(fmt = 'Cramer-Von Mises Distance: %.5f\n') |> 
      cat()
    
    x@dist.kl |>
      sprintf(fmt = 'Kullback-Leibler divergence: %.5f\n') |> 
      cat()

  }
  
  return(invisible(x))
}














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
#' @description 
#' Log-likelihood of an \linkS4class{fmx} object.
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @returns 
#' Function [logLik.fmx()] returns a \link[stats]{logLik} object.
#' 
#' @keywords internal
#' @importFrom stats logLik
#' @export logLik.fmx
#' @export
logLik.fmx <- function(object, ...) object@logLik



#' @title Number of Observations in \linkS4class{fmx} Object
#' 
#' @description 
#' Number of observations in an \linkS4class{fmx} object.
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @details 
#' Function [nobs.fmx()] finds the sample size of `@data` slot of
#' an \linkS4class{fmx} object.
#' 
#' @returns 
#' Function [nobs.fmx()] returns an \link[base]{integer} scalar.
#' 
#' @keywords internal
#' @importFrom stats nobs
#' @export nobs.fmx
#' @export
nobs.fmx <- function(object, ...) length(object@data)





