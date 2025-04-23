

#' @title Density, Distribution and Quantile of Finite Mixture Distribution
#' 
#' @description 
#' 
#' Density function, distribution function, quantile function and random generation for a finite mixture distribution 
#' with normal or Tukey \eqn{g}-&-\eqn{h} components.
#' 
#' @param x,q \link[base]{numeric} \link[base]{vector}, quantiles, `NA_real_` value(s) allowed.
#' 
#' @param p \link[base]{numeric} \link[base]{vector}, probabilities.
#' 
#' @param n \link[base]{integer} scalar, number of observations.
#' 
#' @param dist \linkS4class{fmx} object, a finite mixture distribution
#' 
#' @param log,log.p \link[base]{logical} scalar. 
#' If `TRUE`, probabilities are given as \eqn{\log(p)}.
#' 
#' @param lower.tail \link[base]{logical} scalar. 
#' If `TRUE` (default), probabilities are \eqn{Pr(X\le x)}, otherwise, \eqn{Pr(X>x)}.
#' 
#' @param interval \link[base]{length}-2 \link[base]{numeric} \link[base]{vector}, interval for root finding, see \link[rstpm2]{vuniroot}
#' 
#' @param distname,K,pars,w auxiliary parameters, whose default values are determined by argument `dist`.
#' The user-specified \link[base]{vector} of `w` does not need to sum up to 1; `w/sum(w)` will be used internally.
#' 
#' @param ... additional parameters
#' 
#' @details 
#' 
#' A computational challenge in function [dfmx()] is when mixture density is very close to 0,
#' which happens when the per-component log densities are negative with big absolute values.  
#' In such case, we cannot compute the log densities (i.e., `-Inf`).
#' 
#' Function [qfmx()] gives the quantile function, by numerically solving [pfmx].
#' One major challenge when dealing with the finite mixture of Tukey \eqn{g}-&-\eqn{h} family distribution
#' is that Brent–Dekker's method needs to be performed in both \link[TukeyGH77]{pGH} and [qfmx] functions, 
#' i.e. *two layers* of root-finding algorithm.
#' 
#' 
#' @returns 
#' 
#' Function [dfmx()] returns a \link[base]{numeric} \link[base]{vector} of probability density values of an \linkS4class{fmx} object at specified quantiles `x`.
#' 
#' Function [pfmx()] returns a \link[base]{numeric} \link[base]{vector} of cumulative probability values of an \linkS4class{fmx} object at specified quantiles `q`.
#' 
#' Function [qfmx()] returns an unnamed \link[base]{numeric} \link[base]{vector} of quantiles of an \linkS4class{fmx} object, based on specified cumulative probabilities `p`.
#' 
#' Function [rfmx()] generates random deviates of an \linkS4class{fmx} object.
#' 
#' @note
#' Function \link[stats]{qnorm} returns an unnamed \link[base]{vector} of quantiles, 
#' although \link[stats]{quantile} returns a named \link[base]{vector} of quantiles.
#' 
#' @keywords internal
#' @name dfmx
#' @import stats
#' @importFrom sn dsn psn qsn rsn dst pst qst rst
#' @importFrom TukeyGH77 dGH
#' @export
dfmx <- function(x, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, ..., log = FALSE) {
  
  if (K == 1L) { # no mixture required!!
    switch(distname, 
           gamma = return(dgamma(x = x, shape = pars[,1L], scale = pars[,2L], log = log)),
           nbinom = return(dnbinom(x = x, size = pars[,1L], prob = pars[,2L], log = log)),
           norm = return(dnorm(x = x, mean = pars[,1L], sd = pars[,2L], log = log)), 
           GH = return(dGH(x = x, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], log = log, ...)),
           sn = return(dsn(x = x, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], log = log)),
           st = return(dst(x = x, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L], log = log)),
           stop('I do not have `d', distname, '` function'))
  }
  
  xm <- tcrossprod(rep(1, times = K), x)
  
  lds <- switch(distname, # `lds` is per-component log-densities
                gamma = dgamma(x = xm, shape = pars[,1L], scale = pars[,2L], log = TRUE),
                nbinom = dnbinom(x = xm, size = pars[,1L], prob = pars[,2L], log = TRUE),
                norm = dnorm(x = xm, mean = pars[,1L], sd = pars[,2L], log = TRUE), 
                GH = {
                  # ?TukeyGH77::dGH only takes scalar A/B/g/h
                  K |>
                    seq_len() |>
                    lapply(FUN = \(i) dGH(x = x, A = pars[i,1L], B = pars[i,2L], g = pars[i,3L], h = pars[i,4L], log = TRUE)) |>
                    do.call(what = rbind)
                },
                sn = {
                  # ?sn::dsn does not respect `attr(x, 'dim')`
                  # to make things worse, ?sn::dsn does not even handle \link[base]{matrix} `x` correctly!!
                  # array(dsn(x = xm, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], log = TRUE), dim = dim(xm), dimnames = dimnames(xm))
                  # this is wrong!!
                  # have to go the stupid way!!
                  do.call(rbind, args = lapply(seq_len(K), FUN = \(i) dsn(x = x, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], log = TRUE)))
                }, 
                st = {
                  # ?sn::dst gives error on vector `nu`
                  do.call(rbind, args = lapply(seq_len(K), FUN = \(i) dst(x = x, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], nu = pars[i,4L], log = TRUE)))
                },
                stop('I do not have `d', distname, '` function'))
  if (any(is.infinite(lds))) {
    #tmp <<- dist; x <<- x;
    #stop('per-component log-density should not be -Inf (unless `sd` or `B` is 0)') # is.infinite(NA) is `FALSE`
    # `sd = 0` or `B = 0` may happen 
  }
  
  d <- c(crossprod(w, exp(lds)))
  
  if (!log) {
    attr(d, which = 'posterior') <- w * exp(lds)
    # range(colSums(attr(d, which = 'posterior')) - d) # super small
    return(d)
  }
  #any(is.infinite(log(d))) # could happen
  
  nx <- length(x)
  wlds <- log(w) + lds # weighted-lds
  id_max <- max.col(t.default(wlds))  # dont want to import ?matrixStats::colMaxs
  id_max_seq <- cbind(id_max, seq_len(nx))
  logd <- wlds[id_max_seq] + log(.colSums(tcrossprod(w, 1/w[id_max]) * exp(t.default(t.default(lds) - lds[id_max_seq])), m = K, n = nx, na.rm = FALSE))
  # if (!isTRUE(all.equal.numeric(exp(logd), d))) stop('new dfmx wrong?')
  return(logd)
  
}





# not compute intensive
#' @rdname dfmx
#' @importFrom TukeyGH77 pGH gh2z
#' @export
pfmx <- function(q, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, ..., lower.tail = TRUE, log.p = FALSE) { # not compute-intensive
  if (K == 1L) {
    
    p <- switch(distname, 
           gamma = pgamma(q = q, shape = pars[,1L], scale = pars[,2L]),
           nbinom = pnbinom(q = q, size = pars[,1L], prob = pars[,2L]),
           norm = pnorm(q = q, mean = pars[,1L], sd = pars[,2L]),
           GH = pGH(q = q, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], ...),
           sn = {
             psn(x = q, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]) # first parameter is `x`, not `q`
           },
           st = {
             pst(x = q, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L]) # first parameter is `x`, not `q`
           },
           stop('I do not have `p', distname, '` function'))
    
  } else {
    
    qM_naive <- tcrossprod(rep(1, times = K), q)
    ps <- switch(distname, gamma = {
      pgamma(qM_naive, shape = pars[,1L], scale = pars[,2L])
    }, nbinom = {
      pnbinom(qM_naive, size = pars[,1L], prob = pars[,2L])
    }, sn = {
      # ?sn::psn does not respect `attr(x, 'dim')`, but do handle \link[base]{matrix} `x` correctly
      # packageDate('sn') 2023-04-04: *sometimes* get error by using matrix `x`, dont know why
      # tmp2 <- array(psn(qM_naive, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]), dim = dim(qM_naive))
      do.call(rbind, args = lapply(seq_len(K), FUN = \(i) psn(q, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L])))
    }, st = {
      # ?sn::pst does not respect `attr(x, 'dim')`, and do not handle \link[base]{matrix} `x` correctly!!
      do.call(rbind, args = lapply(seq_len(K), FUN = \(i) pst(q, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], nu = pars[i,4L])))
      # tmp2 <- array(psn(qM_naive, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L]), dim = dim(qM_naive))
      # range(tmp - tmp2) # not the same!!
    }, norm = {
      zM <- tcrossprod(1/pars[,2L], q) - pars[,1L]/pars[,2L]
      pnorm(q = zM)
    }, GH = {
      zM <- tcrossprod(1/pars[,2L], q) - pars[,1L]/pars[,2L]
      g <- pars[,3L]
      h <- pars[,4L]
      for (i in seq_len(K)) zM[i,] <- gh2z(q = zM[i,], g = g[i], h = h[i], ...)
      pnorm(q = zM)
    }, stop('I do not have `p', distname, '` function'))
    
    p <- c(crossprod(w, ps))
    
  }
  
  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}


#' @title Obtain `interval` for \link[rstpm2]{vuniroot} for Function [qfmx()] 
#' 
#' @param dist \linkS4class{fmx} object
#' 
#' @param p \link[base]{length}-2 \link[base]{numeric} \link[base]{vector}
#' 
#' @param distname,K,pars,w (optional) ignored if `dist` is provided
#' 
#' @param ... additional parameters, currently not used
#' 
#' @returns 
#' Function [qfmx_interval()] returns a \link[base]{length}-2 \link[base]{numeric} \link[base]{vector}.
#' 
#' @keywords internal
#' @export
qfmx_interval <- function(dist, p = c(1e-6, 1-1e-6), distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, ...) {
  qfun <- paste0('q', distname)
  y_ls <- lapply(seq_len(K), FUN = \(i) {# single component
    iw <- w[i]
    ip <- p
    ip[1L] <- min(p[1L]/iw, .05)
    ip[2L] <- max(1-(1-p[2L])/iw, .95)
    #ip[1L] <- max(p[1L]/iw, .001)
    #ip[2L] <- min(1-(1-p[2L])/iw, .999)
    out <- do.call(qfun, args = c(list(p = ip), as.list.default(pars[i,])))
    if (any(is.infinite(out))) return(invisible()) # very likely to be a malformed estimate
    return(out)
  })
  y <- unlist(y_ls, use.names = FALSE)
  if (anyNA(y)) stop(sQuote(qfun), ' returns NA_real_')
  switch(distname, gamma = {
    c(0, max(y)) # important!  for distribution with lower bound 0
  }, c(min(y), max(y)))
}



#' @rdname dfmx
#' @importFrom TukeyGH77 qGH gh2z
#' @importFrom rstpm2 vuniroot
#' @export
qfmx <- function(p, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, interval = qfmx_interval(dist = dist), ..., lower.tail = TRUE, log.p = FALSE) {
  
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  if (K == 1L) {
    switch(
      EXPR = distname, 
      gamma = return(qgamma(p, shape = pars[,1L], scale = pars[,1L])),
      nbinom = return(qnbinom(p, size = pars[,1L], prob = pars[,2L])),
      norm = return(qnorm(p, mean = pars[,1L], sd = pars[,2L])),
      GH = return(qGH(p, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L])),
      sn = {
        return(qsn(p, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]))
      },
      st = {
        return(qst(p, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L]))
      },
      stop('I do not have `q', distname, '` function'))
  }
  
  t_w <- t.default(w)
  seqid <- seq_len(K) # 'GH' and 'gamma' will use this
  switch(distname, norm =, GH = {
    sdinv <- 1 / pars[,2L] # constant, save time in vuniroot algorithm
    eff <- pars[,1L] * sdinv # effect size
  })
  
  ones <- rep(1, times = K)
  # `f` in \link[rstpm2]{vuniroot} is essentially \link[TukeyGH77]{pfmx} !!
  f <- switch(distname, gamma = {
    shape <- pars[,1L]
    scale <- pars[,2L]
    function(q) {
      c(t_w %*% pgamma(tcrossprod(ones, q), shape = shape, scale = scale)) - p
    }
    # R does not have a function for https://en.wikipedia.org/wiki/Incomplete_gamma_function
    # thus this calculation cannot be speed up, as for now..
    # see ?stats::pgamma for more details
  }, nbinom = {
    size <- pars[,1L]
    prob <- pars[,2L]
    function(q) {
      c(t_w %*% pnbinom(tcrossprod(ones, q), size = size, prob = prob)) - p
    }
  }, sn = {
    xi <- pars[,1L]
    omega <- pars[,2L]
    alpha <- pars[,3L]
    function(q) {
      # ?sn::psn does not respect `attr(q, 'dim')`, but do handle \link[base]{matrix} `x` correctly
      #ps <- array(psn(tcrossprod(ones, q), xi = xi, omega = omega, alpha = alpha), dim = c(K, length(q)))
      ps <- do.call(rbind, args = lapply(seq_len(K), FUN = \(i) psn(q, xi = xi[i], omega = omega[i], alpha = alpha[i])))
      c(t_w %*% ps) - p
    }
  }, st = {
    xi <- pars[,1L]
    omega <- pars[,2L]
    alpha <- pars[,3L]
    nu <- pars[,4L]
    function(q) {
      # ?sn::psn does not respect `attr(q, 'dim')`, and do not handle \link[base]{matrix} `x` correctly !!
      ps <- do.call(rbind, args = lapply(seq_len(K), FUN = \(i) pst(q, xi = xi[i], omega = omega[i], alpha = alpha[i], nu = nu[i])))
      c(t_w %*% ps) - p
    }
  }, norm = {
    function(q) {
      zM <- tcrossprod(sdinv, q) - eff
      c(t_w %*% pnorm(q = zM)) - p
    }
  }, GH = {
    g <- pars[,3L]
    h <- pars[,4L]
    function(q) {
      z <- q0 <- tcrossprod(sdinv, q) - eff
      for (i in seqid) z[i,] <- gh2z(q = q0[i,], g = g[i], h = h[i], ...)
      c(t_w %*% pnorm(q = z)) - p
    }
  }, stop('I do not have `q', distname, '` function'))
  
  return(vuniroot(f = f, lower = interval[1L], upper = interval[2L])[[1L]])
}






#' @rdname dfmx
#' @importFrom sn rsn rst
#' @importFrom TukeyGH77 rGH
#' @export
rfmx <- function(n, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w) {
  if (!is.integer(n) || anyNA(n) || length(n) != 1L || n <= 0L) stop('sample size must be len-1 positive integer (e.g., use 100L instead of 100)')
  id <- if (K == 1L) {
    rep(1L, times = n) 
    # important! do not disturb the random seed when `K = 1L`
  } else sample.int(n = K, size = n, replace = TRUE, prob = w)
  d2 <- cbind(pars, n = tabulate(id, nbins = K)) # 'matrix'
  r <- paste0('r', distname)
  
  K |>
    seq_len() |>
    lapply(FUN = \(i) {
      d2[i, ] |>
        as.list.default() |> 
        do.call(what = r)
    }) |> 
    unlist(use.names = FALSE)

}














