

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
#' @param interval length two \link[base]{numeric} \link[base]{vector}, interval for root finding, see \link[TukeyGH77]{vuniroot2} and \link[rstpm2]{vuniroot}
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
#' In such case, we cannot compute the log mixture densities (i.e., `-Inf`), 
#' for the log-likelihood using function [logLik.fmx()].
#' Our solution is to replace these `-Inf` log mixture densities by 
#' the weighted average (using the mixing proportions of `dist`) 
#' of the per-component log densities.
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
#' 
#' @examples 
#' library(ggplot2)
#' 
#' (e1 = fmx('norm', mean = c(0,3), sd = c(1,1.3), w = c(1, 1)))
#' curve(dfmx(x, dist = e1), xlim = c(-3,7))
#' ggplot() + geom_function(fun = dfmx, args = list(dist = e1)) + xlim(-3,7)
#' ggplot() + geom_function(fun = pfmx, args = list(dist = e1)) + xlim(-3,7)
#' hist(rfmx(n = 1e3L, dist = e1), main = '1000 obs from e1')
#' 
#' x = (-3):7
#' round(dfmx(x, dist = e1), digits = 3L)
#' round(p1 <- pfmx(x, dist = e1), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p1, dist = e1), x, tol = 1e-4))
#' 
#' (e2 = fmx('GH', A = c(0,3), g = c(.2, .3), h = c(.2, .1), w = c(2, 3)))
#' ggplot() + geom_function(fun = dfmx, args = list(dist = e2)) + xlim(-3,7)
#' 
#' round(dfmx(x, dist = e2), digits = 3L)
#' round(p2 <- pfmx(x, dist = e2), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p2, dist = e2), x, tol = 1e-4))
#' 
#' (e3 = fmx('GH', g = .2, h = .01)) # one-component Tukey
#' ggplot() + geom_function(fun = dfmx, args = list(dist = e3)) + xlim(-3,5)
#' set.seed(124); r1 = rfmx(1e3L, dist = e3); 
#' set.seed(124); r2 = TukeyGH77::rGH(n = 1e3L, g = .2, h = .01)
#' stopifnot(identical(r1, r2)) # but ?rfmx has much cleaner code
#' round(dfmx(x, dist = e3), digits = 3L)
#' round(p3 <- pfmx(x, dist = e3), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p3, dist = e3), x, tol = 1e-4))
#' 
#' 
#' a1 = fmx('GH', A = c(7,9), B = c(.8, 1.2), g = c(.3, 0), h = c(0, .1), w = c(1, 1))
#' a2 = fmx('GH', A = c(6,9), B = c(.8, 1.2), g = c(-.3, 0), h = c(.2, .1), w = c(4, 6))
#' library(ggplot2)
#' (p = ggplot() + 
#'  geom_function(fun = pfmx, args = list(dist = a1), mapping = aes(color = 'g2=h1=0')) + 
#'  geom_function(fun = pfmx, args = list(dist = a2), mapping = aes(color = 'g2=0')) + 
#'  xlim(3,15) + 
#'  scale_y_continuous(labels = scales::percent) +
#'  labs(y = NULL, color = 'models') +
#'  coord_flip())
#' p + theme(legend.position = 'none')
#' 
#' 
#' @name dfmx
#' @import stats
#' @importFrom sn dsn psn qsn rsn dst pst qst rst
#' @importFrom TukeyGH77 .dGH
#' @importFrom VGAM dgenpois1 pgenpois1 qgenpois1 rgenpois1
#' @export
dfmx <- function(x, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, ..., log = FALSE) {
  if (K == 1L) { # no mixture required!!
    switch(distname, 
           gamma = return(dgamma(x = x, shape = pars[,1L], scale = pars[,2L], log = log)),
           genpois1 = return(dgenpois1(x = x, meanpar = pars[,1L], dispind = pars[,2L], log = log)),
           nbinom = return(dnbinom(x = x, size = pars[,1L], prob = pars[,2L], log = log)),
           norm = return(dnorm(x = x, mean = pars[,1L], sd = pars[,2L], log = log)), 
           GH = return(.dGH(x = x, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], log = log, ...)),
           sn = return(dsn(x = x, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], log = log)),
           st = return(dst(x = x, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L], log = log)),
           stop('I do not have `d', distname, '` function'))
  }
  
  xm <- tcrossprod(rep(1, times = K), x)
  
  lds <- switch(distname, # `lds` is per-component log-densities
                gamma = dgamma(x = xm, shape = pars[,1L], scale = pars[,2L], log = TRUE),
                genpois1 = dgenpois1(x = xm, meanpar = pars[,1L], dispind = pars[,2L], log = TRUE),
                nbinom = dnbinom(x = xm, size = pars[,1L], prob = pars[,2L], log = TRUE),
                norm = dnorm(x = xm, mean = pars[,1L], sd = pars[,2L], log = TRUE), 
                GH = .dGH(x = xm, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], log = TRUE, ...),
                sn = {
                  # ?sn::dsn does not respect `attr(x, 'dim')`
                  # to make things worse, ?sn::dsn does not even handle \link[base]{matrix} `x` correctly!!
                  # array(dsn(x = xm, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], log = TRUE), dim = dim(xm), dimnames = dimnames(xm))
                  # this is wrong!!
                  # have to go the stupid way!!
                  do.call(rbind, args = lapply(seq_len(K), FUN = function(i) dsn(x = x, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], log = TRUE)))
                }, 
                st = {
                  # ?sn::dst gives error on vector `nu`
                  do.call(rbind, args = lapply(seq_len(K), FUN = function(i) dst(x = x, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], nu = pars[i,4L], log = TRUE)))
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
#' @importFrom TukeyGH77 pGH GH2z
#' @export
pfmx <- function(q, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, ..., lower.tail = TRUE, log.p = FALSE) { # not compute-intensive
  if (K == 1L) {
    switch(distname, 
           gamma = return(pgamma(q = q, shape = pars[,1L], scale = pars[,2L], lower.tail = lower.tail, log.p = log.p)),
           genpois1 = {
             # ?VGAM::pgenpois1 do not have `log.p` argument
             p <- pgenpois1(q = q, meanpar = pars[,1L], dispind = pars[,2L], lower.tail = lower.tail)
             if (!log.p) return(p)
             return(log(p))
           },
           nbinom = return(pnbinom(q = q, size = pars[,1L], prob = pars[,2L], lower.tail = lower.tail, log.p = log.p)),
           norm = return(pnorm(q = q, mean = pars[,1L], sd = pars[,2L], lower.tail = lower.tail, log.p = log.p)),
           GH = return(pGH(q = q, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], lower.tail = lower.tail, log.p = log.p, ...)),
           sn = {
             p <- psn(x = q, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]) # first parameter is `x`, not `q`
             if (!lower.tail) p <- 1 - p
             if (log.p) return(log(p))
             return(p)
           },
           st = {
             p <- pst(x = q, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L]) # first parameter is `x`, not `q`
             if (!lower.tail) p <- 1 - p
             if (log.p) return(log(p))
             return(p)
           },
           stop('I do not have `p', distname, '` function'))
  }
  
  qM_naive <- tcrossprod(rep(1, times = K), q)
  ps <- switch(distname, gamma = {
    pgamma(qM_naive, shape = pars[,1L], scale = pars[,2L], lower.tail = lower.tail)
  }, genpois1 = {
    pgenpois1(qM_naive, meanpar = pars[,1L], dispind = pars[,2L], lower.tail = lower.tail)
  }, nbinom = {
    pnbinom(qM_naive, size = pars[,1L], prob = pars[,2L], lower.tail = lower.tail)
  }, sn = {
    # ?sn::psn does not respect `attr(x, 'dim')`, but do handle \link[base]{matrix} `x` correctly
    # packageDate('sn') 2023-04-04: *sometimes* get error by using matrix `x`, dont know why
    # tmp2 <- array(psn(qM_naive, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]), dim = dim(qM_naive))
    tmp <- do.call(rbind, args = lapply(seq_len(K), FUN = function(i) psn(q, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L])))
    if (!lower.tail) 1 - tmp else tmp
  }, st = {
    # ?sn::pst does not respect `attr(x, 'dim')`, and do not handle \link[base]{matrix} `x` correctly!!
    tmp <- do.call(rbind, args = lapply(seq_len(K), FUN = function(i) pst(q, xi = pars[i,1L], omega = pars[i,2L], alpha = pars[i,3L], nu = pars[i,4L])))
    # tmp2 <- array(psn(qM_naive, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L], nu = pars[,4L]), dim = dim(qM_naive))
    # range(tmp - tmp2) # not the same!!
    if (!lower.tail) 1 - tmp else tmp
  }, norm = {
    zM <- tcrossprod(1/pars[,2L], q) - pars[,1L]/pars[,2L]
    pnorm(q = zM, lower.tail = lower.tail)
  }, GH = {
    zM <- tcrossprod(1/pars[,2L], q) - pars[,1L]/pars[,2L]
    gs <- pars[,3L]
    hs <- pars[,4L]
    for (i in seq_len(K)) zM[i,] <- GH2z(q0 = zM[i,], g = gs[i], h = hs[i], ...)
    pnorm(q = zM, lower.tail = lower.tail)
  }, stop('I do not have `p', distname, '` function'))
  
  p <- c(crossprod(w, ps))
  if (log.p) return(log(p)) 
  return(p)
}


#' @title Obtain `interval` for \link[TukeyGH77]{vuniroot2} for Function [qfmx()] 
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
  y_ls <- lapply(seq_len(K), FUN = function(i) {# single component
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
#' @importFrom TukeyGH77 qGH GH2z vuniroot2
#' @export
qfmx <- function(p, dist, distname = dist@distname, K = dim(pars)[1L], pars = dist@pars, w = dist@w, interval = qfmx_interval(dist = dist), ..., lower.tail = TRUE, log.p = FALSE) {
  # if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval)) stop('masked to save time')
  if (log.p) p <- exp(p)
  
  if (K == 1L) {
    switch(
      EXPR = distname, 
      gamma = return(qgamma(p, shape = pars[,1L], scale = pars[,1L], lower.tail = lower.tail, log.p = FALSE)),
      genpois1 = {
        # VGAM::qgenpois1 do not have arguments `lower.tail` and `log.p`
        if (log.p) p <- exp(p)
        if (!lower.tail) p <- 1 - p
        return(qgenpois1(p, meanpar = pars[,1L], dispind = pars[,1L]))
      },
      nbinom = return(qnbinom(p, size = pars[,1L], prob = pars[,2L], lower.tail = lower.tail, log.p = FALSE)),
      norm = return(qnorm(p, mean = pars[,1L], sd = pars[,2L], lower.tail = lower.tail, log.p = FALSE)),
      GH = return(qGH(p, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], lower.tail = lower.tail, log.p = FALSE)),
      sn = {
        if (lower.tail) p <- 1 - p
        return(qsn(p, xi = pars[,1L], omega = pars[,2L], alpha = pars[,3L]))
      },
      st = {
        if (lower.tail) p <- 1 - p
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
  # `f` in \link[TukeyGH77]{vuniroot2} is essentially \link[TukeyGH77]{pfmx} !!
  f <- switch(distname, gamma = {
    shape <- pars[,1L]
    scale <- pars[,2L]
    function(q) c(t_w %*% pgamma(tcrossprod(ones, q), shape = shape, scale = scale, lower.tail = lower.tail))
    # R does not have a function for https://en.wikipedia.org/wiki/Incomplete_gamma_function
    # thus this calculation cannot be speed up, as for now..
    # see ?stats::pgamma for more details
  }, genpois1 = {
    meanpar <- pars[,1L]
    dispind <- pars[,2L]
    function(q) c(t_w %*% pgenpois1(tcrossprod(ones, q), meanpar = meanpar, dispind = dispind, lower.tail = lower.tail))
  }, nbinom = {
    size <- pars[,1L]
    prob <- pars[,2L]
    function(q) c(t_w %*% pnbinom(tcrossprod(ones, q), size = size, prob = prob, lower.tail = lower.tail))
  }, sn = {
    xi <- pars[,1L]
    omega <- pars[,2L]
    alpha <- pars[,3L]
    function(q) {
      # ?sn::psn does not respect `attr(q, 'dim')`, but do handle \link[base]{matrix} `x` correctly
      #ps <- array(psn(tcrossprod(ones, q), xi = xi, omega = omega, alpha = alpha), dim = c(K, length(q)))
      ps <- do.call(rbind, args = lapply(seq_len(K), FUN = function(i) psn(q, xi = xi[i], omega = omega[i], alpha = alpha[i])))
      if (!lower.tail) ps <- 1 - ps
      c(t_w %*% ps)
    }
  }, st = {
    xi <- pars[,1L]
    omega <- pars[,2L]
    alpha <- pars[,3L]
    nu <- pars[,4L]
    function(q) {
      # ?sn::psn does not respect `attr(q, 'dim')`, and do not handle \link[base]{matrix} `x` correctly !!
      ps <- do.call(rbind, args = lapply(seq_len(K), FUN = function(i) pst(q, xi = xi[i], omega = omega[i], alpha = alpha[i], nu = nu[i])))
      if (!lower.tail) ps <- 1 - ps
      c(t_w %*% ps)
    }
  }, norm = {
    function(q) {
      zM <- tcrossprod(sdinv, q) - eff
      c(t_w %*% pnorm(q = zM, lower.tail = lower.tail))
    }
  }, GH = {
    gs <- pars[,3L]
    hs <- pars[,4L]
    function(q) {
      z <- q0 <- tcrossprod(sdinv, q) - eff
      for (i in seqid) z[i,] <- GH2z(q0 = q0[i,], g = gs[i], h = hs[i], ...)
      c(t_w %*% pnorm(q = z, lower.tail = lower.tail))
    }
  }, stop('I do not have `q', distname, '` function'))
  
  return(vuniroot2(y = p, f = f, interval = interval))
}






#' @rdname dfmx
#' @examples
#' # to use [rfmx] without \pkg{fmx}
#' (d = fmx(distname = 'GH', A = c(-1,1), B = c(.9,1.1), g = c(.3,-.2), h = c(.1,.05), w = c(2,3)))
#' d@@pars
#' set.seed(14123); x = rfmx(n = 1e3L, dist = d)
#' set.seed(14123); x_raw = rfmx(n = 1e3L,
#'  distname = 'GH', K = 2L,
#'  pars = rbind(
#'   c(A = -1, B = .9, g = .3, h = .1),
#'   c(A = 1, B = 1.1, g = -.2, h = .05)
#'  ), 
#'  w = c(.4, .6)
#' )
#' stopifnot(identical(x, x_raw))
#' 
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
  r_fn <- paste0('r', distname)
  xs <- lapply(seq_len(K), FUN = function(i) {
    do.call(what = r_fn, args = as.list.default(d2[i, ]))
  })
  out <- unlist(xs, use.names = FALSE)
  return(out)
}














