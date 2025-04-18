
#' @title Multinomial Probabilities & Logits
#' 
#' @description 
#' 
#' Transformation between the \link[base]{vector}s of multinomial probabilities and logits.
#' 
#' @param p \link[base]{double} \link[base]{vector}, multinomial probabilities
#' 
#' @param q \link[base]{double} \link[base]{vector}, multinomial logits
#' 
#' @details 
#' 
#' Function [pmlogis()] takes a length \eqn{k-1} \link[base]{double} \link[base]{vector} of 
#' multinomial logits \eqn{q} and convert them into length \eqn{k} multinomial probabilities \eqn{p}, 
#' regarding the *first* category as reference.
#' 
#' Function [qmlogis()] takes a length \eqn{k} \link[base]{double} \link[base]{vector} of 
#' multinomial probabilities \eqn{p} and convert them into length \eqn{k-1} multinomial logits \eqn{q}, 
#' regarding the *first* category as reference.
#' 
#' @note
#' This transformation is a generalization of functions \link[stats]{plogis} (i.e., logit to probability)
#' and \link[stats]{qlogis} (i.e., probability to logit).
#' 
#' 
#' @returns 
#' 
#' Function [pmlogis()] returns a \link[base]{double} \link[base]{vector} of multinomial probabilities \eqn{p}.
#' 
#' Function [qmlogis()] returns a \link[base]{double} \link[base]{vector} of multinomial logits \eqn{q}.
#'
#' @examples
#' c(2,5,3) |> qmlogis() |> pmlogis()
#' 
#' # various exceptions
#' c(1, 0) |> qmlogis() |> pmlogis()
#' c(0, 1) |> qmlogis() |> pmlogis()
#' c(0, 0, 1) |> qmlogis() |> pmlogis()
#' c(1, 0, 0, 0) |> qmlogis() |> pmlogis()
#' c(0, 1, 0, 0) |> qmlogis() |> pmlogis()
#' c(0, 0, 1, 0) |> qmlogis() |> pmlogis()
#'    
#' @keywords internal
#' @name mlogis
#' @export
qmlogis <- function(p) {
  lp <- log(p)
  lp[-1L] - lp[1L]
}




#' @rdname mlogis
#' @export
pmlogis <- function(q) {
  #if (anyNA(q)) stop('we *do* allow NA_real_ in logits')
  eq <- exp(q)
  id <- is.infinite(eq)
  ret <- if (any(id)) {
    if (all(id) && all(q < 0)) {
      c(1, q*0)
    } else if (sum(id) > 1L || q[id] < 0) {
      stop('must be single positive Inf')
    } else c(0, id)
  } else {
    y0 <- c(1, eq)
    y0/sum(y0)
  }
  if (!all(is.finite(ret))) {
    #print(q)
    stop('these multinomial logits creates NA or Inf proportions?')
  }
  return(ret)
}



