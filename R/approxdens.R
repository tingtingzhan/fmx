
######################################
## Modify in tzhInternal package!!
######################################




#' @title Empirical Density Function
#' 
#' @description ..
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations
#' 
#' @param ... additional parameters of \link[stats]{density.default}
#' 
#' @details 
#' \link[stats]{approx} inside \link[stats]{density.default}
#' 
#' another 'layer' of \link[stats]{approxfun}
#' 
#' @returns 
#' Function [approxdens()] returns a \link[base]{function}.
#' 
#' @examples 
#' x = rnorm(1e3L)
#' f = approxdens(x)
#' f(x[1:3])
#' 
#' @importFrom stats density.default approxfun
#' @export
approxdens <- function(x, ...) {
  if (anyNA(x)) stop('input must not contain missing value')
  kern <- density.default(x, ...)
  approxfun(x = kern$x, y = kern$y, method = 'linear', rule = 1) 
}
