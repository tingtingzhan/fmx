


#' @title Maximum a Posteriori clustering
#' 
#' @description ..
#' 
#' @param dist an \linkS4class{fmx} object
#' 
#' @param x \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @returns 
#' Function [MaP()] returns an \link[base]{integer} \link[base]{vector}.
#' 
#' @examples 
#' x = rnorm(1e2L, sd = 2)
#' m = fmx('norm', mean = c(-1.5, 1.5), w = c(1, 2))
#' library(ggplot2)
#' ggplot() + geom_function(fun = dfmx, args = list(dist = m)) + 
#'   geom_point(mapping = aes(x = x, y = .05, color = factor(MaP(x, dist = m)))) + 
#'   labs(color = 'Maximum a Posteriori\nClustering')
#' 
#' @keywords internal
#' @export
MaP <- function(x, dist, ...) {
  d <- dfmx(x = x, dist = dist)
  ret <- t.default(attr(d, which = 'posterior', exact = TRUE)) / d
  max.col(ret)
}

# https://www3.cs.stonybrook.edu/~has/CSE594/Notes/(3)%20Clustering%20and%20Prediction%20(all%20slides).pdf