
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

