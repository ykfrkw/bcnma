#' Forest plot generic
#'
#' Generic function for forest plots. Dispatches to class-specific methods.
#' For \code{bcnma} objects, see \code{\link{forest.bcnma}}.
#'
#' @param x an object
#' @param ... additional arguments passed to methods
#' @export
forest <- function(x, ...) UseMethod("forest")
