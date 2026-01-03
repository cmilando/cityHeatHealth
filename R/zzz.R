#' spatial_plot base class
#'
#' @param x
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
spatial_plot <- function(x, ...) {
  UseMethod("spatial_plot")
}


#' forest_plot base class
#'
#' @param x
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
forest_plot <- function(x, ...) {
  UseMethod("forest_plot")
}
