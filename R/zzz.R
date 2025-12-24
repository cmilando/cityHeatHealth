#' checks for cmdstan
#' copied from https://github.com/rok-cesnovar/misc/blob/master/democmdstanr/R/zzz.R
#' @param ...
#'
#' @returns
#'
#' @examples
.onLoad <- function(...) {
  cmdstan_version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  if (is.null(cmdstan_version)) {
    stop("No CmdStan installation found. Run cmdstanr::install_cmdstan() to install.", call. = FALSE)
  }
}


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
