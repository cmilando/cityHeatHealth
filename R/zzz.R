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
