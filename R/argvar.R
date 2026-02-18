#' Internal function to check argvar
#'
#' @param argvar
#' @param this_exp
#'
#' @returns
#'
#' @examples
check_argvar <- function(argvar, this_exp) {

  if(is.null(argvar)) {

    x_knots = quantile(this_exp, probs = c(0.5, 0.9))
    this_argvar <- list(fun = 'ns', knots = x_knots)

  } else {

    # this one is a little tricky because its based on the local data
    if(argvar$fun == 'ns') {

      stopifnot(all(argvar$pct < 1) & all(argvar$pct > 0))
      x_knots = quantile(this_exp, probs = argvar$pct)
      this_argvar <- list(fun = argvar$fun, knots = x_knots)

    } else if(argvar$fun == 'bs') {

      stopifnot(all(argvar$pct < 1) & all(argvar$pct > 0))
      stopifnot(argvar$degree %in% c(2:4))
      x_knots = quantile(this_exp, probs = argvar$pct)
      this_argvar <- list(fun = argvar$fun, knots = x_knots, degree = argvar$degree)

    } else if(argvar$fun == 'strata') {

      stopifnot(
        min(argvar$breaks) >= min(this_exp),
        max(argvar$breaks) <= max(this_exp)
      )
      this_argvar <- list(fun = argvar$fun, breaks = argvar$breaks)

    } else {

      stop("argvar fun that isn't `ns`, `bs`, or `strata` not implemented yet")

    }
  }
  return(this_argvar)
}
