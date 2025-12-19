#' Run a conditional poisson model for a set of parameters
#'
#'
#' @import data.table
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom gnm gnm
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes a data.table of outcomes, created by `make_outcome_table`
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#'
#' @returns
#' @export
#'
#' @examples
single_run <- function(exposure_matrix, outcomes_tbl,
                        argvar = NULL, arglag = NULL, maxlag = NULL,
                       min_n = 50) {

  ## ******************
  ## TODO
  ## - I want to have defaults for argvar and arglag
  ## - tmaxC has to not be hard-coded, you'll figure this tou
  ## - TOWN20 has to be renamed, and flexible to what the strata are
  ## - and keep in mind you'll be working on different resolutions
  ## - give a warning that the strata are fixed at day, month, year
  ##   and area unit. this is the default unless you really want to
  ##   do things differently
  ## ******************

  warning("check that both inputs are the right class of variables")

  # well
  warning("check geo_unit is the same for both")

  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }


  warning("set minN and check it for this single_run")
  # you also need a place to define MinN? is it here or not at all, just
  # a warning since that is what is defined later


  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure

  if(is.null(maxlag)) {
    maxlag = 5
  }

  if(is.null(argvar)) {
    x_knots = quantile(exposure_matrix[, get(exposure_col)], probs = c(0.5, 0.9))
    argvar <- list(fun = 'ns', knots = x_knots)
  }

  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = dlnm::logknots(maxlag, nk = 2))
  }

  ## get the columns you need
  xcols <- c(exposure_col, paste0('explag',1:maxlag))
  x_mat <- exposure_matrix[, ..xcols]

  ## probably should make sure that exposure_matrix and outcome_tbl
  ## are the same size, at least
  stopifnot(dim(x_mat)[1] == dim(outcomes_tbl)[1])



  ## if you are safe to proceed, make the x_mat
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  ## ******************************************
  ## if using GNM, you get COEF and VCOV as part of the model objects
  m_sub <- gnm(daily_deaths ~ cb,
               data = outcomes_tbl,
               family = quasipoisson,
               eliminate = factor(strata),
               subset = keep == 1)

  m_coef <- coef(m_sub)
  m_vcov <- vcov(m_sub)
  ## ******************************************

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = mean(exposure_matrix[, get(exposure_col)]),
                  by = 0.1)

  cen = cp$predvar[which.min(cp$allRRfit)]

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = cen,
                  by = 0.1)

  return(cp)


}
