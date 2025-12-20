#' Run a conditional poisson model for a single geographic unit
#'
#'
#' @import data.table
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom dlnm crossreduce
#' @importFrom gnm gnm
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes a data.table of outcomes, created by `make_outcome_table`
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#'
#' @returns
#' @export
#'
#' @examples
condPois_single <- function(exposure_matrix, outcomes_tbl,
                        argvar = NULL, arglag = NULL, maxlag = NULL,
                       min_n = NULL) {

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcomes_tbl))

  ## Check 1.5 -- there should be just one geo_unit
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  exp_geo_unit <- unlist(unique(exposure_matrix[, get(exp_geo_unit_col)]))
  out_geo_unit <- unlist(unique(outcomes_tbl[, get(out_geo_unit_col)]))

  stopifnot(identical(exp_geo_unit, out_geo_unit))
  stopifnot(length(out_geo_unit) == 1)

  ## Check 2
  ## probably should make sure that exposure_matrix and outcome_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  outcome_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  setorderv(
    exposure_matrix,
    c(exp_date_col, exp_geo_unit_col)
  )

  setorderv(
    outcomes_tbl,
    c(outcome_date_col, outcome_geo_unit_col)
  )

  stopifnot(dim(exposure_matrix)[1] == dim(outcomes_tbl)[1])
  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcomes_tbl[, get(outcome_date_col)]))

  # CHECK 4 geo_unit is the same for both"
  stopifnot(all(outcomes_tbl[, get(outcome_geo_unit_col)] %in%
                  exposure_matrix[, get(exp_geo_unit_col)]))

  # CHECK 5
  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }

  ## CHECK 6 - minN
  if(is.null(min_n)) {
    min_n = 50
  }
  outcome_col <- attributes(outcomes_tbl)$column_mapping$outcome
  stopifnot(sum(outcomes_tbl[, get(outcome_col)]) >= min_n)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CREATE CROSSBASIS for this single zone
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////
  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure

  if(is.null(maxlag)) {
    maxlag = 5
  } else {
    warning("check that these are valid")
  }

  if(is.null(argvar)) {
    x_knots = quantile(exposure_matrix[, get(exposure_col)], probs = c(0.5, 0.9))
    argvar <- list(fun = 'ns', knots = x_knots)
  } else {
    # this one is a little tricky because its based on the local data
    if(argvar$fun == 'ns') {
      stopifnot(all(argvar$pct < 1) & all(argvar$pct > 0))
      x_knots = quantile(exposure_matrix[, get(exposure_col)], probs = argvar$pct)
      argvar <- list(fun = argvar$fun, knots = x_knots)
    } else if(argvar$fun == 'bs') {
      stopifnot(all(argvar$pct < 1) & all(argvar$pct > 0))
      stopifnot(argvar$degree %in% c(2:4))
      x_knots = quantile(exposure_matrix[, get(exposure_col)], probs = argvar$pct)
      argvar <- list(fun = argvar$fun, knots = x_knots, degree = argvar$degree)
    } else {
      stop("argvar that isn't `ns` or `bs` not implemented yet")
    }
  }

  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = dlnm::logknots(maxlag, nk = 2))
  } else {
    warning("check that these are valid")
  }

  ## get the columns you need
  xcols <- c(exposure_col, paste0('explag',1:maxlag))
  x_mat <- exposure_matrix[, ..xcols]

  ## if you are safe to proceed, make the x_mat
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  # there should be no NAs here
  if(any(is.na(cb))) stop("crossbasis has NULL, something went wrong")

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' RUN GNM
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## if using GNM, you get COEF and VCOV as part of the model objects
  m_sub <- gnm(daily_deaths ~ cb,
               data = outcomes_tbl,
               family = quasipoisson,
               eliminate = factor(strata),
               subset = keep == 1)

  m_coef <- coef(m_sub)
  m_vcov <- vcov(m_sub)

  # there should be no NAs
  if(any(is.na(m_coef))) stop("coef has NULL, something went wrong")
  if(any(is.na(m_vcov))) stop("vcov has NULL, something went wrong")

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT OBJECTS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit
  geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp

  this_geo_unit <- unlist(unique(outcomes_tbl[, get(geo_unit_col)]))
  this_geo_unit_grp <- unlist(unique(outcomes_tbl[, get(geo_unit_grp_col)]))

  exp_mean = mean(exposure_matrix[, get(exposure_col)])
  exp_IQR = IQR(exposure_matrix[, get(exposure_col)])

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = exp_mean,
                  by = 0.1)

  cen = cp$predvar[which.min(cp$allRRfit)]

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = cen,
                  by = 0.1)

  cr <- crossreduce(cb,
                    coef = m_coef,
                    vcov = m_vcov,
                    model.link = "log",
                    cen = cen,     # not necessary but avoids messages
                    by = 0.1)

  return(list(geo_unit = this_geo_unit,
              geo_unit_grp = this_geo_unit_grp,
              cr = cr,
              exp_mean = exp_mean,
              exp_IQR = exp_IQR))


}
