#' Run a conditional poisson model for a single zone

#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr join_by
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom gnm gnm
#' @importFrom lubridate month
#' @importFrom lubridate year
#' @importFrom lubridate wday
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes a dataframe of date, and outcome
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#'
#' @returns
#' @export
#'
#' @examples
single_zone <- function(exposure_matrix, outcomes,
                        argvar = NULL, arglag = NULL, maxlag = NULL) {

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

  if(is.null(maxlag)) {
    maxlag = 5
  }

  if(is.null(argvar)) {
    x_knots = quantile(exposure_matrix$tmax_C, probs = c(0.5, 0.9))
    argvar <- list(fun = 'ns', knots = x_knots)
  }

  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = dlnm::logknots(maxlag, nk = 2))
  }

  ## keep only the columns you need
  x_mat <- exposure_matrix[, c('tmax_C', paste0('Templag',1:maxlag))]

  ## create the crossbasis
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  ## create the strata
  outcomes$dow   <- lubridate::wday(outcomes$date, label = T)
  outcomes$month <- lubridate::month(outcomes$date, label = T)
  outcomes$year  <- lubridate::year(outcomes$date)

  outcomes$strata <- paste0(outcomes$TOWN20, ":",
                            outcomes$year, ":",
                            outcomes$month, ":",
                            outcomes$dow)

  ## get rid of 0 strata
  ## you make sure there are no empty strata
  outcomes_agg <- outcomes %>%
    group_by(strata) %>%
    summarize(
      .groups = 'keep',
      total_daily_deaths = sum(daily_deaths)
    ) %>%
    mutate(keep = ifelse(total_daily_deaths > 0, 1, 0))

  outcomes_comb <- left_join(outcomes, outcomes_agg, by = join_by(strata))

  ## ******************************************
  ## if using GNM, you get COEF and VCOV as part of the model objects
  m_sub <- gnm(daily_deaths ~ cb,
               data = outcomes_comb,
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
                  cen = mean(exposure_matrix$tmax_C),
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
