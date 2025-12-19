#' Scripts to clean the exposure data
#'
#' things to be aware of
# - there may be missing data, so fill in with na.approx
# - you need to include the full year so you can do xbasis matrix
# - you may need to know the mapping so you can average up
#' @param data
#' @importFrom dplyr lag
#' @importFrom dplyr expand_grid
#' @importFrom dplyr left_join
#' @importFrom dplyr join_by
#' @importFrom dplyr expand_grid
#' @importFrom zoo na.approx
#' @importFrom lubridate month
#' @importFrom lubridate year
#' @returns
#' @export
#'
#' @examples
make_exposure_matrix <- function(data, warm_season_months = 5:9) {

  # *************
  # THINGS TO FIX
  # - switch to data.table not dplyr
  # - names -- need to have the mapping be for TOWN20, date, tmax_C
  # - don't use lubridate for month
  # - warning if the dates are not the whole year because then
  #   there will be NAs
  # - dont hardcode the number of lags
  # - don't hardcode the name of the lagcolumn
  # - give a warning if NAs are filled and the gap is large
  # - Do the spatial averaging if you want to map to a larger unit
  #   >> read this from data, which will have an attribute that says
  #      `geo_unit` and `geo_unit_grp`, so the default is `geo_unit`
  # *************

  # get years froom data
  years <- sort(unique(lubridate::year(data$date)))

  # make the skeleton you need later
  get_dt <- function(yy) {
    st = make_date(yy, 1, 1)
    ed = make_date(yy, 12, 31)
    dt = seq.Date(st, ed, by = 'day')
    return(dt)
  }

  all_dt <- lapply(years, get_dt)

  all_dt <- do.call(c, all_dt)

  # this is where xgrid comes in
  unique_areas <- unique(data$TOWN20)
  xgrid <- expand_grid(date = all_dt, TOWN20 = unique_areas)

  # and now left_join
  xgrid <- left_join(xgrid, data, by = join_by(date, TOWN20))

  # need to fill in NA values
  exposure1_l <- split(xgrid, f = xgrid$TOWN20)

  exposure2_l <- vector("list", length(exposure1_l))

  for(i in 1:length(exposure1_l)) {

    x <- exposure1_l[[i]]

    ev1 <- is.na(x$tmax_C)

    # give a warning here
    if(all(ev1)) {
      warning("all exposures for TOWN20 X were blank -- skipping")
      next
    }

    # fix the edge case where either the first or last one
    i = 1
    cont = F
    while(is.na(x$tmax_C[i])) {
      i = i+1
      cont = T
    }
    if(cont) {
      fillx <- x$tmax_C[i]
      x$tmax_C[1:(i - 1)] <- fillx
    }

    i = nrow(x)
    cont = F
    while(is.na(x$tmax_C[i])) {
      i = i-1
      cont = T
    }
    if(cont) {
      fillx <- x$tmax_C[i]
      x$tmax_C[(i + 1):nrow(x)] <- fillx
    }

    # this should have caught the edge cases, cause a failure if not
    if(is.na(x$tmax_C[1])) stop()
    if(is.na(x$tmax_C[nrow(x)])) stop()

    # # has to be arranged by date so the lags make sense
    x <- x[order(x$date), ]

    # fill in any NAs
    x$tmax_C  <- zoo::na.approx(x$tmax_C, maxgap = 5)
    if(any(is.na(x$tmax_C))) stop()

    # have to do this here so it doesn't bleed over
    # you can hard-code this so it doesn't look weird
    x$Templag1 <- dplyr::lag(x$tmax_C, 1)
    x$Templag2 <- dplyr::lag(x$tmax_C, 2)
    x$Templag3 <- dplyr::lag(x$tmax_C, 3)
    x$Templag4 <- dplyr::lag(x$tmax_C, 4)
    x$Templag5 <- dplyr::lag(x$tmax_C, 5)
    x$Templag6 <- dplyr::lag(x$tmax_C, 6)
    x$Templag7 <- dplyr::lag(x$tmax_C, 7)
    x$Templag8 <- dplyr::lag(x$tmax_C, 8)

    exposure2_l[[i]] <- x
  }

  # rbind them all together
  exposure2 <- do.call(rbind, exposure2_l)

  # get just the warm season months
  rr <- lubridate::month(exposure2$date) %in% warm_season_months
  warm_season_exposure <- exposure2[rr, ]

  # set the class as an exposure
  class(warm_season_exposure) <- c(class(warm_season_exposure), "exposure")

  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(warm_season_exposure))) {
    warning("some NAs persist, investigate")
  }

  return(warm_season_exposure)

}
