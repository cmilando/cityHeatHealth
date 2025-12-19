#' Script to clean the exposure data
#'
#' things to be aware of
#' - there may be missing data, so fill in with na.approx
#' - you need to include the full year so you can do xbasis matrix
#' - you may need to know the mapping so you can average up
#'
#' @param data a dataset of exposures
#' @param column_mapping a named list that indicates relevant columns
#' @param warm_season_months the warm season months for this region, default is Northern Hemisphere's 5 through 9
#' @param maxgaps the maximum allowable gap, to be passed to zoo::na.approx
#' @import data.table
#' @importFrom dplyr lag
#' @importFrom tidyr expand_grid
#' @importFrom zoo na.approx
#' @importFrom lubridate make_date
#' @returns
#' @export
#'
#' @examples
make_exposure_matrix <- function(data, column_mapping,
                                 warm_season_months = 5:9,
                                 maxgap = 5) {

  # *************
  # THINGS TO FIX

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

  # column types
  col_types <- c('date', "exposure", 'geo_unit', 'geo_unit_grp')

  # check that all the types are valid
  if(!all(names(column_mapping) %in% col_types))
    stop('Names of column mapping is not one of the valid types:
          date, exposure, geo_unit, geo_unit_grp')

  # check that all values are correct
  if(!all(column_mapping %in% names(data)))
    stop('Values of column mapping not matched with column names of data.
          Look for a typo')

  # get years from data
  years <- sort(unique(year(data[, column_mapping$date])))

  # make the skeleton you need later
  get_dt <- function(yy) {
    st = lubridate::make_date(yy, 1, 1)
    ed = lubridate::make_date(yy, 12, 31)
    dt = seq.Date(st, ed, by = 'day')
    return(dt)
  }

  all_dt <- lapply(years, get_dt)
  all_dt <- do.call(c, all_dt)

  # this is where xgrid comes in
  unique_areas <- unique(data[, column_mapping$geo_unit])
  xgrid <- tidyr::expand_grid(date = all_dt, geo_unit = unique_areas)
  names(xgrid) = c(column_mapping$date, column_mapping$geo_unit)

  # join back in county
  town_mapping <- unique(data[, c(column_mapping$geo_unit,
                                  column_mapping$geo_unit_grp)])

  # old tidyverse code
  # xgrid <- left_join(xgrid, town_mapping,
  # by = join_by(!!sym(column_mapping$geo_unit)))

  # new data.table code
  setDT(xgrid)
  setDT(town_mapping)

  join_col <- column_mapping$geo_unit

  xgrid <- town_mapping[
    xgrid,
    on = setNames(join_col, join_col)
  ]

  # and now join in data
  # xgrid <- left_join(xgrid, data, by = join_by(!!sym(column_mapping$date),
  #                                              !!sym(column_mapping$geo_unit),
  #                                              !!sym(column_mapping$geo_unit_grp)))
  setDT(data)

  join_cols <- c(
    column_mapping$date,
    column_mapping$geo_unit,
    column_mapping$geo_unit_grp
  )

  xgrid <- data[
    xgrid,
    on = setNames(join_cols, join_cols)
  ]

  # need to fill in NA values
  exposure1_l <- split(xgrid, f = xgrid[, ..join_col])
  length(exposure1_l)

  exposure_col <- column_mapping$exposure
  exposure2_l <- vector("list", length(exposure1_l))

  for(i in 1:length(exposure1_l)) {

    x <- exposure1_l[[i]]

    # check for NAs
    ev1 <- is.na(x[, ..exposure_col])

    # give a warning here
    if(all(ev1)) {
      warning(paste0("all exposures for geo_unit ", x[1, ..join_col],
                     " were blank -- skipping"))
      next
    }

    # fix the edge case where either the first or last one
    j = 1
    cont = F
    while(is.na(x[j, ..exposure_col])) {
      j = j + 1
      cont = T
    }
    if(cont) {
      fillx <- x[j, ..exposure_col]
      x[1:(j - 1), ..exposure_col] <- fillx
    }

    j = nrow(x)
    cont = F
    while(is.na(x[j, ..exposure_col])) {
      j = j-1
      cont = T
    }
    if(cont) {
      fillx <- x[j, ..exposure_col]
      x[(j + 1):nrow(x), ..exposure_col] <- fillx
    }

    # this should have caught the edge cases, cause a failure if not
    if(is.na(x[1, ..exposure_col]))
      stop(paste0("first value is still NA despite cleaning for geo_unit ", x[1, ..join_col]))

    if(is.na(x[nrow(x), ..exposure_col]))
      stop(paste0("last value is still NA despite cleaning for geo_unit ", x[1, ..join_col]))

    # # has to be arranged by date so the lags make sense
    x <- x[order(x$date), ]

    # fill in any NAs using zoo::na.approx
    x[, (exposure_col) := zoo::na.approx(get(exposure_col),maxgap = maxgap)]
    if(any(is.na(x[, ..exposure_col]))) stop("zoo::na.approx() was not able to remove all NAs")

    # have to do this here so it doesn't bleed over
    # you can hard-code this so it doesn't look weird
    x$explag1 <- dplyr::lag(x[, ..exposure_col], 1)
    x$explag2 <- dplyr::lag(x[, ..exposure_col], 2)
    x$explag3 <- dplyr::lag(x[, ..exposure_col], 3)
    x$explag4 <- dplyr::lag(x[, ..exposure_col], 4)
    x$explag5 <- dplyr::lag(x[, ..exposure_col], 5)
    x$explag6 <- dplyr::lag(x[, ..exposure_col], 6)
    x$explag7 <- dplyr::lag(x[, ..exposure_col], 7)
    x$explag8 <- dplyr::lag(x[, ..exposure_col], 8)

    exposure2_l[[i]] <- x

  }

  # rbind them all together
  exposure2 <- do.call(rbind, exposure2_l)

  # get just the warm season months
  rr <- month(exposure2$date) %in% warm_season_months
  warm_season_exposure <- exposure2[rr, ]

  # set the class as an exposure
  class(warm_season_exposure) <- c(class(warm_season_exposure), "exposure")

  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(warm_season_exposure))) {
    warning("some NAs persist, investigate")
  }

  return(warm_season_exposure)

}
