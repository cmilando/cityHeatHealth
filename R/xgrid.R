#' An internal function to make the xgrid
#'
#' @param data
#' @param column_mapping
#' @param months
#' @param dt_by either by day or by week
#' @import data.table
#' @importFrom lubridate make_date
#' @importFrom tidyr expand_grid
#' @returns
#'
#' @examples
make_xgrid <- function(data, column_mapping, months = 1:12,
                       dt_by = 'day') {

  #
  setDT(data)

  stopifnot(dt_by %in% c('day', 'week'))

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' GET ALL DATES
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  date_col <- column_mapping$date
  years   <- sort(unique(data[, year(get(date_col))]))

  # make the skeleton you need later
  # this is one of the first key stumbling blocks
  get_dt <- function(yy) {
    st = make_date(yy, 1, 1)
    ed = make_date(yy, 12, 31)
    dt = seq.Date(st, ed, by = dt_by)
    mn = month(dt)
    return(as.IDate(dt[mn %in% months]))
  }

  all_dt <- lapply(years, get_dt)
  all_dt <- do.call(c, all_dt)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CREATE XGRID SKELETON
  #'
  #' For dates, you need all the dates
  #'
  #' For outcomes, you also need all the dates within strata
  #' that have SOME data so you can add 0s to low data places
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  geo_col <- column_mapping$geo_unit
  unique_areas <- unlist(unique(data[, get(geo_col)]))

  if(any(names(column_mapping) == 'factor')) {

    stopifnot(length(which(names(column_mapping) == 'factor')) == 1)

    factor_col <- column_mapping$factor

    unique_fcts <- unlist(unique(data[, get(factor_col)]))

    xgrid <- tidyr::expand_grid(date = all_dt,
                                geo_unit = unique_areas,
                                fct = unique_fcts)

    names(xgrid) = c(column_mapping$date,
                     column_mapping$geo_unit,
                     column_mapping$factor)

  } else {

    xgrid <- tidyr::expand_grid(date = all_dt, geo_unit = unique_areas)

    names(xgrid) = c(column_mapping$date, column_mapping$geo_unit)

  }

  setDT(xgrid)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' JOIN WITH GROUP DATA
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  geo_cols <- c(
    column_mapping$geo_unit,
    column_mapping$geo_unit_grp
  )

  geo_unit_mapping <- unique(data[, ..geo_cols])

  # if any geo_unit appears more than once then its not 1:1
  if(length(unique_areas) != nrow(geo_unit_mapping)) {
    stop("need to enforce that there is a 1:1 mapping of geo_unit to geo_unit_grp.
         That is, a single geo_unit cannot belong to more than one geo_unit_grp")
  }

  setDT(geo_unit_mapping)

  join_col <- column_mapping$geo_unit

  xgrid <- geo_unit_mapping[
    xgrid,
    on = setNames(join_col, join_col)
  ]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' FINALLY - JOIN WITH DATA
  #'
  #' This shouldn't add any rows
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(any(names(column_mapping) == 'factor')) {
    join_cols <- c(
      column_mapping$date,
      column_mapping$geo_unit,
      column_mapping$geo_unit_grp,
      column_mapping$factor
    )
  } else {
    join_cols <- c(
      column_mapping$date,
      column_mapping$geo_unit,
      column_mapping$geo_unit_grp
    )
  }

  xgrid <- data[
    xgrid,
    on = setNames(join_cols, join_cols)
  ]

  return(xgrid)

}
