#' Make xgrid
#'
#' @param data
#' @param column_mapping
#' @param months
#' @import data.table
#' @importFrom lubridate make_date
#' @importFrom tidyr expand_grid
#' @returns
#'
#' @examples
make_xgrid <- function(data, column_mapping, months = 1:12) {

  setDT(data)
  date_col <- column_mapping$date
  years   <- sort(unique(data[, year(get(date_col))]))

  # make the skeleton you need later
  # this is one of the first key stumbling blocks
  get_dt <- function(yy) {
    st = make_date(yy, 1, 1)
    ed = make_date(yy, 12, 31)
    dt = seq.Date(st, ed, by = 'day')
    mn = month(dt)
    return(dt[mn %in% months])
  }

  all_dt <- lapply(years, get_dt)
  all_dt <- do.call(c, all_dt)

  # ******
  # this is where xgrid comes in
  geo_col <- column_mapping$geo_unit
  unique_areas <- unlist(unique(data[, get(geo_col)]))
  xgrid <- tidyr::expand_grid(date = all_dt, geo_unit = unique_areas)
  names(xgrid) = c(column_mapping$date, column_mapping$geo_unit)
  head(xgrid)

  # ******
  # join back in group level
  geo_cols <- c(
    column_mapping$geo_unit,
    column_mapping$geo_unit_grp
  )

  geo_unit_mapping <- unique(data[, ..geo_cols])

  # i think here, if any geo_unit appears more than once
  # this should enforce
  if(length(unique_areas) != nrow(geo_unit_mapping)) {
    stop("need to enforce that there is a 1:1 mapping of geo_unit to geo_unit_grp.
         That is, a single geo_unit cannot belong to more than one geo_unit_grp")
  }

  # old tidyverse code
  # xgrid <- left_join(xgrid, geo_unit_mapping,
  # by = join_by(!!sym(column_mapping$geo_unit)))

  # new data.table code
  setDT(xgrid)
  setDT(geo_unit_mapping)

  join_col <- column_mapping$geo_unit

  xgrid <- geo_unit_mapping[
    xgrid,
    on = setNames(join_col, join_col)
  ]

  # ******
  # and now join in data
  # xgrid <- left_join(xgrid, data,
  # by = join_by(!!sym(column_mapping$date),
  #              !!sym(column_mapping$geo_unit),
  #              !!sym(column_mapping$geo_unit_grp)))
  join_cols <- c(
    column_mapping$date,
    column_mapping$geo_unit,
    column_mapping$geo_unit_grp
  )

  xgrid <- data[
    xgrid,
    on = setNames(join_cols, join_cols)
  ]

  return(xgrid)

}
