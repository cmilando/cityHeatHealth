#' Function to create the outcome table
#'
#' @param data
#' @param column_mapping a named list that indicates relevant columns in `data`. for the exposure
#' data table, these need to be one of: c('date', "outcome",'factor, 'geo_unit', 'geo_unit_grp')
#' @param warm_season_months the warm season months for this region, default is Northern Hemisphere's
#' May through September (5 through 9)
#' @param collapse_to which factors to collapse across
#' @param grp_level whether to summarize to the group level or not (default)
#'
#' @import data.table
#' @importFrom lubridate days_in_month
#'
#' @returns
#' @export
#'
#' @examples
make_outcome_table <- function(data,
                               column_mapping,
                               warm_season_months = 5:9,
                               collapse_to = NULL,
                               grp_level = FALSE) {

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ##
  setDT(data)

  # column types
  col_types <- c('date', 'factor', 'outcome', 'geo_unit', 'geo_unit_grp')

  # check that all the types are valid
  if(!all(names(column_mapping) %in% col_types))
    stop('Names of column mapping is not one of the valid types:
          date, exposure, geo_unit, geo_unit_grp')

  # check that all values are correct
  if(!all(column_mapping %in% names(data)))
    stop('Values of column mapping not matched with column names of data.
          Look for a typo')

  # type checks
  stopifnot(
    inherits(data[[column_mapping$date]], "Date"),
    is.integer(data[[column_mapping$outcome]]),
    is.character(data[[column_mapping$geo_unit]]),
    is.character(data[[column_mapping$geo_unit_grp]]),
    is.character(data[[column_mapping$factor]])
  )

  # at the beginning there shouldn't be any outcomes < 0
  outcome_col <- column_mapping$outcome
  if(any(data[, get(outcome_col)] < 0)) {
    stop("some outcomes < 0, investigate")
  }

  if(any(is.na(data))) {
    warning("check about any NA")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' COLLAPSE AND SUMMARIZE
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # **************
  ## first collapse to by summing
  ## both by collapse to and by group
  date_col         = column_mapping$date
  geo_unit_col     = column_mapping$geo_unit
  geo_unit_grp_col = column_mapping$geo_unit_grp
  outcome_col      = column_mapping$outcome

  # Next check about collapsing across factors
  if(is.null(collapse_to)) {

    collapse_to = 'ALL'

    if(grp_level == FALSE) {

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_col),
                get(geo_unit_grp_col))]

      names(data) <- c(date_col, geo_unit_col, geo_unit_grp_col, outcome_col)

      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_col,
        "geo_unit_grp" = geo_unit_grp_col
      )

    } else {

      warning("make type checks  (e.g., so Date == Date),
         for some reason this doesn't work in some cases? but ok in others?")

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_grp_col))]

      names(data) <- c(date_col, geo_unit_grp_col, outcome_col)

      data$spatial_grp <- 'ALL'

      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'spatial_grp'
      )
    }

  } else {

    factor_cols <- which(names(column_mapping) == 'factor')
    factor_cols <- unlist(column_mapping[factor_cols])
    stopifnot(collapse_to %in% factor_cols)

    if(grp_level == FALSE) {

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_col),
                get(geo_unit_grp_col),
                get(collapse_to))]

      names(data) <- c(date_col, geo_unit_col, geo_unit_grp_col,
                       collapse_to, outcome_col)

      # update the properties here
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_col,
        "geo_unit_grp" = geo_unit_grp_col,
        "factor" = collapse_to
      )
    } else {

      warning("make type checks  (e.g., so Date == Date),
         for some reason this doesn't work in some cases? but ok in others?")

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_grp_col),
                get(collapse_to))]

      names(data) <- c(date_col, geo_unit_grp_col,
                       collapse_to, outcome_col)

      data$spatial_grp <- 'ALL'

      # update the properties here
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'spatial_grp',
        "factor" = collapse_to
      )


    }
  }

  # overwrite date
  data[, (column_mapping$date) :=
         as.Date(
           as.integer(get(column_mapping$date)),
           origin = "1970-01-01"
         )
  ]

  geo_unit_col = column_mapping$geo_unit
  geo_unit_grp_col = column_mapping$geo_unit_grp

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' MAKE XGRID AND STRATA
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## fill in the blanks with 0s
  ## so make xgrid again
  xgrid <- make_xgrid(data, column_mapping, warm_season_months)

  # **************
  ## ADD ZEROS back in given the now full calendar
  rr <- which(is.na(xgrid[[outcome_col]]))
  xgrid[rr, (outcome_col) := 0]

  # **************
  ## create the strata
  dow <- wday(xgrid[, get(date_col)])
  mn  <- month(xgrid[, get(date_col)])
  yr  <- year(xgrid[, get(date_col)])
  if(any(names(column_mapping) == 'factor')) {
    fct <- xgrid[, get(column_mapping$factor)]
    xgrid$strata <- paste0(xgrid[, get(geo_unit_col)],
                           ":yr",yr, ":mn",mn, ":dow",dow,":",
                           column_mapping$factor, fct)

  } else {
    xgrid$strata <- paste0(xgrid[, get(geo_unit_col)],
                           ":yr",yr, ":mn",mn, ":dow",dow)
  }


  # **************
  # Label strata that have no cases, these will be removed later
  # Extract outcome column name programmatically
  group_col   <- "strata"

  # 1. Aggregate by group and create the 'keep' flag
  xgrid_agg <- xgrid[, .(
    strata_total = round(sum(get(outcome_col)))
  ), by = group_col]

  # 2. Join back to the original data (left join)
  xgrid_comb <- xgrid[xgrid_agg, on = group_col]

  # remove any empty strata
  # again this probably could be done better by only adding 0s to
  # strata that have data, thus avoiding the bump in memory, but
  # i don't think that it will matter .... right ?
  # TODO
  rr <- which(xgrid_comb$strata_total > 0)
  xgrid_comb <- xgrid_comb[rr, ,drop = FALSE]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # reset the order
  join_col <- column_mapping$geo_unit
  date_col <- column_mapping$date
  setorderv(
    xgrid_comb,
    c(date_col, join_col)
  )

  # set the class as an exposure
  class(xgrid_comb) <- c(class(xgrid_comb), "outcome")

  # set attributes
  attr(xgrid_comb, "column_mapping") <- column_mapping

  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(xgrid_comb))) {
    warning("some NAs persist, investigate and submit a Github issue :) ")
  }

  return(xgrid_comb)
}


