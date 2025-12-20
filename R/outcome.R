#' Function to create the outcome table
#'
#' @param data
#' @param column_mapping a named list that indicates relevant columns in `data`. for the exposure
#' data table, these need to be one of: c('date', "outcome",'factor, 'geo_unit', 'geo_unit_grp')
#' @param warm_season_months the warm season months for this region, default is Northern Hemisphere's
#' May through September (5 through 9)
#' @param collapse_to which factors to collapse across
#' @param grp_level whether to summarize to the group level or not (default)
#' @import data.table
#' @importFrom lubridate days_in_month
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

  # ********
  ## TODO
  #' Collapse outcomes to a given factor level, or overall
  #' This is a good first step and you can deal here with missing data in the record
  #' like adding zeros and what not
  #'
  #'
  #' *********

  # column types
  col_types <- c('date',
                 'factor',
                 'outcome',
                 'geo_unit',
                 'geo_unit_grp')

  warning("need to add checks of column type")


  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' COLLAPSE AND SUMMARIZE
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ##
  setDT(data)

  # **************
  ## first collapse to by summing
  ## both by collapse to and by group
  date_col = column_mapping$date
  geo_unit_col = column_mapping$geo_unit
  geo_unit_grp_col = column_mapping$geo_unit_grp
  outcome_col <- column_mapping$outcome

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

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_grp_col))]

      names(data) <- c(date_col, geo_unit_grp_col, outcome_col)

      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'ALL'
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

      names(data) <- c(date_col, geo_unit_col, geo_unit_grp_col, outcome_col,
                       collapse_to)

      # update the properties here
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_col,
        "geo_unit_grp" = geo_unit_grp_col,
        "factor" = collapse_to
      )
    } else {

      data <- data[,.(
        xoutcome = sum(get(outcome_col))
      ), by = .(get(date_col),
                get(geo_unit_grp_col),
                get(collapse_to))]

      names(data) <- c(date_col, geo_unit_grp_col, outcome_col,
                       collapse_to)

      # update the properties here
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'ALL',
        "factor" = collapse_to
      )


    }
  }

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

  xgrid$strata <- paste0(xgrid[, get(geo_unit_col)],
                        ":yr",yr, ":mn",mn, ":dow",dow, ":", collapse_to)

  # **************
  # Label strata that have no cases, these will be removed later
  # Extract outcome column name programmatically
  group_col   <- "strata"

  # 1. Aggregate by group and create the 'keep' flag
  xgrid_agg <- xgrid[, .(
    keep = as.integer(sum(get(outcome_col), na.rm = TRUE) > 0)
  ), by = group_col]

  # 2. Join back to the original data (left join)
  xgrid_comb <- xgrid[xgrid_agg, on = group_col]

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
    warning("some NAs persist, investigate")
  }

  return(xgrid_comb)
}


