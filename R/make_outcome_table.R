#' Function to create the outcome table
#'
#' @param data
#' @param column_mapping a named list that indicates relevant columns in `data`. for the exposure
#' data table, these need to be one of: c('date', "outcome",'factor, 'geo_unit', 'geo_unit_grp')
#' @param collapse_to which factors to collapse across
#' @param grp_level whether to summarize to the group level or not (default)
#' @import data.table
#' @returns
#' @export
#'
#' @examples
make_outcome_table <- function(data,
                               column_mapping,
                               collapse_to = NULL,
                               grp_level = FALSE) {

  # ********
  ## TODO
  #' Collapse outcomes to a given factor level, or overall
  #' This is a good first step and you can deal here with missing data in the record
  #' like adding zeros and what not
  #'
  #' you also need a place to define MinN? is it here or not at all, just
  #' a warning since that is what is defined later
  #' *********


  # column types
  col_types <- c('date',
                 'factor',
                 'outcome',
                 'geo_unit',
                 'geo_unit_grp')

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

  # **************
  ## create the strata
  dow <- wday(data$date)
  mn  <- month(data$date)
  yr  <- year(data$date)

  data$strata <- paste0(data[, get(geo_unit_col)],
                        ":yr",yr, ":mn",mn, ":dow",dow, ":", collapse_to)

  # **************
  # Label strata that have no cases, these will be removed later

  # Extract outcome column name programmatically
  group_col   <- "strata"

  # 1. Aggregate by group and create the 'keep' flag
  data_agg <- data[, .(
    keep = as.integer(sum(get(outcome_col), na.rm = TRUE) > 0)
  ), by = group_col]

  # 2. Join back to the original data (left join)
  data_comb <- data[data_agg, on = group_col]

  # **************
  # prepare for output

  # set the class as an exposure
  class(data_comb) <- c(class(data_comb), "outcoome")

  # set attributes
  attr(data_comb, "column_mapping") <- column_mapping


  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(data_comb))) {
    warning("some NAs persist, investigate")
  }

  return(data_comb)
}


