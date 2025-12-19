#' Function to create the outcome table
#'
#' @param data
#' @param column_mapping
#' @param collapse_to
#' @import data.table
#' @returns
#' @export
#'
#' @examples
make_outcome_table <- function(data,
                               column_mapping,
                               collapse_to = NULL) {

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
                 'geo_unit_grp',
                 'strata_unit')

  ##
  setDT(data)

  # **************
  ## first collapse to by summing
  date_col = column_mapping$date
  geo_unit_col = column_mapping$geo_unit
  geo_unit_grp_col = column_mapping$geo_unit_grp
  outcome_col <- column_mapping$outcome
  join_col <- column_mapping$strata_unit
  stopifnot(join_col %in% c(geo_unit_col, geo_unit_grp_col))

  if(is.null(collapse_to)) {

    collapse_to = 'ALL'

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
      "geo_unit_grp" = geo_unit_grp_col,
      'strata_unit' = join_col
    )

  } else {

    factor_cols <- which(names(column_mapping) == 'factor')
    factor_cols <- unlist(column_mapping[factor_cols])
    stopifnot(collapse_to %in% factor_cols)

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
      'strata_unit' = join_col,
      "factor" = collapse_to
    )


  }

  # **************
  ## create the strata
  dow <- wday(data$date)
  mn  <- month(data$date)
  yr  <- year(data$date)

  data$strata <- paste0(data[, get(join_col)],
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


