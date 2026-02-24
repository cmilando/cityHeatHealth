#' Function to create the outcome table
#' TO DO : EDIT XGRID
#' @param data
#' @param column_mapping a named list that indicates relevant columns in `data`. for the exposure
#' data table, these need to be one of: c('date', "outcome",'factor, 'geo_unit', 'geo_unit_grp')
#' @param months_subset the warm season months for this region, default is Northern Hemisphere's
#' May through September (5 through 9)
#' @param collapse_to which factors to collapse across
#' @param collapse_is_spatial is collapse a spatial variable
#' @param collapse_is_temporal is collapse a spatial variable
#' @param grp_level whether to summarize to the group level or not (default)
#' @param keep_unit_outcomes if grp_level is true, whether to keep original unit-level outcomes
#' @param dt_by is it daily data, or weekly or ...
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
                               months_subset = 5:9,
                               dt_by = 'day',
                               collapse_to = NULL,
                               collapse_is_spatial = FALSE,
                               collapse_is_temporal = FALSE,
                               grp_level = FALSE,
                               keep_unit_outcomes = FALSE) {

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ##
  stopifnot(nrow(data) > 0)
  setDT(data)

  #
  stopifnot(typeof(column_mapping) == 'list')

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
    is.character(data[[column_mapping$geo_unit_grp]])
  )
  if("factor" %in% names(column_mapping)) {
    stopifnot(is.character(data[[column_mapping$factor]]))
  }

  # overwrite date
  data[, (column_mapping$date) := as.IDate(get(column_mapping$date))]

  # at the beginning there shouldn't be any outcomes < 0
  outcome_col <- column_mapping$outcome
  if(any(data[, get(outcome_col)] < 0)) {
    stop("some outcomes < 0, investigate")
  }

  if(any(is.na(data))) {
    warning("check about any NA")
  }

  # check that all are unique 1:1
  geo_cols <- c(
    column_mapping$geo_unit,
    column_mapping$geo_unit_grp
  )
  unique_geos <- unique(data[, get(column_mapping$geo_unit)])
  unique_geos_and_grps <- unique(data[, ..geo_cols])
  if(length(unique_geos) != nrow(unique_geos_and_grps)) {
    stop("`geo_unit` repeated across multiple `grps`")
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

      # collapse to = NULL --> so this collapses across factors
      # grp_level = FALSE  --> and doesn't summarize to the group level

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

      # collapse to = NULL --> so this collapses across factors
      # grp_level = TRUE  --> and does summarize to the group level
      #
      # warning("make type checks  (e.g., so Date == Date),
      #    for some reason this doesn't work in some cases? but ok in others?")


      if(keep_unit_outcomes == FALSE) {

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

      } else {

        data <- data[,.(
          xoutcome = sum(get(outcome_col))
        ), by = .(get(date_col),
                  get(geo_unit_col),
                  get(geo_unit_grp_col))]

        names(data) <- c(date_col, geo_unit_col, geo_unit_grp_col, outcome_col)

        # and SKIP column re-mapping until later
        # except remove factor
        rr <- which(names(column_mapping) == 'factor')
        column_mapping[rr] <- NULL

      }

    }

  } else {

    factor_cols <- which(names(column_mapping) == 'factor')
    factor_cols <- unlist(column_mapping[factor_cols])
    stopifnot(collapse_to %in% factor_cols)

    if(grp_level == FALSE) {

      # collapse to = NOT NULL
      # grp_level = FALSE

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

      # collapse to = NOT NULL
      # grp_level = TRUE
      #
      # warning("make type checks  (e.g., so Date == Date),
      #    for some reason this doesn't work in some cases? but ok in others?")

      if(keep_unit_outcomes == FALSE) {
        data <- data[,.(
          xoutcome = sum(get(outcome_col))
        ), by = .(get(date_col),
                  get(geo_unit_grp_col),
                  get(collapse_to))]

        names(data) <- c(date_col, geo_unit_grp_col,
                         collapse_to, outcome_col)

        #
        data$spatial_grp <- 'ALL'

        # update the properties here
        column_mapping <- list(
          "date" = date_col,
          "outcome" = outcome_col,
          "geo_unit" = geo_unit_grp_col,
          "geo_unit_grp" = 'spatial_grp',
          "factor" = collapse_to
        )
      } else {
        data <- data[,.(
          xoutcome = sum(get(outcome_col))
        ), by = .(get(date_col),
                  get(geo_unit_col),
                  get(geo_unit_grp_col),
                  get(collapse_to))]

        names(data) <- c(date_col, geo_unit_col, geo_unit_grp_col,
                         collapse_to, outcome_col)

        # just over-write the collapse_to
        rr <- which(names(column_mapping) == 'factor')
        column_mapping[rr] <- NULL
        column_mapping[['factor']] <- collapse_to

      }

    }
  }

  # overwrite date
  data[, (column_mapping$date) := as.IDate(get(column_mapping$date))]

  geo_unit_col = column_mapping$geo_unit
  geo_unit_grp_col = column_mapping$geo_unit_grp

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' MAKE XGRID AND STRATA
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## fill in the blanks with 0s
  ## so make xgrid again
  xgrid <- make_xgrid(data = data,
                      column_mapping = column_mapping,
                      months_subset = months_subset,
                      dt_by = dt_by,
                      collapse_is_spatial = collapse_is_spatial,
                      collapse_is_temporal = collapse_is_temporal)

  # **************
  ## remove missing values
  ## --> you probably want to set these to 0 I think instead
  rr <- which(is.na(xgrid[[outcome_col]]))
  if(length(rr) > 0) {
    message("Missing values in outcome xgrid were set to 0")

    # default
    xgrid[rr, (outcome_col) := 0]

    # --> could also switch to just removing them altogether
    # xgrid <- xgrid[-rr, ]
  }

  if(any(is.na(xgrid))) {
    message('some NA in outcome xgrid, investigate')
    return(xgrid)
  }

  # **************
  ## create the strata
  dow <- wday(xgrid[, get(date_col)])
  mn  <- month(xgrid[, get(date_col)])
  yr  <- year(xgrid[, get(date_col)])
  xgrid$strata <- paste0(xgrid[, get(geo_unit_col)],
                         ":yr", yr,
                         ":mn", sprintf("%02i", mn),
                         ":dow", sprintf("%02i", dow))

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

  # now that you have a final version, make match_strata
  # also make match strata
  xgrid_comb$match_strata = paste0(
    xgrid_comb[, get(column_mapping$geo_unit)], ":",
    xgrid_comb[, get(date_col)])

  # **************
  # Lastly, re-set column mapping if group-level = TRUE but keep_orig = TRUE
  if(grp_level & keep_unit_outcomes) {
    xgrid_comb$spatial_grp <- 'ALL'
    if(collapse_to == 'ALL') {
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'spatial_grp'
      )
    } else {
      column_mapping <- list(
        "date" = date_col,
        "outcome" = outcome_col,
        "geo_unit" = geo_unit_grp_col,
        "geo_unit_grp" = 'spatial_grp',
        "factor" = collapse_to
      )
    }
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # reset the order
  date_col <- column_mapping$date
  setorderv(
    xgrid_comb,
    'match_strata'
  )

  # set the class as an exposure
  class(xgrid_comb) <- c(class(xgrid_comb), "outcome")

  # set attributes
  attr(xgrid_comb, "column_mapping") <- column_mapping

  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(xgrid_comb))) {
    stop("some NAs persist, investigate and submit a Github issue :) ")
  }

  return(xgrid_comb)
}


