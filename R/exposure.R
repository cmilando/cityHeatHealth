#' Function to clean and prepare the exposure data matrix
#'
#' @param data a dataset of exposures
#' @param column_mapping a named list that indicates relevant columns in `data`. for the exposure
#' data table, these need to be one of: c('date', "exposure", 'geo_unit', 'geo_unit_grp')
#' @param months_subset the warm season months for this region, default is Northern Hemisphere's
#' May through September (5 through 9)
#' @param maxgaps the maximum allowable missing exposure data gap, to be passed to zoo::na.approx (default is 5)
#' @param maxlag the number of lags for the exposure variable (default is 5)
#' @param grp_level whether to summarize to the group level or not (default)
#' @param dt_by is it daily data, or weekly or ...
#'
#' @import data.table
#' @importFrom zoo na.approx
#'
#' @returns a data.table of class("exposure")
#' @export
#'
#' @examples
make_exposure_matrix <- function(data,
                                 column_mapping,
                                 months_subset = 5:9,
                                 dt_by = 'day',
                                 maxgap = 5,
                                 maxlag = 5,
                                 grp_level = FALSE) {

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  #
  setDT(data)

  # validation block
  # set some arbitrary limits on these but users could always make a local
  # copy and override if they really want to
  stopifnot(all(months_subset %in% 1:12) &
              length(months_subset) == length(unique(months_subset)))
  stopifnot(length(maxgap) == 1 & maxgap %in% 1:10)
  stopifnot(length(maxlag) == 1 & maxlag %in% 1:10)
  stopifnot(length(grp_level) == 1 & grp_level %in% c(T, F))

  #
  stopifnot(typeof(column_mapping) == 'list')

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

  #
  if(any(is.na(data))) {
    warning("check about any NA, some corrections for this later,
            but only in certain columns")
  }

  # type checks
  stopifnot(
    inherits(data[[column_mapping$date]], "Date"),
    is.numeric(data[[column_mapping$exposure]]),
    is.character(data[[column_mapping$geo_unit]]),
    is.character(data[[column_mapping$geo_unit_grp]])
  )

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

  # overwrite date
  data[, (column_mapping$date) := as.IDate(get(column_mapping$date))]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CREATE XGRID
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # build xgrid
  # no month subset here, do not change from 1:12!
  xgrid <- make_xgrid(data, column_mapping, months_subset = 1:12, dt_by)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' FILL NA VALUES
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # Loop 1: need to fill in NA values
  exposure_col <- column_mapping$exposure
  join_col <- column_mapping$geo_unit
  date_col <- column_mapping$date

  exposure1_l <- split(xgrid, f = xgrid[, get(join_col)])
  stopifnot(length(exposure1_l) > 0)

  exposure2_l <- vector("list", length(exposure1_l))

  for(i in 1:length(exposure1_l)) {

    x <- exposure1_l[[i]]

    # check for NAs
    ev1 <- is.na(x[, get(exposure_col)])

    # give a warning here
    if(all(ev1)) {
      warning(paste0("all exposures for geo_unit ", x[1, get(join_col)],
                     " were blank -- skipping"))
      next
    }

    # fix the edge case where either the first or last one
    j = 1
    cont = F
    while(is.na(x[j, get(exposure_col)])) {
      j = j + 1
      cont = T
    }
    if(cont) {
      fillx <- x[j, get(exposure_col)]
      x[1:(j - 1), (exposure_col):=fillx]
    }

    j = nrow(x)
    cont = F
    while(is.na(x[j, get(exposure_col)])) {
      j = j - 1
      cont = T
    }
    if(cont) {
      fillx <- x[j, get(exposure_col)]
      x[(j + 1):nrow(x), (exposure_col):= fillx]
    }

    # this should have caught the edge cases, cause a failure if not
    if(is.na(x[1, get(exposure_col)]))
      stop(paste0("first value is still NA despite cleaning for geo_unit ",
                  x[1, get(join_col)]))

    if(is.na(x[nrow(x), get(exposure_col)]))
      stop(paste0("last value is still NA despite cleaning for geo_unit ",
                  x[1, get(join_col)]))

    # # has to be arranged by date so the lags make sense
    x <- x[order(x[, get(date_col)]), ]

    # fill in any NAs using zoo::na.approx
    # hmm this fails if years are not continuous
    # because the gap is too large
    # should i let this fail and then it gets caught later?
    x[, (exposure_col) := zoo::na.approx(get(exposure_col),
                                         maxgap = maxgap)]
    # if(any(is.na(x[, get(exposure_col)])))
    #   stop(paste0("zoo::na.approx() was not able to remove all NAs in ",
    #               x[1, get(join_col)]))

    exposure2_l[[i]] <- x
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' GROUP LEVEL SUMMARY?
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # grp_level summary?
  if(grp_level) {

    geo_grp_col <- column_mapping$geo_unit_grp
    unique_grps <- unlist(unique(data[, get(geo_grp_col)]))
    n_grps <- length(unique_grps)
    exposure2 <- do.call(rbind, exposure2_l)
    setDT(exposure2)

    # 3. Aggregate by group
    # Here I'm assuming you want the mean of exposure columns; adjust as needed
    exposure2 <- exposure2[, lapply(.SD, mean, na.rm = TRUE),
                           by = .(get(geo_grp_col), get(date_col)),
                           .SDcols = exposure_col]

    names(exposure2)[1:2] <- c(column_mapping$geo_unit_grp, column_mapping$date)

    # and make the dates ok again
    # warning("make type checks  (e.g., so Date == Date),
    #      for some reason this doesn't work in some cases? but ok in others?")

    exposure2[, (column_mapping$date) := as.IDate(get(column_mapping$date))]

    # then in order to set up lags you need to split again
    exposure2_l <- split(exposure2, f = exposure2[, get(geo_grp_col)])
    length(exposure2_l)

    # and update column_mapping
    column_mapping$geo_unit <- column_mapping$geo_unit_grp
    column_mapping$geo_unit_grp <- 'ALL'

  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' ADD LAGS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  for(i in 1:length(exposure2_l)) {

    x <-  exposure2_l[[i]]
    setDT(x)

    if(nrow(x) > 0) {

      # Now expand to get lags
      for (k in seq_len(maxlag)) {
        lag_name <- paste0("explag", k)
        x[, (lag_name) := shift(get(exposure_col), k)]
      }

      exposure2_l[[i]] <- x
    }

  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # rbind them all together
  exposure2 <- do.call(rbind, exposure2_l)

  # get just the warm season months
  rr <- month(exposure2[, get(column_mapping$date)]) %in% months_subset
  stopifnot(length(rr) > 1)
  exposure_subset <- exposure2[rr, ]

  # reset the order
  join_col <- column_mapping$geo_unit
  date_col <- column_mapping$date
  setorderv(
    exposure_subset,
    c(join_col, date_col)
  )

  # set the class as an exposure
  class(exposure_subset) <- c(class(exposure_subset), "exposure")

  # set attributes
  attr(exposure_subset, "column_mapping") <- column_mapping

  # at the end there shouldn't be any NAs, so give a warning to investigate
  if(any(is.na(exposure_subset))) {
    stop("some NAs persist, investigate (might be zoo::na_approx) and submit a Github issue :)")
  }

  return(exposure_subset)

}
