#' Input validation
#'
#' @param exposure_matrix
#' @param outcomes_tbl
#'
#' @returns
#' @export
#'
#' @examples
input_validation <- function(exposure_matrix, outcomes_tbl) {

  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcomes_tbl))

  # ## Check 1.5 -- there should be just one geo_unit
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit
  #
  # # added sorting
  # exp_geo_unit <- sort(unlist(unique(exposure_matrix[, get(exp_geo_unit_col)])))
  # out_geo_unit <- sort(unlist(unique(outcomes_tbl[, get(out_geo_unit_col)])))
  # if(!identical(exp_geo_unit, out_geo_unit)) {
  #   print(exp_geo_unit)
  #   print(out_geo_unit)
  #   stop('geo units not identical')
  # }

  ## Check 2
  ## probably should make sure that exposure_matrix and outcomes_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  # subset so its a complete match based on DATE and GEO_UNIT
  # Update this actually to be based on STRATA
  exp_strata = unique(exposure_matrix$match_strata)
  out_strata = unique(outcomes_tbl$match_strata)
  all_strata = intersect(exp_strata, out_strata)

  # exposure subset
  rr <- which(exposure_matrix$match_strata %in% all_strata)
  exposure_matrix <- exposure_matrix[rr, ]

  # and you also have to remove outcomes for which there are no exposures
  ss <- which(outcomes_tbl$match_strata %in% all_strata)
  outcomes_tbl <- outcomes_tbl[ss, ]

  # get the order correct
  setorderv(
    exposure_matrix,
    "match_strata"
  )

  setorderv(
    outcomes_tbl,
    "match_strata"
  )

  # check that it worked
  if(!(dim(exposure_matrix)[1] == dim(outcomes_tbl)[1])) {
    stop("Dim of exposure matrix != dim of outcomes tbl. Usually means
         there is a bug in the match_strata")
  }

  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcomes_tbl[, get(outcome_date_col)]))

  # they should be line by line equal
  stopifnot(identical(exposure_matrix$strata, outcomes_tbl$strata))
  stopifnot(identical(exposure_matrix$match_strata, outcomes_tbl$match_strata))

  if(!identical(exposure_matrix[, get(exp_geo_unit_col)],
                      outcomes_tbl[, get(out_geo_unit_col)])) {

    rr_df <- data.frame(expmat = exposure_matrix[, get(exp_geo_unit_col)],
                        outmat = outcomes_tbl[, get(out_geo_unit_col)])

    rr <- which(rr_df[, 1] != rr_df[, 2])

    print(unique(rr_df[rr, ]))

    stop("geo unit col mismatch")

  }

  # CHECK 4 geo_unit is the same for both"
  stopifnot(all(outcomes_tbl[, get(out_geo_unit_col)] %in%
                  exposure_matrix[, get(exp_geo_unit_col)]))

  # CHECK 4
  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }

  return(list(outcomes_tbl = outcomes_tbl,
              exposure_matrix = exposure_matrix))

}
