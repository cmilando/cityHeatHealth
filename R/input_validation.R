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

  ## Check 1.5 -- there should be just one geo_unit
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  # added sorting
  exp_geo_unit <- sort(unlist(unique(exposure_matrix[, get(exp_geo_unit_col)])))
  out_geo_unit <- sort(unlist(unique(outcomes_tbl[, get(out_geo_unit_col)])))
  stopifnot(identical(exp_geo_unit, out_geo_unit))

  ## Check 2
  ## probably should make sure that exposure_matrix and outcomes_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  # subset so its a complete match based on DATE and GEO_UNIT
  orig_exp_mapping <- attributes(exposure_matrix)$column_mapping
  exposure_matrix <- exposure_matrix[
    outcomes_tbl,
    on = setNames(
      c(outcome_date_col, out_geo_unit_col),
      c(exp_date_col,    exp_geo_unit_col)
    ),
    nomatch = 0L, drop = F
  ]
  attributes(exposure_matrix)$column_mapping <- orig_exp_mapping

  # get the order correct
  setorderv(
    exposure_matrix,
    c(exp_geo_unit_col, exp_date_col)
  )

  setorderv(
    outcomes_tbl,
    c(out_geo_unit_col, outcome_date_col)
  )

  # check that it worked
  stopifnot(dim(exposure_matrix)[1] == dim(outcomes_tbl)[1])

  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcomes_tbl[, get(outcome_date_col)]))

  stopifnot(identical(exposure_matrix[, get(exp_geo_unit_col)],
                      outcomes_tbl[, get(out_geo_unit_col)]))

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
