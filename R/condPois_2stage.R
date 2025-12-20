#' Run a conditional poisson model for across many geographic units
#'
#' MixMeta
#' nope its this: https://www.sciencedirect.com/science/article/pii/S0140673614621140?via%3Dihub
#' avgtmax and range and random effect for county
#' unclear, Gasp does both here https://github.com/gasparrini/RiskExtrapolation/blob/3caa5c46a748d6a789021554f435e53f9960d4c1/24_MetaRegression.R
#' and it goes counrty/city
#' https://github.com/gasparrini/MCC-SO2/blob/1a158b31bf83f3bbca3449b553a6570280306299/Rcode/06.nonlin.R#L47
#' nlmeta <- mixmeta(nlcoefall, nlvcovall, random=~1|country/city, data=cities,
#'                   control=list(showiter=T, igls.inititer=10))
#' basically, did Gasp 2015 (temp_mean, temp_iqr, random effect for county)
#' and included the hierarchicacy structure for ids within countries
#' https://github.com/gasparrini/2015_gasparrini_EHP_Rcodedata
#' https://github.com/gasparrini/2019_sera_StatMed_Rcode/blob/master/03.MultilevelMA.R
#'
#' @param exposure_matrix
#' @param outcome_tbl
#' @param verbose whether to print geo_units as they are run
#' @param global_cen a global centering point
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#' @importFrom mixmeta mixmeta
#' @importFrom mixmeta blup
#' @returns
#' @export
#'
#' @examples
condPois_2stage <- function(exposure_matrix,
                            outcome_tbl,
                            global_cen = NULL,
                            verbose = F,
                            argvar = NULL,
                            arglag = NULL,
                            maxlag = NULL,
                            min_n = NULL) {


  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATION
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcome_tbl))

  ## Check 1.5 -- every outcome geo_unit should have data
  ## NOTE - this is sligtly different from the check for pois_single
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcome_tbl)$column_mapping$geo_unit

  exp_geo_units <- unlist(unique(exposure_matrix[, get(exp_geo_unit_col)]))
  out_geo_units <- unlist(unique(outcome_tbl[, get(out_geo_unit_col)]))

  stopifnot(all(out_geo_units %in% exp_geo_units))

  # subset so its a complete match
  rr <- which(exposure_matrix[, get(exp_geo_unit_col)] %in% out_geo_units)
  exposure_matrix <- exposure_matrix[rr, ,drop = FALSE]

  ## Check 2
  ## probably should make sure that exposure_matrix and outcome_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcome_tbl)$column_mapping$date

  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  outcome_geo_unit_col <- attributes(outcome_tbl)$column_mapping$geo_unit

  setorderv(
    exposure_matrix,
    c(exp_geo_unit_col, exp_date_col)
  )

  setorderv(
    outcome_tbl,
    c(outcome_geo_unit_col, outcome_date_col)
  )

  stopifnot(dim(exposure_matrix)[1] == dim(outcome_tbl)[1])
  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcome_tbl[, get(outcome_date_col)]))

  # CHECK 4 geo_unit is the same for both"
  stopifnot(all(outcome_tbl[, get(outcome_geo_unit_col)] %in%
                  exposure_matrix[, get(exp_geo_unit_col)]))

  # CHECK 5
  if("factor" %in% names(attributes(outcome_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 1
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  outcome_columns <- attributes(outcome_tbl)$column_mapping

  geo_cols <- c(
    outcome_columns$geo_unit,
    outcome_columns$geo_unit_grp
  )
  unique_geos <- unique(outcome_tbl[, ..geo_cols])
  n_geos      <- nrow(unique_geos)

  # get things you need
  cp_list     <- vector("list", n_geos)
  vcov_model  <- vector("list", n_geos);
  exp_mean    <- vector("list", n_geos);
  exp_IQR     <- vector("list", n_geos);
  vcov_list   <- vector("list", n_geos);
  coef_list   <- vector("list", n_geos);
  names(vcov_model) <- unique_geos
  names(coef_list) <- unique_geos

  # loop through geos
  for(i in seq_along(cp_list)) {

    # get the name, which you know exists in both datasets
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]

    if(verbose) {
      cat(this_geo, '\t')
    }

    # this cities exposure matrix
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]

    rr <- outcome_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcome_tbl = outcome_tbl[rr, ,drop = FALSE]

    # run the model once
    cp_list[[i]] <- condPois_single(single_exposure_matrix, single_outcome_tbl)

    # get coef and vcov
    coef_list[[i]] <- coef(cp_list[[i]]$cr)
    vcov_list[[i]] <- vcov(cp_list[[i]]$cr)

    # other things for mixmeta scaling
    exp_mean[[i]]   <- cp_list[[i]]$exp_mean
    exp_IQR[[i]]    <- cp_list[[i]]$exp_IQR
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' MIXMETA
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # convert to a matrix
  coef_matrix <- do.call(rbind, coef_list)

  # add in the mean and IQR of exposure in each geo_unit
  unique_geos$exp_mean <- do.call(c, exp_mean)
  unique_geos$exp_IQR <- do.call(c, exp_IQR)

  # create the random effect formula
  rf = as.formula(paste0("~ 1 | ", outcome_columns$geo_unit_grp))

  # see function description for references for this
  meta_fit <- mixmeta(coef_matrix ~ exp_mean + exp_IQR,
                      random = rf,
                      S = vcov_list,
                      data = unique_geos,
                      control = list(showiter=verbose),
                      na.action = 'na.exclude')


  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 2 - BLUP
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # have to think about this
  warning("global_cen")

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  warning("probably want to output a whole model here?
          or a way to access individual results")


}
