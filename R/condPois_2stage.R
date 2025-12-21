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
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes_tbl a data.table of outcomes, created by `make_outcome_table`
#' @param verbose whether to print geo_units as they are run
#' @param global_cen a global centering point
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#' @importFrom mixmeta mixmeta
#' @importFrom mixmeta blup
#' @importFrom dlnm onebasis
#' @importFrom dlnm crosspred
#' @returns
#' @export
#'
#' @examples
condPois_2stage <- function(exposure_matrix,
                            outcomes_tbl,
                            global_cen = NULL,
                            argvar = NULL,
                            arglag = NULL,
                            maxlag = NULL,
                            min_n = NULL,
                            strata_min = 0,
                            verbose = 0) {


  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcomes_tbl))

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' IF the outcomes_tbl has a FACTOR, enter a recursive loop
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////
  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {

    factor_col <- attributes(outcomes_tbl)$column_mapping$factor

    unique_fcts <- unlist(unique(outcomes_tbl[, get(factor_col)]))

    fct_outlist <- vector("list", length(unique_fcts))

    for(fct_i in seq_along(fct_outlist)) {

      if(verbose > 0) {
        cat("<",factor_col,":", unique_fcts[fct_i], ">\n")
      }

      rr <- which(outcomes_tbl[, get(factor_col)] == unique_fcts[fct_i])
      subset_outcomes_tbl <- outcomes_tbl[rr, , drop = FALSE]
      attributes(subset_outcomes_tbl)$column_mapping$factor <- NULL

      # re-call the function, but with just one subset, it should work now
      # and you take the first element here to get rid of having
      # another "_"
      fct_outlist[[fct_i]] <- condPois_2stage(exposure_matrix,
                                              subset_outcomes_tbl,
                                              global_cen = global_cen,
                                              argvar = argvar,
                                              arglag = arglag,
                                              maxlag = maxlag,
                                              min_n = min_n,
                                              strata_min = strata_min,
                                              verbose = verbose)[[1]]


    }

    names(fct_outlist) = unique_fcts

    return(fct_outlist)


  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATION
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## Check 1.5 -- every outcome geo_unit should have data
  ## NOTE - this is sligtly different from the check for pois_single
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  exp_geo_units <- unlist(unique(exposure_matrix[, get(exp_geo_unit_col)]))
  out_geo_units <- unlist(unique(outcomes_tbl[, get(out_geo_unit_col)]))

  stopifnot(all(out_geo_units %in% exp_geo_units))

  # subset so its a complete match
  rr <- which(exposure_matrix[, get(exp_geo_unit_col)] %in% out_geo_units)
  exposure_matrix <- exposure_matrix[rr, ,drop = FALSE]

  ## Check 2
  ## probably should make sure that exposure_matrix and outcomes_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  outcome_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  setorderv(
    exposure_matrix,
    c(exp_geo_unit_col, exp_date_col)
  )

  setorderv(
    outcomes_tbl,
    c(outcome_geo_unit_col, outcome_date_col)
  )

  stopifnot(dim(exposure_matrix)[1] == dim(outcomes_tbl)[1])
  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcomes_tbl[, get(outcome_date_col)]))

  # CHECK 4 geo_unit is the same for both"
  stopifnot(all(outcomes_tbl[, get(outcome_geo_unit_col)] %in%
                  exposure_matrix[, get(exp_geo_unit_col)]))

  # CHECK5
  if(!is.null(global_cen)) {
    stopifnot(is.numeric(global_cen))
  }

  if(verbose > 0) {
    cat("** Validation passed **\n")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 1
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("** STAGE 1 **\n")
  }

  outcome_columns <- attributes(outcomes_tbl)$column_mapping

  geo_cols <- c(
    outcome_columns$geo_unit,
    outcome_columns$geo_unit_grp
  )
  unique_geos <- unique(outcomes_tbl[, ..geo_cols])
  n_geos      <- nrow(unique_geos)

  # get things you need
  cp_list     <- vector("list", n_geos);
  vcov_model  <- vector("list", n_geos);
  argvar_list <- vector("list", n_geos);
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

    if(verbose > 1) {
      cat(this_geo, '\t')
    }

    # this cities exposure matrix
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]

    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]

    # run the model once
    # pass all the main arguments forward
    cp_list[[i]] <- condPois_single(exposure_matrix = single_exposure_matrix,
                                    outcomes_tbl = single_outcomes_tbl,
                                    argvar = argvar, arglag = arglag,
                                    maxlag = maxlag, min_n = min_n,
                                    strata_min = strata_min)

    # get coef and vcov
    coef_list[[i]] <- coef(cp_list[[i]]$cr)
    vcov_list[[i]] <- vcov(cp_list[[i]]$cr)

    # get argvar_list
    argvar_list[[i]]   <- cp_list[[i]]$argvar

    # other things for mixmeta scaling
    exp_mean[[i]]   <- cp_list[[i]]$exp_mean
    exp_IQR[[i]]    <- cp_list[[i]]$exp_IQR
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' MIXMETA
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("** MIXMETA **\n")
  }

  # convert to a matrix
  coef_matrix <- do.call(rbind, coef_list)

  # add in the mean and IQR of exposure in each geo_unit
  unique_geos$exp_mean <- do.call(c, exp_mean)
  unique_geos$exp_IQR <- do.call(c, exp_IQR)

  # create the random effect formula
  # ~ 1 | geo_unit_grp / geo_unit
  # TODO will this work if geo_unit_grp is 'ALL' ?
  #      easy to create a switch if not
  rf = as.formula(paste0("~ 1 | ", outcome_columns$geo_unit_grp, " / ",
                         outcome_columns$geo_unit))

  # see function description for references for this
  meta_fit <- mixmeta(coef_matrix ~ exp_mean + exp_IQR,
                      random = rf,
                      S = vcov_list,
                      data = unique_geos,
                      control = list(showiter=(verbose > 1)),
                      na.action = 'na.exclude')

  # get BLUP
  # there is a `level` argument here

  blup_geo <- blup(meta_fit, vcov = T)
  names(blup_geo) = unique_geos

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 2 - BLUP
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("** STAGE 2 - BLUP **\n")
  }

  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure

  blup_cp   <- vector("list", n_geos);

  # loop through geos
  for(i in seq_along(cp_list)) {

    #
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]

    # get x
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    # get argvar
    this_argvar <- argvar_list[[i]]

    # (1) get onebasis
    basis_x <- do.call("onebasis", modifyList(this_argvar,
                                              list(x=this_exp,
                                                   Boundary.knots = x_b)))

    # *********
    # (2) Center basis
    # either MMT or GLOBAL CEN
    # you need boundary knots because the centerpoint is almost
    # always outside of the percentiles in this work
    # so this creates a full range to test over
    if(!is.null(global_cen)) {
      cen = global_cen
      stopifnot(global_cen >= x_b[1] & global_cen <= x_b[2])
      basis_mmt <- do.call("onebasis", modifyList(this_argvar,
                                                  list(x=global_cen,
                                                       Boundary.knots = x_b)))
    } else {
      cen = cp_list[[i]]$cr$cen
      basis_mmt <- do.call("onebasis",
                           modifyList(this_argvar,
                                      list(x=cen, Boundary.knots = x_b)))
    }

    # *********

    # (3) Center and scale
    basis_cen <- scale(basis_x, center = basis_mmt, scale = FALSE)

    # Apply to a grid
    grid <- seq(from =  x_b[1], to = x_b[2], by = 0.1)

    # get the cross-pred object
    # cen is passed forward from before
    blup_cp[[i]] <- crosspred(basis_cen,
                           cen = cen,
                           coef = blup_geo[[i]]$blup,
                           vcov = blup_geo[[i]]$vcov,
                           model.link = "log",
                           at = grid)

    #
    attr(blup_cp[[i]], "geo_unit") = this_geo
  }

  # set names
  names(blup_cp) <- unique_geos

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # yes, a list in a list, this will make sense later
  # aka, in the recursive call to this function that happens above
  outlist = list(list(meta_fit = meta_fit, blup_cp = blup_cp))
  names(outlist) = "_"
  class(outlist) = 'condPois_2stage'

  return(outlist)

}

