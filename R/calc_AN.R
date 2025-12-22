#' Calculate attributable number and attributable rates
#'
#' @param model a model object of class `condPois_2stage`, `condPois_single`, or `condPois_bayes`
#' @param outcomes_tbl a table of outcomes, of class `outcomes`
#' @param pop_data population data
#' @param agg_type what is the spatial resolution you are aggregating to
#' @param join_cols how should you join population data to the outcome table
#' @param nsim number of simulations required for calculation of empirical CI (default = 300)
#' @param verbose 0 = no printing, 1 = headers, 2 = detailed
#' @import data.table
#' @returns
#' @export
#'
#' @examples
calc_AN <- function(model, outcomes_tbl, pop_data,
                    agg_type, join_cols, nsim = 300, verbose = 0) {


  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("outcome" %in% class(outcomes_tbl))

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' RECURSIVE CALL IF FACTOR
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

      # model subset
      sub_model <- model[[unique_fcts[fct_i]]]

      # outcomes subset
      rr <- which(outcomes_tbl[, get(factor_col)] == unique_fcts[fct_i])
      sub_outcomes_tbl <- outcomes_tbl[rr, , drop = FALSE]
      attributes(sub_outcomes_tbl)$column_mapping$factor <- NULL

      # pop_data subset
      rr <- which(pop_data[, get(factor_col)] == unique_fcts[fct_i])
      sub_pop_data <- pop_data[rr, , drop = FALSE]


      # re-call the function, but with just one subset, it should work now
      # and you take the first element here to get rid of having
      # another "_"
      fct_outlist[[fct_i]] <- calc_AN(sub_model,        # <<< **Updated
                                      sub_outcomes_tbl, # <<< **Updated
                                      sub_pop_data,     # <<< **Updated
                                      agg_type,
                                      join_cols,
                                      nsim = 300,
                                      verbose = 0)

      fct_outlist[[fct_i]]$factor_col <- factor_col
      fct_outlist[[fct_i]]$factor_val <- unique_fcts[fct_i]


    }

    names(fct_outlist) = unique_fcts

    class(fct_outlist) = 'calcAN'

    return(fct_outlist)


  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # check agg_type
  stopifnot(agg_type %in% c(geo_unit_col, geo_unit_grp_col, 'all'))

  # check jointype
  stopifnot('population' %in% names(pop_data))
  stopifnot(all(join_cols %in% names(pop_data)))
  setDT(pop_data)

  if(verbose > 0) {
    cat("-- validation passed\n")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' Part 1
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- estimate in each geo_unit\n")
  }

  # starting here
  x <- model$`_`

  # get the blup object
  n_geo_units <- length(x$blup_out)
  AN <- vector("list", n_geo_units)

  #
  outcomes_col     <- attributes(outcomes_tbl)$column_mapping$outcome
  exposure_col     <- x$blup_out[[1]]$exposure_col
  date_col         <- attributes(outcomes_tbl)$column_mapping$date
  geo_unit_col     <- attributes(outcomes_tbl)$column_mapping$geo_unit
  geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp

  for(i in 1:n_geo_units) {

    if(verbose > 1) {
      if(i %% 5 == 0) cat(i, '\t')
    }

    # things you need
    this_geo  <- x$blup_out[[i]]$geo_unit
    basis_cen <- x$blup_out[[i]]$basis_cen
    xcoef     <- x$blup_out[[i]]$coef
    xvcov     <- x$blup_out[[i]]$vcov
    outcomes  <- x$blup_out[[i]]$outcomes
    this_exp  <- x$blup_out[[i]]$this_exp
    cen       <- x$blup_out[[i]]$cen

    # and the outcome database, which should match
    rr <- which(outcomes_tbl[, get(geo_unit_col)] == this_geo)
    single_outcomes_tbl <- outcomes_tbl[rr, ,drop = FALSE]
    stopifnot(identical(single_outcomes_tbl[, get(outcomes_col)], outcomes))

    # convert to matrix
    bvar_mat <- as.matrix(basis_cen)

    # check that AF is always positive
    af_updated <- (1-exp(-bvar_mat %*% xcoef))
    if(any(af_updated < -0.001)) {
      print(summary(af_updated))
      stop(paste0("Attributable Fraction (AF) < -0.001 for ", this_geo,
                  ", which means that centering was likely not done
                  correctly in an earlier step"))
    }

    # create the simulated coefs
    coefsim <- MASS::mvrnorm(nsim,xcoef,xvcov)

    AN_SIM <- vector("list", nsim)

    for(s in seq(nsim)) {

      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE OUTCOMES
      AN_SIM[[s]] <- (1-exp(-bvar_mat %*% coefsim[s,])) * outcomes
    }

    AN_SIM_mat <- data.table(do.call(cbind, AN_SIM))

    colnames(AN_SIM_mat) <- paste0("X", 1:nsim)

    AN_SIM_mat$this_exp = this_exp
    AN_SIM_mat$cen      = cen

    AN_SIM_mat[, (colnames(single_outcomes_tbl)) := single_outcomes_tbl]

    AN_SIM_mat[, year := year(AN_SIM_mat[, get(date_col)])]

    AN[[i]] <- AN_SIM_mat

  }

  if(verbose > 1) {
    cat('\n')
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' Part 2
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- summarize by simulation\n")
  }

  # collapse population data
  pop_data_collapse <- pop_data[, .(
    population = sum(population)),
    by = join_cols
  ]

  AN_ANNUAL <- vector("list", nsim)

  for(xi in 1:nsim) {

    if(verbose > 1) {
      if(xi %% 5 == 0) cat(xi, '\t')
    }

    this_col <- paste0("X", xi)

    AN_sub <- vector("list", n_geo_units)

    get_cols <- c(geo_unit_col, geo_unit_grp_col,
                  'this_exp', 'year', 'cen', this_col)

    for(ani in 1:n_geo_units) {

      setDT(AN[[ani]])
      x1 <- AN[[ani]][, ..get_cols]
      x1$nsim <- xi
      AN_sub[[ani]] <- x1
    }

    AN_sub_all <- do.call(rbind, AN_sub)
    stopifnot(nrow(AN_sub_all) == nrow(outcomes_tbl))
    names(AN_sub_all)[length(get_cols)] <- "attributable_number"

    ## get just those above the centering point
    rr <- which(AN_sub_all$this_exp > AN_sub_all$cen)
    AN_sub_all <- AN_sub_all[rr, ,drop = FALSE]

    ## do an initial fine grain summary
    group_cols = c(
      geo_unit_col, geo_unit_grp_col, 'year', 'nsim'
    )

    AN_ANNUAL[[xi]] <- AN_sub_all[, .(
      annual_AN = round(sum(attributable_number))
    ), by = group_cols]

    n_rows_orig <- nrow( AN_ANNUAL[[xi]] )

    ## and join pop
    stopifnot(all(join_cols %in% names(pop_data_collapse)))
    AN_ANNUAL[[xi]] <- pop_data_collapse[
      AN_ANNUAL[[xi]],
      on = setNames(join_cols, join_cols)
    ]

    if(!(nrow(AN_ANNUAL[[xi]]) == n_rows_orig)) {
      stop("rows after merge is not the same, some error is happening")
    }


    # and then summarize to the annual level per group based on agg_type
    # this is OK to do within an `xi` because you know its all the data
    if(agg_type == geo_unit_col) {

      invisible(1)

    } else if(agg_type == geo_unit_grp_col) {

      group_cols = c(
        geo_unit_col, 'year', 'nsim'
      )

      AN_ANNUAL[[xi]] = AN_sub_all[,.(
        annual_AN = sum(annual_AN),
        population = sum(population)
      ), by = group_cols]

    } else {

      group_cols = c(
        'year', 'nsim'
      )

      AN_ANNUAL[[xi]] = AN_sub_all[,.(
        annual_AN = sum(annual_AN),
        population = sum(population)
      ), by = group_cols]

    }

  }

  AN_ANNUAL <- do.call(rbind, AN_ANNUAL)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  c1 <- which(! (names(AN_ANNUAL) %in% c('population','year', 'nsim','annual_AN') ))
  g1_cols <- c(names(AN_ANNUAL)[c1], "population", "nsim")

  # mean annual
  x1 <- AN_ANNUAL[,.(
    mean_annual_AN = mean(annual_AN)
  ), by = g1_cols]

  # rate column
  x1[, mean_annual_AN_rate := mean_annual_AN / population * 1e5]

  # rate eCI
  g2_cols <- names(AN_ANNUAL)[c1]
  rate_table <- x1[,.(
    mean_annual_attr_rate_est = quantile(mean_annual_AN_rate, 0.50),
    mean_annual_attr_rate_lb = quantile(mean_annual_AN_rate, 0.025),
    mean_annual_attr_rate_ub = quantile(mean_annual_AN_rate, 0.975)
  ), by = g2_cols]

  # number eCI
  number_table <- x1[,.(
    mean_annual_attr_num_est = quantile(mean_annual_AN, 0.50),
    mean_annual_attr_num_lb = quantile(mean_annual_AN, 0.025),
    mean_annual_attr_num_ub = quantile(mean_annual_AN, 0.975)
  ), by = g2_cols]

  outlist <- list(list(rate_table = rate_table, number_table = number_table))
  names(outlist) = "_"
  class(outlist) <- 'calcAN'

  return(outlist)

}
