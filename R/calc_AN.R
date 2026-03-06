#' Calculate attributable number and attributable rates
#'
#' @param model a model object of class `condPois_2stage`, `condPois_1stage`, or `condPois_sb`
#' @param outcomes_tbl a table of outcomes, of class `outcomes`
#' @param pop_data population data
#' @param spatial_agg_type what is the spatial resolution you are aggregating to
#' @param spatial_join_col how should you join population data to the outcome table
#' @param nsim number of simulations required for calculation of empirical CI (default = 300)
#' @param verbose 0 = no printing, 1 = headers, 2 = detailed
#' @param by_year Export annual counts
#' @import data.table
#' @importFrom tidyr expand_grid
#' @returns
#' @export
#'
#' @examples
calc_AN <- function(model, outcomes_tbl, pop_data,
                    spatial_agg_type, spatial_join_col,
                    by_year = FALSE,
                    nsim = 300, verbose = 0) {


  ## Check 1 -- that both inputs are the right class of variables
  stopifnot(class(model) %in%
              c('condPois_1stage', 'condPois_1stage_list',
                'condPois_2stage', 'condPois_2stage_list',
                'condPois_sb',  'condPois_sb_list'))

  stopifnot("outcome" %in% class(outcomes_tbl))

  stopifnot('population' %in% names(pop_data))
  stopifnot(all(spatial_join_col %in% names(pop_data)))

  stopifnot(is.data.table(pop_data))

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' RECURSIVE CALL IF FACTOR
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) &
     grepl("_list", class(model))) {

    factor_col <- attributes(outcomes_tbl)$column_mapping$factor

    unique_fcts <- unlist(unique(outcomes_tbl[, get(factor_col)]))

    fct_outlist <- vector("list", length(unique_fcts))

    for(fct_i in seq_along(fct_outlist)) {

      if(verbose > 0) {
        cat("<",factor_col,":", unique_fcts[fct_i], ">\n")
      }

      # model subset
      sub_model <- model[[ unique_fcts[fct_i] ]]
      stopifnot(class(sub_model) %in% c('condPois_1stage', 'condPois_2stage',
                                        'condPois_sb'))
      stopifnot("_" %in% names(sub_model))

      # outcomes subset
      stopifnot(factor_col %in% names(outcomes_tbl))
      rr <- which(outcomes_tbl[, get(factor_col)] == unique_fcts[fct_i])
      sub_outcomes_tbl <- outcomes_tbl[rr, , drop = FALSE]
      attributes(sub_outcomes_tbl)$column_mapping$factor <- NULL

      # pop_data subset
      stopifnot(factor_col %in% names(pop_data))
      rr <- which(pop_data[, get(factor_col)] == unique_fcts[fct_i])
      sub_pop_data <- pop_data[rr, , drop = FALSE]

      # re-call the function, but with just one subset
      fct_outlist[[fct_i]] <- calc_AN(sub_model,        # <<< **Updated
                                      sub_outcomes_tbl, # <<< **Updated
                                      sub_pop_data,     # <<< **Updated
                                      spatial_agg_type,
                                      spatial_join_col,
                                      nsim = nsim,
                                      verbose = verbose)

      fct_outlist[[fct_i]]$factor_col <- factor_col
      fct_outlist[[fct_i]]$factor_val <- unique_fcts[fct_i]


    }

    names(fct_outlist) = unique_fcts

    class(fct_outlist) = 'calcAN_list'

    return(fct_outlist)


  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  outcomes_col     <- attributes(outcomes_tbl)$column_mapping$outcome
  date_col         <- attributes(outcomes_tbl)$column_mapping$date
  geo_unit_col     <- attributes(outcomes_tbl)$column_mapping$geo_unit
  geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp

  # check spatial_agg_type
  stopifnot(spatial_agg_type %in% c(geo_unit_col, geo_unit_grp_col, 'all'))

  # TODO: check jointype
  stopifnot(length(spatial_join_col) == 1)
  stopifnot(length(spatial_agg_type) == 1)

  # check pop data
  setDT(pop_data)

  stopifnot(all(spatial_join_col %in% names(pop_data)))

  pop_data_collapse <- pop_data[, .(
    population = sum(population)),
    by = spatial_join_col
  ]

  if(any(pop_data_collapse$population == 0)) {
    warning("some pop data are zero")
    pop_data_collapse <- subset(pop_data_collapse, population > 0)
  }

  # subsetting
  safe_geos <- unique(pop_data_collapse[, ..spatial_join_col])
  outcomes_tbl <- safe_geos[outcomes_tbl, on = spatial_join_col, nomatch = 0]


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

  ## >>>>>>>>>>>
  ## TODO
  ## DO SOME EXTRA SUBSETTING because not every x <- model$_ has
  ## every geo_unit I think, so the purpose of this is to remove
  ## the difference between pop_data and model
  ## >>>>>>>>>>>


  # get the blup object
  n_geo_units <- length(x$out)
  if(!(n_geo_units >= 1)) {
    stop("the model object is empty, check that model and outcomes are the same type
         (e.g., that one is not a `_list` type while the other is a standard.")
  }
  stopifnot(n_geo_units >= 1)

  AN <- vector("list", n_geo_units)

  #
  exposure_col     <- x$out[[1]]$exposure_col

  for(i in 1:n_geo_units) {

    if(verbose > 1) {
      if(i %% 5 == 0) cat(i, '\t')
    }

    # things you need
    this_geo     <- x$out[[i]]$geo_unit
    basis_cen    <- x$out[[i]]$basis_cen
    xcoef        <- x$out[[i]]$coef
    xvcov        <- x$out[[i]]$vcov
    outcomes     <- x$out[[i]]$outcomes
    match_strata <- x$out[[i]]$match_strata
    this_exp     <- x$out[[i]]$this_exp
    cen          <- x$out[[i]]$cen
    global_cen   <- x$out[[i]]$global_cen

    # stop if this isn't in pop_data
    if(! (this_geo %in% as.vector(unlist(safe_geos)))) {
      cat("skipping based on no population data: ", this_geo, "\t")
      next
    }

    # and the outcome database, which should match
    rr <- which(outcomes_tbl$match_strata %in% match_strata)
    single_outcomes_tbl <- outcomes_tbl[rr, ,drop = FALSE]
    date_fmt <- single_outcomes_tbl[, get(date_col)]

    if(!identical(single_outcomes_tbl[, get(outcomes_col)], outcomes)) {
      cat("This geo:\n")
      cat(this_geo, '\n')
      cat("Outcomes table being passed in:\n")
      print(head(single_outcomes_tbl))
      cat("Outcomes vector associated with the model object for this geo:\n")
      print(head(outcomes))
      stop("Outcomes not the same - usually this means there is a spatial
           mismatch somewhere in one of the other datasets - in this case
           either outcomes_tbl or the model doesn't have the spatial resolution
           that you are asking for. Check that each one has the spatial units
           and each level that are required for matching.")
    }

    # convert to matrix
    bvar_mat <- as.matrix(basis_cen)

    # check that AF is always positive
    af_updated <- (1-exp(-bvar_mat %*% xcoef))
    if(any(af_updated < -0.001) & is.null(global_cen)) {
      print(summary(af_updated))
      stop(paste0("Attributable Fraction (AF) < -0.001 for ", this_geo,
                  "and global_cen is NULL,
                  which means that centering was likely not done
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
    AN_SIM_mat$date_fmt = date_fmt

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

  ## you have missing rows here because not every
  ## place has values below MMT in every year
  ## so its almost like you need another xgrid
  ## by geo_unit_col, geo_unit_grp_col, 'year', and above_MMT
  geo_cols <- c(
    geo_unit_col,
    geo_unit_grp_col
  )
  unique_geos <- unique(outcomes_tbl[, ..geo_cols])
  unique_years <- unique(year(outcomes_tbl[, get(date_col)]))

  ## now expand_grid for every year and above_MMT = c(T, F)
  ## actually, removing nsim here
  xgrid <- tidyr::expand_grid(data.frame(unique_geos),
                              year = unique_years,
                              above_MMT = c(T, F))
  setDT(xgrid)

  ## join w population

  AN_ANNUAL <- vector("list", nsim)

  for(xi in 1:nsim) {

    if(verbose > 1) {
      if(xi %% 5 == 0) cat(xi, '\t')
    }

    this_col <- paste0("X", xi)

    AN_sub <- vector("list", n_geo_units)

    get_cols <- c(geo_unit_col, geo_unit_grp_col,
                  'this_exp', 'year', 'cen', 'date_fmt', this_col)

    for(ani in 1:n_geo_units) {

      if(!is.null(AN[[ani]])) {

        setDT(AN[[ani]])
        x1 <- AN[[ani]][, ..get_cols]
        x1$nsim <- xi
        AN_sub[[ani]] <- x1

      }
    }

    AN_sub_all <- do.call(rbind, AN_sub)

    # rename
    names(AN_sub_all)[length(get_cols)] <- "attributable_number"

    ## get just those above the centering point
    # rr <- AN_sub_all$this_exp > AN_sub_all$cen)
    # AN_sub_all <- AN_sub_all[rr, ,drop = FALSE]
    AN_sub_all$above_MMT = AN_sub_all$this_exp >= AN_sub_all$cen

    ## do an initial fine grain summary
    group_cols = c(
      geo_unit_col, geo_unit_grp_col, 'year', 'above_MMT'
    )

    AN_ANNUAL[[xi]] <- AN_sub_all[, .(
      annual_AN = round(sum(attributable_number))
    ), by = group_cols]

    sum_AN_orig <- sum(AN_ANNUAL[[xi]]$annual_AN)

    ## merge to xgrid
    ## which way does this go?
    AN_ANNUAL[[xi]] <- AN_ANNUAL[[xi]][
      xgrid,
      on = setNames(group_cols, group_cols)
    ]

    # mannual set xi
    AN_ANNUAL[[xi]]$nsim = xi

    # fill in NAs
    rr <- which(is.na(AN_ANNUAL[[xi]]$annual_AN))
    AN_ANNUAL[[xi]]$annual_AN[rr] <- 0
    if(any(is.na(AN_ANNUAL[[xi]]))) stop("AN expand_grid didn't work correctly")

    ## and join pop
    AN_ANNUAL[[xi]] <- pop_data_collapse[
      AN_ANNUAL[[xi]],
      on = setNames(spatial_join_col, spatial_join_col)
    ]

    sum_AN_new <- sum(AN_ANNUAL[[xi]]$annual_AN)

    if(!(sum_AN_new == sum_AN_orig)) {
      stop("sum_AN after merge is not the same, some error is happening")
    }


    # and then summarize to the annual level per group based on spatial_agg_type
    # this is OK to do within an `xi` because you know its all the data
    if(spatial_agg_type == geo_unit_col) {

      invisible(1)

    } else if(spatial_agg_type == geo_unit_grp_col) {

      group_cols = c(
        geo_unit_grp_col, 'year', 'nsim', 'above_MMT'
      )

      AN_ANNUAL[[xi]] = AN_ANNUAL[[xi]][,.(
        annual_AN = sum(annual_AN, na.rm = T),
        population = sum(population, na.rm = T)
      ), by = group_cols]

    } else {

      AN_ANNUAL[[xi]]$all <- 'ALL'

      group_cols = c(
        'all', 'year', 'nsim', 'above_MMT'
      )

      AN_ANNUAL[[xi]] = AN_ANNUAL[[xi]][,.(
        annual_AN = sum(annual_AN, na.rm = T),
        population = sum(population, na.rm = T)
      ), by = group_cols]

    }

  }

  AN_ANNUAL <- do.call(rbind, AN_ANNUAL)

  if(verbose > 1) {
    cat('\n')
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  c1 <- which(! (names(AN_ANNUAL) %in%
                   c('population', 'year', 'nsim','annual_AN', 'above_MMT') ))

  # remove year and annual_AN
  g1_cols <- c(names(AN_ANNUAL)[c1], "population", "nsim", 'above_MMT')

  # mean annual
  if(by_year == TRUE) {
    x1 <- AN_ANNUAL[,.(
      mean_annual_AN = mean(annual_AN)
    ), by = c(g1_cols, 'year')]
  } else {
    x1 <- AN_ANNUAL[,.(
      mean_annual_AN = mean(annual_AN)
    ), by = g1_cols]
  }


  # final checks
  if(any(is.na(x1$mean_annual_AN))) {
    rr <- which(is.na(x1$mean_annual_AN))
    print(x1[rr, ])
    warning('annual_AN has NA, find out why')
    return(x1)
  }

  if(any(is.na(x1$population))) {
    rr <- which(is.na(x1$population))
    print(x1[rr, ])
    warning('population has NA, find out why')
    return(x1)
  }

  # rate column
  x1[, mean_annual_AN_rate := mean_annual_AN / population * 1e5]

  if(any(is.na(x1$mean_annual_AN_rate))) {
    rr <- which(is.na(x1$mean_annual_AN_rate))
    print(x1[rr, ])
    warning('mean_annual_AN_rate has NA, find out why')
    return(x1)
  }

  # rate eCI
  if(by_year == TRUE) {
    g2_cols <- c(names(AN_ANNUAL)[c1], "year", "population", 'above_MMT')
  } else {
    g2_cols <- c(names(AN_ANNUAL)[c1], "population", 'above_MMT')
  }

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

  outlist <- list(list(
    rate_table = rate_table,
    number_table = number_table))

  names(outlist) = "_"
  class(outlist) <- 'calcAN'

  return(outlist)

}

#' @export
#' print.calcAN
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.calcAN <- function(x) {
  cat("< an object of class `calcAN` >\n")
  invisible(x)
}

#' @export
#' print.calcAN_list
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.calcAN_list <- function(x) {
  cat("< an object of class `calcAN_list`:",
      paste(names(x), collapse = ",")," >\n")
  invisible(x)
}

#' @export
#' plot.calcAN
#'
#' @param x an object of class plot.calcAN
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @param spatial_sub an option argument to subset the geo_units investigated
#' @param override_limit override the built-in plot limit
#' @importFrom ggplot2 ggplot
#' @importFrom scales number
#' @import data.table
#' @returns
#' @export
#'
#' @examples
plot.calcAN <- function(x, table_type, above_MMT, spatial_sub = NULL,
                        override_limit = FALSE) {

  stopifnot(table_type %in% c('rate', 'num'))
  stopifnot(above_MMT %in% c(T, F))
  geo_unit_col <- attributes(x$`_`)$column_mapping
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  if(table_type == 'num') {

    byX_df <- x$`_`$number_table
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    if(length(spatial_sub) > 0) {
      rr <- which(byX_df[[x_col]] %in% spatial_sub)
      byX_df <- byX_df[rr, ]
    }

    if(nrow(byX_df) > 20 & !override_limit) {
      warning("plot elements > 20, subsetting to top 20")
      setorder(byX_df, -mean_annual_attr_num_est)
      byX_df <- byX_df[1:20]
    }

    ggplot(byX_df) +
      coord_cartesian() +
      theme_classic() +
      geom_col(aes(
        y = mean_annual_attr_num_est,
        x = reorder(!!sym(x_col), mean_annual_attr_num_est)),
        linewidth = 0.4, fill = 'lightblue', col = 'black') +
      geom_errorbar(aes(
        ymax = mean_annual_attr_num_ub,
        ymin = mean_annual_attr_num_lb,
        x = reorder(!!sym(x_col), mean_annual_attr_num_est)),
        color = 'grey15', width = 0.3,
        linewidth = 0.4) +
      geom_label(aes(x = reorder(!!sym(x_col), mean_annual_attr_num_est),
                     y = mean_annual_attr_num_est,
                     label = number(round(mean_annual_attr_num_est, 0),
                                    big.mark = ',')),
                 color = 'grey15', linewidth = 0.4, fill = 'white',
                 size = 2) +
      theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust = 0.5),
            axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      xlab(x_col) +
      ylab(paste0("Annual Temp Attr.\nOutcomes\n", ylab_flag, "\n(#)"))

  } else {

    byX_df <- x$`_`$rate_table
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    if(length(spatial_sub) > 0) {
      rr <- which(byX_df[[x_col]] %in% spatial_sub)
      byX_df <- byX_df[rr, ]
    }

    if(nrow(byX_df) > 20 & !override_limit) {
      warning("plot elements > 20, subsetting to top 20")
      setorder(byX_df, -mean_annual_attr_rate_est)
      byX_df <- byX_df[1:20]
    }

    ggplot(byX_df) +
      coord_cartesian() +
      theme_classic() +
      geom_col(aes(
        y = mean_annual_attr_rate_est,
        x = reorder(!!sym(x_col), mean_annual_attr_rate_est)),
        linewidth = 0.4, fill = 'lightblue', col = 'black') +
      geom_errorbar(aes(
        ymax = mean_annual_attr_rate_ub,
        ymin = mean_annual_attr_rate_lb,
        x = reorder(!!sym(x_col), mean_annual_attr_rate_est)),
        color = 'grey15', width = 0.3,
        linewidth = 0.4) +
      geom_label(aes(x = reorder(!!sym(x_col), mean_annual_attr_rate_est),
                     y = mean_annual_attr_rate_est,
                     label = number(round(mean_annual_attr_rate_est, 1),
                                    big.mark = ',')),
                 color = 'grey15', linewidth = 0.4, fill = 'white',
                 size = 2) +
      theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust = 0.5),
            axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      xlab(x_col) +
      ylab(paste0("Annual Temp Attr.\nOutcomes Rate\n",
                  ylab_flag, "\n(# per 100,000)"))

  }

}


#' @export
#' plot.calcAN_list
#'
#' @param x an object of class plot.calcAN_list
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @param spatial_sub an option argument to subset the geo_units investigated
#' @param ordered_levels factor levels
#' @param override_limit override the built-in plot limit
#' @importFrom ggplot2 ggplot
#' @import data.table
#' @returns
#' @export
#'
#' @examples
plot.calcAN_list <- function(x, table_type, above_MMT, spatial_sub = NULL,
                             ordered_levels = NULL, override_limit = FALSE) {

  stopifnot(table_type %in% c('rate', 'num'))
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  fct_names <- names(x)
  fct_lab <- x[[1]]$factor_col

  byX_df  <- vector("list", length(fct_names))

  if(table_type == 'num') {

    for(i in 1:length(byX_df)) {
      byX_df[[i]] <- x[[fct_names[i]]]$`_`$number_table
      byX_df[[i]][, (fct_lab) := fct_names[i]]
      x_col <- names(byX_df[[i]])[1]

      if(length(spatial_sub) > 0) {
        rr <- which(byX_df[[i]][[x_col]] %in% spatial_sub)
        byX_df[[i]] <- byX_df[[i]][rr, ]
      }

      if(nrow(byX_df[[i]]) > 20 & !override_limit) {
        warning("plot elements > 20, subsetting to top 20")
        setorder(byX_df[[i]], -mean_annual_attr_num_est)
        byX_df[[i]] <- byX_df[[i]][1:20]
      }
    }

    byX_df <- do.call(rbind, byX_df)
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    ##
    pd <- position_dodge2(width = 0.9, preserve = "single")

    ##
    if(!is.null(ordered_levels)) {
      byX_df[[fct_lab]] = factor(byX_df[[fct_lab]],
                                 levels  = ordered_levels,
                                 ordered = T)
    }


    ggplot(byX_df) +
      coord_cartesian() +
      theme_classic() +
      geom_col(aes(
        y = mean_annual_attr_num_est,
        fill = !!sym(fct_lab),
        group = !!sym(fct_lab),
        x = reorder(!!sym(x_col), mean_annual_attr_num_est)),
        linewidth = 0.4, col = 'black', position = pd) +
      geom_errorbar(aes(
        group = !!sym(fct_lab),
        ymax = mean_annual_attr_num_ub,
        ymin = mean_annual_attr_num_lb,
        x = reorder(!!sym(x_col), mean_annual_attr_num_est)),
        color = 'grey15', position = pd,
        linewidth = 0.4) +
      geom_label(aes(x = reorder(!!sym(x_col), mean_annual_attr_num_est),
                     y = mean_annual_attr_num_est,
                     group = !!sym(fct_lab),
                     label = number(round(mean_annual_attr_num_est, 0),
                                    big.mark = ',')),
                 color = 'grey15', linewidth = 0.4, fill = 'white',
                 size = 2, position = pd) +
      theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust = 0.5),
            axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      xlab(x_col) +
      ylab(paste0("Annual Temp Attr.\nOutcomes\n", ylab_flag, "\n(#)"))

  } else {

    for(i in 1:length(byX_df)) {
      byX_df[[i]] <- x[[fct_names[i]]]$`_`$rate_table
      byX_df[[i]][, (fct_lab) := fct_names[i]]

      x_col <- names(byX_df[[i]])[1]

      if(length(spatial_sub) > 0) {
        rr <- which(byX_df[[i]][[x_col]] %in% spatial_sub)
        byX_df[[i]] <- byX_df[[i]][rr, ]
      }

      if(nrow(byX_df[[i]]) > 20 & !override_limit) {
        warning("plot elements > 20, subsetting to top 20")
        setorder(byX_df[[i]], -mean_annual_attr_rate_est)
        byX_df[[i]] <- byX_df[[i]][1:20]
      }
    }

    byX_df <- do.call(rbind, byX_df)
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    ##
    pd <- position_dodge2(width = 0.9, preserve = "single")

    ##
    if(!is.null(ordered_levels)) {
      byX_df[[fct_lab]] = factor(byX_df[[fct_lab]],
                                 levels  = ordered_levels,
                                 ordered = T)
    }

    ggplot(byX_df) +
      coord_cartesian() +
      theme_classic() +
      geom_col(aes(
        y = mean_annual_attr_rate_est,
        fill = !!sym(fct_lab),
        group = !!sym(fct_lab),
        x = reorder(!!sym(x_col), mean_annual_attr_rate_est)),
        linewidth = 0.4, col = 'black', position = pd) +
      geom_errorbar(aes(
        group = !!sym(fct_lab),
        ymax = mean_annual_attr_rate_ub,
        ymin = mean_annual_attr_rate_lb,
        x = reorder(!!sym(x_col), mean_annual_attr_rate_est)),
        color = 'grey15', position = pd,
        linewidth = 0.4) +
      geom_label(aes(x = reorder(!!sym(x_col), mean_annual_attr_rate_est),
                     y = mean_annual_attr_rate_est,
                     group = !!sym(fct_lab),
                     label = number(round(mean_annual_attr_rate_est, 0),
                                    big.mark = ',')),
                 color = 'grey15', linewidth = 0.4, fill = 'white',
                 size = 2, position = pd) +
      theme(axis.text.x = element_text(angle= 90, hjust = 1, vjust = 0.5),
            axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      xlab(x_col) +
      ylab(paste0("Annual Temp Attr.\nOutcomes Rate\n",
                  ylab_flag, "\n(# per 100,000)"))

  }

}

#' @export
#' spatial_plot.calcAN
#'
#' @param x an object of class condPois_AN
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @param pal color palette
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
spatial_plot.calcAN <- function(x, shp, table_type, above_MMT, pal = 'Purples') {

  stopifnot(table_type %in% c('rate', 'num'))
  stopifnot(above_MMT %in% c(T, F))
  geo_unit_col <- attributes(x$`_`)$column_mapping
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  if(table_type == 'num') {

    byX_df <- x$`_`$number_table
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    # join to SF
    stopifnot(geo_unit_col %in% names(shp)) # not a bad first check
    shp_w_data <- merge(shp, byX_df)

    return(ggplot(shp_w_data) +
             theme_classic() +
             geom_sf(aes(fill = mean_annual_attr_num_est)) +
             scale_fill_binned(name = paste0("Annual Temp Attr.\nOutcomes\n",
                                             ylab_flag, "\n(#)"),
                               palette = pal) +
             theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.line = element_blank(),
                   axis.title = element_blank()))

  } else {

    byX_df <- x$`_`$rate_table
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]

    # join to SF
    stopifnot(geo_unit_col %in% names(shp)) # not a bad first check
    shp_w_data <- merge(shp, byX_df)

    return(ggplot(shp_w_data) +
             theme_classic() +
             geom_sf(aes(fill = mean_annual_attr_rate_est)) +
             scale_fill_binned(name = paste0("Annual Temp Attr.\nOutcomes Rate\n",
                                             ylab_flag, "\n(# per 100,000)"),
                               palette = pal) +
             theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.line = element_blank(),
                   axis.title = element_blank()))

  }

}

#' @export
#' spatial_plot.calcAN_list
#'
#' @param x an object of class condPois_AN
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal.info
#' @returns
#' @export
#'
#' @examples
spatial_plot.calcAN_list <- function(x, shp, table_type, above_MMT) {

  stopifnot(table_type %in% c('rate', 'num'))
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  fct_names <- names(x)
  fct_lab <- x[[1]]$factor_col

  plt_obj <- vector("list", length(names(x)))

  set.seed(123)
  tot <- row.names(subset(RColorBrewer::brewer.pal.info,
                          category == 'seq'))

  color_pals <- sample(tot, length(plt_obj), replace = F)

  for(i in 1:length(plt_obj)) {

    plt_obj[[i]] <- spatial_plot(x[[names(x)[i]]], shp,
                                 table_type, above_MMT, color_pals[i]) +
      ggtitle(paste0(fct_lab, ": ", fct_names[i]))

  }

  patchwork::wrap_plots(plt_obj, ncol = 1)

}

