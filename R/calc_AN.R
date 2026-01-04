#' Calculate attributable number and attributable rates
#'
#' @param model a model object of class `condPois_2stage`, `condPois_1stage`, or `condPois_sb`
#' @param outcomes_tbl a table of outcomes, of class `outcomes`
#' @param pop_data population data
#' @param agg_type what is the spatial resolution you are aggregating to
#' @param join_cols how should you join population data to the outcome table
#' @param nsim number of simulations required for calculation of empirical CI (default = 300)
#' @param verbose 0 = no printing, 1 = headers, 2 = detailed
#' @import data.table
#' @importFrom tidyr expand_grid
#' @returns
#' @export
#'
#' @examples
calc_AN <- function(model, outcomes_tbl, pop_data,
                    agg_type, join_cols, nsim = 300, verbose = 0) {


  ## Check 1 -- that both inputs are the right class of variables
  stopifnot(class(model) %in%
              c('condPois_1stage', 'condPois_1stage_list',
                'condPois_2stage', 'condPois_2stage_list',
                'condPois_bayes',  'condPois_bayes_list'))

  stopifnot("outcome" %in% class(outcomes_tbl))

  stopifnot('population' %in% names(pop_data))
  stopifnot(all(join_cols %in% names(pop_data)))

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
                                      agg_type,
                                      join_cols,
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

  # check agg_type
  stopifnot(agg_type %in% c(geo_unit_col, geo_unit_grp_col, 'all'))

  # check jointype

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
  n_geo_units <- length(x$out)
  stopifnot(n_geo_units >= 1)
  AN <- vector("list", n_geo_units)

  #
  exposure_col     <- x$out[[1]]$exposure_col

  for(i in 1:n_geo_units) {

    if(verbose > 1) {
      if(i %% 5 == 0) cat(i, '\t')
    }

    # things you need
    this_geo   <- x$out[[i]]$geo_unit
    basis_cen  <- x$out[[i]]$basis_cen
    xcoef      <- x$out[[i]]$coef
    xvcov      <- x$out[[i]]$vcov
    outcomes   <- x$out[[i]]$outcomes
    this_exp   <- x$out[[i]]$this_exp
    cen        <- x$out[[i]]$cen
    global_cen <- x$out[[i]]$global_cen

    # and the outcome database, which should match
    rr <- which(outcomes_tbl[, get(geo_unit_col)] == this_geo)
    single_outcomes_tbl <- outcomes_tbl[rr, ,drop = FALSE]
    if(!identical(single_outcomes_tbl[, get(outcomes_col)], outcomes)) {
      print(head(single_outcomes_tbl))
      print(nrow(single_outcomes_tbl))
      print(head(outcomes))
      print(length(outcomes))
      stop('Outcomes not the same')
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
  stopifnot(all(join_cols %in% names(pop_data)))
  pop_data_collapse <- pop_data[, .(
    population = sum(population)),
    by = join_cols
  ]

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
  xgrid <- tidyr::expand_grid(data.frame(unique_geos),
                       year = unique_years,
                       above_MMT = c(T, F),
                       nsim = 1:nsim)
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
    # rr <- AN_sub_all$this_exp > AN_sub_all$cen)
    # AN_sub_all <- AN_sub_all[rr, ,drop = FALSE]
    AN_sub_all$above_MMT = AN_sub_all$this_exp >= AN_sub_all$cen

    ## do an initial fine grain summary
    group_cols = c(
      geo_unit_col, geo_unit_grp_col, 'year', 'nsim', 'above_MMT'
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

    # fill in NAs
    rr <- which(is.na(AN_ANNUAL[[xi]]$annual_AN))
    AN_ANNUAL[[xi]]$annual_AN[rr] <- 0
    if(any(is.na(AN_ANNUAL[[xi]]))) stop("AN expand_grid didn't work correctly")

    ## and join pop


    AN_ANNUAL[[xi]] <- pop_data_collapse[
      AN_ANNUAL[[xi]],
      on = setNames(join_cols, join_cols)
    ]

    sum_AN_new <- sum(AN_ANNUAL[[xi]]$annual_AN)

    if(!(sum_AN_new == sum_AN_orig)) {
      stop("sum_AN after merge is not the same, some error is happening")
    }


    # and then summarize to the annual level per group based on agg_type
    # this is OK to do within an `xi` because you know its all the data
    if(agg_type == geo_unit_col) {

      invisible(1)

    } else if(agg_type == geo_unit_grp_col) {

      group_cols = c(
        geo_unit_grp_col, 'year', 'nsim', 'above_MMT'
      )

      AN_ANNUAL[[xi]] = AN_ANNUAL[[xi]][,.(
        annual_AN = sum(annual_AN),
        population = sum(population)
      ), by = group_cols]

    } else {

      AN_ANNUAL[[xi]]$all <- 'ALL'

      group_cols = c(
        'all', 'year', 'nsim', 'above_MMT'
      )

      AN_ANNUAL[[xi]] = AN_ANNUAL[[xi]][,.(
        annual_AN = sum(annual_AN),
        population = sum(population)
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
  x1 <- AN_ANNUAL[,.(
    mean_annual_AN = mean(annual_AN)
  ), by = g1_cols]
  if(any(is.na(x1$mean_annual_AN))) {
    rr <- which(is.na(x1$mean_annual_AN))
    print(x1[rr, ])
    stop('annual_AN has NA, find out why')
  }

  # rate column
  x1[, mean_annual_AN_rate := mean_annual_AN / population * 1e5]

  # rate eCI
  g2_cols <- c(names(AN_ANNUAL)[c1], "population", 'above_MMT')
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

#'@export
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

#'@export
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

#'@export
#' plot.calcAN
#'
#' @param x an object of class plot.calcAN
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @importFrom ggplot2 ggplot
#' @importFrom scales number
#' @import data.table
#' @returns
#' @export
#'
#' @examples
plot.calcAN <- function(x, table_type, above_MMT) {

  stopifnot(table_type %in% c('rate', 'num'))
  stopifnot(above_MMT %in% c(T, F))
  geo_unit_col <- attributes(x$`_`)$column_mapping
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  if(table_type == 'num') {

    byX_df <- x$`_`$number_table
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]
    if(nrow(byX_df) > 20) stop("plot elements > 20, not plotting")

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
    if(nrow(byX_df) > 20) stop("plot elements > 20, not plotting")

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


#'@export
#' plot.calcAN_list
#'
#' @param x an object of class plot.calcAN_list
#' @param table_type showing the rate table "rate" or number table "num"
#' @param above_MMT plot attributable numbers above or below the MMT
#' @importFrom ggplot2 ggplot
#' @import data.table
#' @returns
#' @export
#'
#' @examples
plot.calcAN_list <- function(x, table_type, above_MMT) {

  stopifnot(table_type %in% c('rate', 'num'))
  ylab_flag <- if(above_MMT) 'Above MMT' else 'Below MMT'

  fct_names <- names(x)
  fct_lab <- x[[1]]$factor_col

  byX_df  <- vector("list", length(fct_names))

  if(table_type == 'num') {

    for(i in 1:length(byX_df)) {
      byX_df[[i]] <- x[[fct_names[i]]]$`_`$number_table
      byX_df[[i]][, (fct_lab) := fct_names[i]]
    }

    byX_df <- do.call(rbind, byX_df)
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]
    if(nrow(byX_df) > 20 * length(byX_df)) stop("too many plot elements")

    ##
    pd <- position_dodge2(width = 0.9, preserve = "single")


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
    }

    byX_df <- do.call(rbind, byX_df)
    rr <- which(byX_df$above_MMT == above_MMT)
    byX_df <- byX_df[rr, ]
    x_col <- names(byX_df)[1]
    if(nrow(byX_df) > 20 * length(byX_df)) stop("too many plot elements")

    ##
    pd <- position_dodge2(width = 0.9, preserve = "single")

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

#'@export
#' spatial_plot.calcAN
#'
#' @param x an object of class condPois_2stage
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

#'@export
#' spatial_plot.calcAN_list
#'
#' @param x an object of class condPois_2stage_list
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

