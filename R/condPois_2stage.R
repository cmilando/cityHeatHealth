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
#' update as of 2025.12.21, its this https://github.com/gasparrini/Extended2stage/blob/cf82fd2f0d616210e8aa6911d2fb49a732f65795/02.multilevel.R#L68
#'
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes_tbl a data.table of outcomes, created by `make_outcome_table`
#' @param verbose whether to print geo_units as they are run
#' @param global_cen a global centering point
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#' @param strata_min minimum number of cases per strata
#'
#' @importFrom mixmeta mixmeta
#' @importFrom mixmeta blup
#' @importFrom dlnm onebasis
#' @importFrom dlnm crosspred
#' @import data.table
#'
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

    # run this function again for each factor
    for(fct_i in seq_along(fct_outlist)) {

      if(verbose > 0) {
        cat("<",factor_col,":", unique_fcts[fct_i], ">\n")
      }

      rr <- which(outcomes_tbl[, get(factor_col)] == unique_fcts[fct_i])
      subset_outcomes_tbl <- outcomes_tbl[rr, , drop = FALSE]
      attributes(subset_outcomes_tbl)$column_mapping$factor <- NULL

      # re-call the function, but with just one subset
      fct_outlist[[fct_i]] <- condPois_2stage(exposure_matrix,
                                              subset_outcomes_tbl,
                                              global_cen = global_cen,
                                              argvar = argvar,
                                              arglag = arglag,
                                              maxlag = maxlag,
                                              min_n = min_n,
                                              strata_min = strata_min,
                                              verbose = verbose)

      fct_outlist[[fct_i]]$factor_col <- factor_col
      fct_outlist[[fct_i]]$factor_val <- unique_fcts[fct_i]

    }

    # and then make this a "_list" object
    names(fct_outlist) = unique_fcts

    class(fct_outlist) = 'condPois_2stage_list'

    return(fct_outlist)


  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATION
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # Generic validation tests
  validated <- input_validation(exposure_matrix, outcomes_tbl)
  exposure_matrix <- validated$exposure_matrix
  outcomes_tbl    <- validated$outcomes_tbl

  # make objects available
  exp_geo_unit_col     <- attributes(exposure_matrix)$column_mapping$geo_unit
  exp_geo_unit_grp_col <- attributes(exposure_matrix)$column_mapping$geo_unit_grp
  exposure_col         <- attributes(exposure_matrix)$column_mapping$exposure

  out_geo_unit_col     <- attributes(outcomes_tbl)$column_mapping$geo_unit
  out_geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp
  outcome_col          <- attributes(outcomes_tbl)$column_mapping$outcome
  ## CHECK 6 - minN
  if(is.null(min_n)) {
    min_n = 50
  }
  stopifnot(sum(outcomes_tbl[, get(outcome_col)]) >= min_n)

  # CHECK 7
  stopifnot(strata_min >= 0)
  stopifnot(strata_min < min_n)

  # CHECK5
  if(!is.null(global_cen)) {
    stopifnot(is.numeric(global_cen))
  }

  if(verbose > 0) {
    cat("-- validation passed\n")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 1
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- stage 1\n")
  }

  geo_cols <- c(
    out_geo_unit_col,
    out_geo_unit_grp_col
  )
  unique_geos <- unique(outcomes_tbl[, ..geo_cols])
  n_geos      <- nrow(unique_geos)
  stopifnot(n_geos > 1)

  # get things you need
  cp_list     <- vector("list", n_geos);
  vcov_model  <- vector("list", n_geos);
  argvar_list <- vector("list", n_geos);
  cen_list    <- vector("list", n_geos);
  exp_mean    <- vector("list", n_geos);
  exp_IQR     <- vector("list", n_geos);
  vcov_list   <- vector("list", n_geos);
  coef_list   <- vector("list", n_geos);
  outc_list   <- vector("list", n_geos);
  names(vcov_model) <- unique_geos[, get(out_geo_unit_col)]
  names(coef_list) <- unique_geos[, get(out_geo_unit_col)]

  # loop through geos
  for(i in seq_along(cp_list)) {

    # get the name, which you know exists in both datasets
    this_geo <- unique_geos[i, get(out_geo_unit_col)]

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
    d1 <- condPois_1stage(exposure_matrix = single_exposure_matrix,
                          outcomes_tbl = single_outcomes_tbl,
                          argvar = argvar,
                          arglag = arglag,
                          maxlag = maxlag,
                          global_cen = global_cen,
                          min_n = min_n,
                          strata_min = strata_min)

    # get the named list object
    cp_list[[i]] <- d1$`_`$out[[this_geo]]

    # remove cb and basis_cen for space
    # thos are only needed for sb testing and AN calc
    cp_list[[i]]$orig_basis <- NULL
    cp_list[[i]]$basis_cen <- NULL

    # get coef and vcov
    coef_list[[i]] <- coef(cp_list[[i]]$cr)
    vcov_list[[i]] <- vcov(cp_list[[i]]$cr)

    # get argvar_list and cen and outcomes
    argvar_list[[i]]   <- cp_list[[i]]$argvar
    cen_list[[i]]      <- cp_list[[i]]$cen
    outc_list[[i]]     <- cp_list[[i]]$outcomes

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
    cat("-- mixmeta\n")
  }

  # convert to a matrix
  coef_matrix <- do.call(rbind, coef_list)

  # add in the mean and IQR of exposure in each geo_unit
  unique_geos$exp_mean <- do.call(c, exp_mean)
  unique_geos$exp_IQR  <- do.call(c, exp_IQR)
  unique_geos$cen      <- do.call(c, cen_list)

  # create the random effect formula
  # ~ 1 | geo_unit_grp / geo_unit
  # TODO will this work if geo_unit_grp is 'ALL' ?
  #      easy to create a switch if not
  rf = as.formula(paste0("~ 1 | ", out_geo_unit_grp_col, " / ",
                         out_geo_unit_col))

  # see function description for references for this
  meta_fit <- mixmeta(coef_matrix ~ exp_mean + exp_IQR,
                      random = rf,
                      S = vcov_list,
                      data = unique_geos,
                      control = list(showiter=(verbose > 1)),
                      na.action = 'na.exclude')

  # get BLUP
  # there is a `level` argument here
  # pausing this for now, see below
  # blup_geo <- blup(meta_fit, vcov = T)
  # blup_geo_lvl2 <- blup(meta_fit, vcov = T, level = 2)
  # stopifnot(identical(blup_geo, blup_geo_lvl2))
  blup_geo <- blup(meta_fit, vcov = T)
  names(blup_geo) = unique_geos[, get(out_geo_unit_col)]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' detour to calculate blups and RRs for other levels
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # this *should* work as described but it doesn't
  # https://github.com/gasparrini/Extended2stage/issues/1

  # cityblup <- exp(blup(meta_fit, pi=T))
  # stateblup <- unique(blup(meta_fit, level=1))

  # ONE FOR EACH GRP
  group_col <- out_geo_unit_grp_col
  datapred <- unique_geos[, .(
    exp_mean = mean(exp_mean),
    exp_IQR = mean(exp_IQR)
  ), by = group_col]

  ngrp <- nrow(datapred)

  grp_l <- vector("list", ngrp)


  for(grp_i in 1:ngrp) {

    # PREDICT COEF/VCOV
    ## THIS DOESN"T DIRECTLY USE THE BLUP
    ## THIS PREDICTS FROM THE MIXMETA OBJECT
    ## AT THE REGION LEVEL
    pred <- predict(meta_fit, datapred[grp_i, ], vcov = T)

    ## get group_level data
    rr <- which(exposure_matrix[, get(exp_geo_unit_grp_col)] ==
                  datapred[grp_i, get(group_col)])

    grp_dat <- exposure_matrix[rr, ]
    this_exp <- grp_dat[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    ## reset argvar for the group
    grp_argvar <- check_argvar(argvar, this_exp)

    ## get a centered cp
    grp_cp <- get_centered_cp(argvar = grp_argvar,
                              this_exp = this_exp,
                              x_b = x_b,
                              xcoef = pred$fit,
                              xvcov = pred$vcov,
                              global_cen = global_cen,
                              cen = datapred[grp_i, cen])

    grp_l[[grp_i]] <- list(
      cp = grp_cp$cp,
      geo_unit_grp = datapred[grp_i, get(group_col)]
    )

  }

  # conevert list to df
  list_to_df <- function(x) {
    # x = PLOT_l[[1]]
    o <- data.frame(x = x$cp$predvar,
                    RR  =  x$cp$allRRfit,
                    RRlow = x$cp$allRRlow,
                    RRhigh = x$cp$allRRhigh,
                    geo_unit_grp = x$geo_unit_grp)
    return(o)
  }

  plot_df <- do.call(rbind, lapply(grp_l, list_to_df))

  #plot_df

  grp_plt <- ggplot(plot_df) + theme_classic() +
    geom_hline(yintercept = 1, linetype = '11') +
    ##
    geom_ribbon(aes(x = x,
                    ymin = RRlow,
                    ymax = RRhigh),
                fill = 'lightblue', alpha = 0.5) +

    geom_line(aes(x = x,
                  y = RR)) +
    facet_wrap(~geo_unit_grp, axes = 'all') +
    scale_y_continuous(transform = 'log') +
    coord_cartesian(ylim = c(0.75, 1.5)) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = 'bold'))  +
    ylab("Relative Risk") +
    xlab(exposure_col)


  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 2 - BLUP
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- stage 2\n")
  }

  out   <- vector("list", n_geos);

  # loop through geos
  for(i in seq_along(cp_list)) {

    #
    this_geo <- unique_geos[i, get(out_geo_unit_col)]
    this_geo_grp <- unique_geos[i, get(out_geo_unit_grp_col)]

    if(verbose > 1) {
      cat(this_geo, '\t')
    }

    # get x
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    # get centered basis
    blup_cp <- get_centered_cp(argvar = argvar_list[[i]],
                               xcoef = blup_geo[[i]]$blup,
                               xvcov = blup_geo[[i]]$vcov,
                               global_cen = global_cen,
                               cen = cp_list[[i]]$cen,
                               this_exp = this_exp,
                               x_b = x_b)

    ## make the out
    RRdf <- data.frame(
      geo_unit = this_geo,
      geo_unit_grp = this_geo_grp,
      x = blup_cp$cp$predvar,
      RR = blup_cp$cp$allRRfit,
      RRlb = blup_cp$cp$allRRlow,
      RRub = blup_cp$cp$allRRhigh
    )
    names(RRdf)[1] <- out_geo_unit_col
    names(RRdf)[2] <- out_geo_unit_grp_col
    names(RRdf)[3] <- exposure_col
    setDT(RRdf)

    ## attach outcomes vector
    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]
    outcomes_vec = single_outcomes_tbl[, get(outcome_col)]
    stopifnot(identical(outcomes_vec, outc_list[[i]]))

    #
    out[[i]] <- list(
      geo_unit = this_geo,
      basis_cen = blup_cp$basis_cen,
      exposure_col = exposure_col,
      this_exp = this_exp,
      cen = blup_cp$cp$cen,
      global_cen = global_cen,
      outcomes = outc_list[[i]],
      coef = blup_geo[[i]]$blup,
      vcov = blup_geo[[i]]$vcov,
      RRdf = RRdf
    )

  }

  # set names
  names(out) <- unique_geos[, get(out_geo_unit_col)]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # yes, a list in a list, this will make sense later
  # aka, in the recursive call to this function that happens above
  # you could modify this to also output the mixmeta object
  # but not clear that you need that
  outlist = list(list(meta_fit = meta_fit, out = out,
                      grp_plt = grp_plt))
  names(outlist) = "_"
  class(outlist) = 'condPois_2stage'

  return(outlist)

}

#'@export
#' print.condPois_2stage
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_2stage <- function(x) {
  cat("< an object of class `condPois_2stage` >\n")
  invisible(x)
}

#'@export
#' print.condPois_2stage_list
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_2stage_list <- function(x) {
  cat("< an object of class `condPois_2stage_list`:",
      paste(names(x), collapse = ",")," >\n")
  invisible(x)
}

#'@export
#' getRR.condPois_2stage
#'
#' @param x
#' @importFrom data.table setDT
#' @returns
#' @export
#'
#' @examples
getRR.condPois_2stage <- function(x) {
  oo <- do.call(rbind, lapply(x$`_`$out, \(obj) obj$RRdf))
  oo$model_class = class(x)
  setDT(oo)
  return(oo)
}

#'@export
#' getRR.condPois_2stage_list
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
getRR.condPois_2stage_list <- function(x) {

  obj_l <- vector("list", length(names(x)))
  fct_lab <- x[[names(x)[1]]]$factor_col
  exp_lab <- x[[names(x)[1]]]$"_"$out[[1]]$exposure_col

  for(i in seq_along(obj_l)) {
    yy <- do.call(rbind, lapply(x[[names(x)[i]]]$"_"$out, \(obj) obj$RRdf))
    fct <- x[[names(x)[i]]]$factor_val
    yy[[fct_lab]] <- fct
    obj_l[[i]] <- yy
  }

  oo <- do.call(rbind, obj_l)
  oo$model_class = class(x)
  setDT(oo)
  return(oo)

}



#'@export
#' plot.condPois_2stage
#'
#' @param x an object of class condPois_2stage
#' @param geo_unit a geo_unit to investigate
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
plot.condPois_2stage <- function(x, geo_unit,
                                 xlab = NULL, ylab = NULL, title = NULL) {

  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = geo_unit

  obj <- x$`_`$out[[geo_unit]]

  if(is.null(xlab)) xlab = obj$exposure_col

  ggplot(obj$RRdf, aes(x = !!sym(obj$exposure_col), y = RR,
                       ymin = RRlb, ymax = RRub)) +
    geom_hline(yintercept = 1, linetype = '11') +
    scale_y_continuous(transform = 'log') +
    theme_classic() +
    ggtitle(title) +
    geom_ribbon(fill = 'lightblue', alpha = 0.2) +
    geom_line() + xlab(xlab) + ylab(ylab)
}


#'@export
#' plot.condPois_2stage_list
#'
#' @param x an object of class condPois_2stage_list
#' @param geo_unit a geo_unit to investigate
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
plot.condPois_2stage_list <- function(x, geo_unit,
                                 xlab = NULL, ylab = NULL, title = NULL) {

  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = geo_unit

  obj_l <- vector("list", length(names(x)))
  fct_lab <- x[[names(x)[1]]]$factor_col
  exp_lab <- x[[names(x)[1]]]$"_"$out[[geo_unit]]$exposure_col

  for(i in seq_along(obj_l)) {
    yy <- x[[names(x)[i]]]$"_"$out[[geo_unit]]$RRdf
    fct <- x[[names(x)[i]]]$factor_val
    yy[[fct_lab]] <- fct
    obj_l[[i]] <- yy
  }

  obj <- do.call(rbind, obj_l)

  if(is.null(xlab)) xlab = exp_lab

  ggplot(obj) +
    geom_hline(yintercept = 1, linetype = '11') +
    theme_classic() +
    ggtitle(title) +
    geom_ribbon(aes(x = !!sym(exp_lab), y = RR,
                    ymin = RRlb, ymax = RRub,
                    fill = !!sym(fct_lab)), alpha = 0.2) +
    geom_line(aes(x = !!sym(exp_lab), y = RR,
                  color = !!sym(fct_lab))) + xlab(xlab) + ylab(ylab) +
    scale_y_continuous(transform = 'log') +
    scale_color_viridis_d() +
    scale_fill_viridis_d()

}


#'@export
#' forest_plot.condPois_2stage
#'
#' @param x an object of class condPois_2stage
#' @param exposure_val exposure value at which to plot
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
forest_plot.condPois_2stage <- function(x, exposure_val) {

  # get subset of X
  n_geos <- length(x$`_`$out)
  plt_slice <- vector("list", n_geos)
  exposure_col <- x$`_`$out[[1]]$exposure_col
  geo_unit_col <- names(x$`_`$out[[1]]$RRdf)[1]
  geo_unit_grp_col <- names(x$`_`$out[[1]]$RRdf)[2]

  for(i in 1:n_geos) {
    rr <- which(x$`_`$out[[i]]$RRdf[[exposure_col]] == exposure_val)
    if(length(rr) != 1) {
      stop(paste0("Exposure value '", exposure_val, "' not in the
                exposure column, try values (with one decimal) between:",
                  min(x$`_`$out[[i]]$RRdf[[exposure_col]]),
                  " and ",
                  max(x$`_`$out[[i]]$RRdf[[exposure_col]])))
    }

    plt_slice[[i]] <- x$`_`$out[[i]]$RRdf[rr, ]
  }

  plt_slice <- do.call(rbind, plt_slice)

  # forest_plot
  if(n_geos > 20) {
    warning("plotting by group since n_geos > 20")
    ggplot(plt_slice, aes(x = RR, xmin = RRlb, xmax = RRub,
                          y = reorder(!!sym(geo_unit_grp_col), RR))) +
      geom_vline(xintercept = 1.0, linetype = '11') +
      ylab(geo_unit_grp_col) +
      theme_classic() +
      scale_x_continuous(transform = 'log') +
      geom_pointrange(position = position_jitterdodge()) +
      ggtitle(paste0(exposure_col, " = ", exposure_val))
  } else {
    ggplot(plt_slice, aes(x = RR, xmin = RRlb, xmax = RRub,
                          y = reorder(!!sym(geo_unit_col), RR))) +
      geom_vline(xintercept = 1.0, linetype = '11') +
      theme_classic() +
      scale_x_continuous(transform = 'log') +
      geom_pointrange() +
      ggtitle(paste0(exposure_col, " = ", exposure_val))
  }

}


#'@export
#' spatial_plot.condPois_2stage
#'
#' @param x an object of class condPois_2stage
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param exposure_val exposure value at which to plot
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_2stage <- function(x, shp, exposure_val,
                                         RRlimits = NULL) {

  n_geos <- length(x$`_`$out)
  plt_slice <- vector("list", n_geos)
  exposure_col <- x$`_`$out[[1]]$exposure_col
  geo_unit_col <- names(x$`_`$out[[1]]$RRdf)[1]
  geo_unit_grp_col <- names(x$`_`$out[[1]]$RRdf)[2]

  for(i in 1:n_geos) {
    rr <- which(x$`_`$out[[i]]$RRdf[[exposure_col]] == exposure_val)
    if(length(rr) != 1) {
      stop(paste0("Exposure value '", this_x, "' not in the
                  exposure column, try values (with one decimal) between:",
                  min(x$`_`$out[[i]]$RRdf[[exposure_col]]),
                  " and ",
                  max(x$`_`$out[[i]]$RRdf[[exposure_col]])))
    }

    plt_slice[[i]] <- x$`_`$out[[i]]$RRdf[rr, ]
  }

  plt_slice <- do.call(rbind, plt_slice)

  # join to sf object
  stopifnot(geo_unit_col %in% names(shp)) # not a bad first check
  shp_w_data <- merge(shp, plt_slice)

  if(is.null(RRlimits)) RRlimits = range(shp_w_data$RR)


  return(ggplot(shp_w_data) +
    theme_classic() +
    geom_sf(aes(fill = log(RR))) +
    scale_fill_gradient2(
      midpoint = 0,
      limits = log(RRlimits),
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      name = "RR",
      labels = function(x) round(exp(x), 2)
    ) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank()) +
    ggtitle(paste0(exposure_col, " = ", exposure_val)))


}

#'@export
#' spatial_plot.condPois_2stage_list
#'
#' @param x an object of class condPois_2stage_list
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param exposure_val exposure value at which to plot
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_2stage_list <- function(x, shp, exposure_val) {

  obj_l <- vector("list", length(names(x)))
  fct_lab <- x[[names(x)[1]]]$factor_col

  for(i in 1:length(names(x))) {

    yy <- x[[names(x)[i]]]$`_`$out

    n_geos <- length(yy)
    plt_slice <- vector("list", n_geos)
    exposure_col <- yy[[1]]$exposure_col
    geo_unit_col <- names(yy[[1]]$RRdf)[1]
    geo_unit_grp_col <- names(yy[[1]]$RRdf)[2]

    for(j in 1:n_geos) {
      rr <- which(yy[[j]]$RRdf[[exposure_col]] == exposure_val)
      if(length(rr) != 1) {
        stop(paste0("Exposure value '", this_x, "' not in the
                  exposure column, try values (with one decimal) between:",
                    min(yy[[j]]$RRdf[[exposure_col]]),
                    " and ",
                    max(yy[[j]]$RRdf[[exposure_col]])))
      }

      plt_slice[[j]] <- yy[[j]]$RRdf[rr, ]
    }

    plt_slice <- do.call(rbind, plt_slice)
    plt_slice[[fct_lab]] <- x[[names(x)[i]]]$factor_val

    obj_l[[i]] <- plt_slice
  }

  obj <- do.call(rbind, obj_l)
  RRlimits = range(obj$RR)

  obj_split <- split(obj, f = obj[[fct_lab]])

  plt_obj <- vector("list", length(names(x)))
  for(i in 1: length(names(x))) {
    plt_obj[[i]] <- spatial_plot(x[[names(x)[i]]], shp,
                                 exposure_val, RRlimits = RRlimits) +
      ggtitle(paste0(fct_lab,": ", names(x)[i],"; ",
                     exposure_col, " = ", exposure_val))
  }

  patchwork::wrap_plots(plt_obj, ncol = 1,guides = 'collect')

}


#'@export
#' forest_plot.condPois_2stage_list
#'
#' @param x an object of class condPois_2stage_list
#' @param exposure_val exposure value at which to plot
#' @returns
#' @export
#'
#' @examples
forest_plot.condPois_2stage_list <- function(x, exposure_val) {

  obj_l <- vector("list", length(names(x)))
  fct_lab <- x[[names(x)[1]]]$factor_col

  for(i in 1:length(names(x))) {
    yy <- x[[names(x)[i]]]$"_"$out
    n_geos <- length(yy)
    plt_slice <- vector("list", n_geos)
    fct <- x[[names(x)[i]]]$factor_val
    exposure_col <- yy[[1]]$exposure_col
    geo_unit_col <- names(yy[[1]]$RRdf)[1]
    geo_unit_grp_col <- names(yy[[1]]$RRdf)[2]

    for(j in 1:n_geos) {
      rr <- which(yy[[j]]$RRdf[[exposure_col]] == exposure_val)
      if(length(rr) != 1) {
        stop(paste0("Exposure value '", exposure_val, "' not in the
                exposure column, try values (with one decimal) between:",
                    min(yy[[j]]$RRdf[[exposure_col]]),
                    " and ",
                    max(yy[[j]]$RRdf[[exposure_col]])))
      }
      plt_slice[[j]] <- yy[[j]]$RRdf[rr, ]
    }

    plt_slice <- do.call(rbind, plt_slice)
    plt_slice[[fct_lab]] <- fct

    obj_l[[i]] <- plt_slice
  }

  obj <- do.call(rbind, obj_l)

  # forest_plot
  if(n_geos > 20) {
    warning("plotting by group since n_geos > 20")
    ggplot(obj, aes(x = RR, xmin = RRlb, xmax = RRub,
                          color = !!sym(fct_lab),
                          y = reorder(!!sym(geo_unit_grp_col), RR))) +
      geom_vline(xintercept = 1.0, linetype = '11') +
      ylab(geo_unit_grp_col) +
      theme_classic() +
      scale_color_viridis_d() +
      scale_x_continuous(transform = 'log') +
      geom_pointrange(position = position_jitterdodge()) +
      ggtitle(paste0(exposure_col, " = ", exposure_val))
  } else {
    ggplot(obj, aes(x = RR, xmin = RRlb, xmax = RRub,
                          y = reorder(!!sym(geo_unit_col), RR))) +
      geom_vline(xintercept = 1.0, linetype = '11') +
      theme_classic() +
      scale_x_continuous(transform = 'log') +
      geom_pointrange() +
      scale_color_viridis_d() +
      ggtitle(paste0(exposure_col, " = ", exposure_val))
  }
}
