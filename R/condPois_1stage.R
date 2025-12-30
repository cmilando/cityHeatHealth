#' Run a conditional poisson model for a single geographic unit
#'
#' @import data.table
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom dlnm crossreduce
#' @importFrom dlnm logknots
#' @importFrom gnm gnm
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes_tbl a data.table of outcomes, created by `make_outcome_table`
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#' @param strata_min an integer describing the minimum number of cases for a single strata
#' @param global_cen global centering point
#' @param multi_zone are multiple strata being used.
#'
#' @returns
#' @export
#'
#' @examples
condPois_1stage <- function(exposure_matrix, outcomes_tbl,
                        argvar = NULL, arglag = NULL, maxlag = NULL,
                       min_n = NULL, strata_min = 0, global_cen = NULL,
                       multi_zone = FALSE) {

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

      # cat("<",factor_col,":", unique_fcts[fct_i], ">\n")
      rr <- which(outcomes_tbl[, get(factor_col)] == unique_fcts[fct_i])
      subset_outcomes_tbl <- outcomes_tbl[rr, , drop = FALSE]
      attributes(subset_outcomes_tbl)$column_mapping$factor <- NULL

      # re-call the function, but with just one subset
      fct_outlist[[fct_i]] <- condPois_1stage(exposure_matrix,
                                              subset_outcomes_tbl,
                                              global_cen = global_cen,
                                              argvar = argvar,
                                              arglag = arglag,
                                              maxlag = maxlag,
                                              min_n = min_n,
                                              strata_min = strata_min,
                                              multi_zone = multi_zone)

      fct_outlist[[fct_i]]$factor_col <- factor_col
      fct_outlist[[fct_i]]$factor_val <- unique_fcts[fct_i]

    }

    names(fct_outlist) = unique_fcts

    class(fct_outlist) = 'condPois_1stage_list'

    return(fct_outlist)


  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## Check 1.5 -- there should be just one geo_unit
  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  out_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  exp_geo_unit <- unlist(unique(exposure_matrix[, get(exp_geo_unit_col)]))
  out_geo_unit <- unlist(unique(outcomes_tbl[, get(out_geo_unit_col)]))

  stopifnot(identical(exp_geo_unit, out_geo_unit))

  # check if multizone
  if(!multi_zone) stopifnot(length(out_geo_unit) == 1)

  ## Check 2
  ## probably should make sure that exposure_matrix and outcomes_tbl
  ## are the same size, at least
  ## and have the same dates
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  # subset so its a complete match
  rr <- which(exposure_matrix[, get(exp_date_col)] %in%
                outcomes_tbl[, get(outcome_date_col)])
  exposure_matrix <- exposure_matrix[rr, ,drop = FALSE]

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

  # CHECK 4
  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }

  ## CHECK 6 - minN
  if(is.null(min_n)) {
    min_n = 50
  }

  outcome_col <- attributes(outcomes_tbl)$column_mapping$outcome
  stopifnot(sum(outcomes_tbl[, get(outcome_col)]) >= min_n)

  # CHECK 7
  stopifnot(strata_min >= 0)
  stopifnot(strata_min < min_n)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CREATE CROSSBASIS for this single zone
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////
  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure

  # maxlag
  if(is.null(maxlag)) {
    maxlag = 5
  } else {
    stopifnot(maxlag %in% 2:10)
  }

  # argvar
  this_exp = exposure_matrix[, get(exposure_col)]
  argvar <- check_argvar(argvar, this_exp)

  # arglag
  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = logknots(maxlag, nk = 2))
  } else {
    warning("check that arglag is valid")
  }

  ## get the columns you need
  xcols <- c(exposure_col, paste0('explag',1:maxlag))
  x_mat <- exposure_matrix[, ..xcols]

  ## if you are safe to proceed, make the x_mat
  ## since you are passing in a matrix, you dont need to do
  ## group = year or location, because all the data you need
  ## are in each row, and this also means the order doesn't matter
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  # there should be no NAs here
  if(any(is.na(cb))) stop("crossbasis has NULL, something went wrong")

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' RUN GNM
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## if using GNM, you get COEF and VCOV as part of the model objects
  ## instead of doing keep == 1, i updated to be a strata total variable
  ## so you can raise the floor if you want to.
  ##
  ## Note -- you also don't need to do offset = log(poplation)
  ## because the population is not changing within the strata,
  ## and since you are doing conditional poisson
  ## it all cancels out in the numerator and denominator of the
  ## multi-nomial distribution
  ff = as.formula(paste(outcome_col, "~ cb"))

  m_sub <- gnm(formula = ff,
               data = outcomes_tbl,
               family = quasipoisson,
               eliminate = factor(strata),
               subset = strata_total > strata_min)

  m_coef <- coef(m_sub)
  m_vcov <- vcov(m_sub)

  # there should be no NAs
  if(any(is.na(m_coef))) stop("coef has NULL, something went wrong")
  if(any(is.na(m_vcov))) stop("vcov has NULL, something went wrong")

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CROSSPRED and CROSSREDUCE
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit
  geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp

  exp_mean = mean(exposure_matrix[, get(exposure_col)])
  exp_IQR = IQR(exposure_matrix[, get(exposure_col)])

  # the crossreduce coefficients are not affected by the centering point
  # but it does make a message if you dont put something there
  # so i center on the min to avoid that fate
  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = exp_mean,
                  by = 0.1)

  if(!is.null(global_cen)) {
    cen = global_cen
  } else {
    cen = cp$predvar[which.min(cp$allRRfit)]
  }

  # now apply to cr and export
  cr <- crossreduce(cb,
                    coef = m_coef,
                    vcov = m_vcov,
                    model.link = "log",
                    cen = cen,
                    by = 0.1)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT OBJECTS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  outcome_columns <- attributes(outcomes_tbl)$column_mapping

  geo_cols <- c(
    outcome_columns$geo_unit,
    outcome_columns$geo_unit_grp
  )
  unique_geos <- unique(outcomes_tbl[, ..geo_cols])

  oo_list <- vector("list", length(out_geo_unit))
  names(oo_list) <- unique_geos[, get(out_geo_unit_col)]

  for(i in 1:nrow(unique_geos)) {

    # get the name, which you know exists in both datasets
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]
    this_geo_grp <- unique_geos[i, get(outcome_columns$geo_unit_grp)]

    # this cities exposure matrix
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))
    this_exp_mean = mean(single_exposure_matrix[, get(exposure_col)])
    this_exp_IQR = IQR(single_exposure_matrix[, get(exposure_col)])

    # this cities cb
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    this_cb <- cb[rr, ]

    # this cities outcome
    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]
    outcomes <- single_outcomes_tbl[, get(outcome_col)]

    # and get centered
    basis1z <- get_centered_cp(argvar = argvar,
                               xcoef = coef(cr),
                               xvcov = vcov(cr),
                               global_cen = global_cen,
                               cen = cen,
                               this_exp = this_exp,
                               x_b = x_b)

    # each of these things you need for BLUP and MIXMETA later
    oo_list[[i]] <- list(geo_unit = this_geo,     ## individual
               geo_unit_grp = this_geo_grp,       ## individual
               basis_cen = basis1z$basis_cen,
               strata_vec = single_outcomes_tbl$strata,
               orig_basis = this_cb,
               orig_coef = m_coef,
               orig_vcov = m_vcov,
               cr = cr,                           ## whole group
               coef = coef(cr),                ## whole group
               vcov = vcov(cr),                ## whole group
               exposure_col = exposure_col,       ## whole group
               this_exp = this_exp,               ## for the individual
               outcomes = outcomes,               ## for the individual
               cen = cen,                         ## whole group
               global_cen = global_cen,           ## whole group
               argvar = argvar,                   ## whole group
               exp_mean = exp_mean,               ## for the whole group
               exp_IQR = exp_IQR)                 ## for the whole group

  }

  outlist = list(list(out = oo_list))
  names(outlist) = "_"
  class(outlist) <- 'condPois_1stage'

  return(outlist)

}



#'@export
#' print.condPois_1stage
#'
#' @param x an object of class condPois_1stage
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_1stage <- function(x) {
  cat("< an object of class `condPois_1stage` >\n")
  invisible(x)
}


#'@export
#' print.condPois_1stage_list
#'
#' @param x an object of class condPois_1stage_list
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_1stage_list <- function(x) {
  cat("< an object of class `condPois_1stage_list`:",
      paste(names(x), collapse = ",")," >\n")
}


#'@export
#' plot.condPois_1stage
#'
#' @param x an object of class condPois_1stage
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @import ggplot2
#' @returns
#' @export
#'
#' @examples
plot.condPois_1stage <- function(x, xlab = NULL, ylab = NULL, title = NULL) {

  n_geo_names <- paste0(names(x$`_`$out), collapse = ":")

  # these will all be the same, so just pick 1
  plot_cp = data.frame(
    x = x$`_`$out[[1]]$cr$predvar,
    RR = x$`_`$out[[1]]$cr$RRfit,
    RRlow = x$`_`$out[[1]]$cr$RRlow,
    RRhigh = x$`_`$out[[1]]$cr$RRhigh
  )

  if(is.null(xlab)) xlab = x$`_`$out[[1]]$exposure_col
  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = n_geo_names

  ggplot(plot_cp, aes(x = x, y = RR, ymin = RRlow, ymax = RRhigh)) +
    geom_hline(yintercept = 1, linetype = '11') +
    theme_classic() +
    ggtitle(title) +
    scale_y_continuous(transform = 'log') +
    geom_ribbon(fill = 'lightblue', alpha = 0.2) +
    geom_line() + xlab(xlab) + ylab(ylab)
}


#'@export
#' plot.condPois_1stage
#'
#' @param x an object of class condPois_1stage
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @import ggplot2
#' @returns
#' @export
#'
#' @examples
plot.condPois_1stage_list <- function(x, xlab = NULL, ylab = NULL, title = NULL) {

  fct_names <- names(x)

  plot_cl_l <- vector("list", length(names(x)))

  for(i in 1:length(names(x))) {

    # these will all be the same so just pick the first one
    plot_cl_l[[i]] = data.frame(
      x = x[[names(x)[i]]]$`_`$out[[1]]$cr$predvar,
      fct = names(x)[i],
      RR = x[[names(x)[i]]]$`_`$out[[1]]$cr$RRfit,
      RRlow = x[[names(x)[i]]]$`_`$out[[1]]$cr$RRlow,
      RRhigh = x[[names(x)[i]]]$`_`$out[[1]]$cr$RRhigh
    )

    n_geo_names <- paste0(names(x[[names(x)[i]]]$`_`$out), collapse = ":")
    factor_col <- x[[names(x)[i]]]$factor_col
    names(plot_cl_l[[i]])[2] <- factor_col

  }

  plot_cp <- do.call(rbind, plot_cl_l)

  if(is.null(xlab)) xlab = x[[1]]$exposure_col
  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = n_geo_names

  ggplot(plot_cp,
         aes(x = x, y = RR, ymin = RRlow, ymax = RRhigh)) +
    geom_hline(yintercept = 1, linetype = '11') +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme_classic() +
    ggtitle(title) +
    scale_y_continuous(transform = 'log') +
    geom_ribbon(aes(fill = !!sym(factor_col)), alpha = 0.2) +
    geom_line(aes(color = !!sym(factor_col))) +
    xlab(xlab) + ylab(ylab)

}



#'@export
#' spatial_plot.condPois_1stage
#'
#' @param x an object of class condPois_1stage
#' @param ... other elements passed to spatial_plot
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_1stage <- function(x, ...) {
  warning("`spatial_plot` method not implemented for objects of class `condPois_1stage`,
      since there is only one 1_stage relative risk curve so all plot
      values would be the same. 1stage attributable number results will change
      over space, so those can be viewed instead by running `spatial_plot` on the
      output of `calcAN` for a 1stage model!")
}

#'@export
#' spatial_plot.condPois_1stage_list
#'
#' @param x an object of class condPois_1stage
#' @param ... other elements passed to spatial_plot
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_1stage_list <- function(x, ...) {
  warning("`spatial_plot` method not implemented for objects of class `condPois_1stage_list`,
      since there is only one 1_stage relative risk curve so all plot
      values would be the same. 1stage attributable number results will change
      over space, so those can be viewed instead by running `spatial_plot` on the
      output of `calcAN` for a 1stage model!")
}
