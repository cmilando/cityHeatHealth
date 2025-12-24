#' Run a conditional poisson model with Spatial Bayesian inference
#'
#' https://mc-stan.org/docs/2_20/functions-reference/multinomial-distribution.html
#'
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes_tbl a data.table of outcomes, created by `make_outcome_table`
#' @param shp_sf a shapefile object describing the geo_units
#' @param stan_opts a list of parameters to send to stan
#' @param verbose whether to print geo_units as they are run
#' @param global_cen a global centering point
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for a single region
#' @param strata_min minimum number of cases per strata
#' @importFrom mixmeta mixmeta
#' @importFrom mixmeta blup
#' @importFrom dlnm onebasis
#' @importFrom dlnm crosspred
#' @import data.table
#' @import spdep
#' @import rstan
#' @import sf
#' @returns
#' @export
#'
#' @examples
condPois_sb <- function(exposure_matrix,
                        outcomes_tbl,
                        shp_sf,
                        stan_opts,
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

      # re-call the function, but with just one subset
      fct_outlist[[fct_i]] <- condPois_sb(exposure_matrix,
                                          subset_outcomes_tbl,
                                          shp,
                                          stan_opts,
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

    names(fct_outlist) = unique_fcts

    class(fct_outlist) = 'condPois_sb_list'

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
    cat("-- validation passed\n")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 1
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- prepare inputs\n")
  }

  outcome_columns <- attributes(outcomes_tbl)$column_mapping
  exposure_columns <- attributes(exposure_matrix)$column_mapping

  geo_cols <- c(
    outcome_columns$geo_unit,
    outcome_columns$geo_unit_grp
  )
  unique_geos <- unique(outcomes_tbl[, ..geo_cols])
  n_geos      <- nrow(unique_geos)
  stopifnot(n_geos > 1)

  # get things you need
  cb_list     <- vector("list", n_geos);
  outc_list   <- vector("list", n_geos);
  cen_list    <- vector("list", n_geos);
  argvar_list <- vector("list", n_geos);
  coef_list   <- vector("list", n_geos);

  # loop through geos
  for(i in 1:n_geos) {

    # get the name, which you know exists in both datasets
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]

    if(verbose > 1) {
      cat(this_geo, '\t')
    }

    # this cities exposure matrix
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_columns$exposure)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]

    # run the model once
    # just using this to get cb so i don't repeat that code twice
    # because it also has a lot of nice validation built into it
    # min_n is set to 0 because you want every region
    # same with strata_min
    local_cp <- condPois_1stage(exposure_matrix = single_exposure_matrix,
                                    outcomes_tbl = single_outcomes_tbl,
                                    argvar = argvar, arglag = arglag,
                                    maxlag = maxlag, min_n = 1,
                                    strata_min = 0)

    blup_cp <- get_centered_cp(argvar = local_cp$argvar,
                               xcoef = local_cp$cr_coef,
                               xvcov = local_cp$cr_vcov,
                               global_cen = global_cen,
                               cen = local_cp$cen,
                               this_exp = this_exp,
                               x_b = x_b)

    # get cb and outcomes lists
    coef_list[[i]]   <- local_cp$cr_coef
    cb_list[[i]]     <- blup_cp$basis_cen
    outc_list[[i]]   <- single_outcomes_tbl
    cen_list[[i]]    <- local_cp$cen
    argvar_list[[i]] <- local_cp$argvar

    rm(local_cp)
  }

  #' ////////////////////////////////////////////////////////////////////////////
  #' ============================================================================
  #' SETUP inputs
  #' ============================================================================
  #' #' /////////////////////////////////////////////////////////////////////////

  # n regions
  J = as.integer(n_geos)

  # nrows --> this has to be the same per region right
  n_l <- lapply(outc_list, \(x) nrow(x))
  n_l <- unique(do.call(c, n_l))
  stopifnot(length(n_l) == 1)
  N = as.integer(n_l)

  # beta values, withouth the intercept
  k_l <- lapply(cb_list, \(x) ncol(as.matrix(x)))
  k_l <- unique(do.call(c, k_l))
  stopifnot(length(k_l) == 1)
  K = as.integer(k_l)
  if(K > 5) warning("K > 5, times may be slow")

  # include the intercept
  X = array(dim = c(dim(cb_list[[1]]), J))
  for(j in 1:J) X[,,j] = as.matrix(cb_list[[j]])

  # outcome in two regions
  y = array(dim = c(n_l, J))
  for(j in 1:J) y[,j] = outc_list[[j]][, get(outcome_columns$outcome)]

  # create S matrix
  getSmat <- function(strata_vector, include_self = T) {

    # strata_vector <- data$stratum
    strata_matrix <- matrix(as.integer(strata_vector),
                            nrow = length(strata_vector),
                            ncol = length(strata_vector),
                            byrow = T)

    for(i in 1:length(strata_vector)) {
      strata_matrix[i, ] = 1*(strata_matrix[i, ] == as.integer(strata_vector[i]))
      if(!include_self) {
        strata_matrix[i, i] = 0
      }
    }

    return(strata_matrix)
  }

  #
  S <- getSmat(factor(outc_list[[1]]$strata), include_self = T)
  stopifnot(any(S == 1))

  # get strata vars
  n_strata <- as.integer(length(unique(outc_list[[1]]$strata))) # 72, cool
  S_list <- apply(S, 1, function(x) which(x == 1))
  max_in_strata <- max(sapply(S_list, length))
  S_list <- lapply(S_list, function(l) {
    if(length(l) == max_in_strata) {
      return(l)
    } else {
      diff_n = max_in_strata - length(l)
      return(c(l, rep(0, times = diff_n)))
    }
  })
  S_condensed <- unique(do.call(rbind, S_list))
  dim(S_condensed)

  #
  stratum_id = as.integer(factor(outc_list[[1]]$strata))
  stratum_id

  ## *** THIS BECOMES THE J MATRIX
  # ok now do the same as include self but with J matrix
  # start simple and you can generalize later
  # this is an adjacency matrix that DOES NOT include itself
  # shp_sf <- read_sf(shp)

  shp_sf_safe <- subset(shp_sf, COUNTY20 %in% unique(outcomes_tbl$COUNTY20))
  shp_sf_safe
  # BUT YOU DON'T KNOW THE ORDER OF THEM .... ..... ..... .....

  shp_sf_safe <- shp_sf_safe[order(shp_sf_safe$COUNTY20), ]
  stopifnot(identical(shp_sf_safe$COUNTY20, unique_geos$COUNTY20))

  # Generate a list of spatial structure from the shapefile
  list_neig <- nb2listw(poly2nb(shp_sf_safe), zero.policy = TRUE)
  stopifnot(all(dim(list_neig)))


  # 1st-order neighbors
  nb1 <- poly2nb(shp_sf_safe)

  # 2nd-order neighbors (neighbors of neighbors)
  nb_lags <- nblag(nb1, 2)

  # nb_lags[[1]] = distance 1
  # nb_lags[[2]] = distance 2

  nb12 <- vector("list", n_geos)
  for (i in seq_along(nb12)) {
    nb12[[i]] <- list(n1 = nb_lags[[1]][[i]],
                      n2 = nb_lags[[2]][[i]])
  }

  nb12




  neighbors <- lapply(list_neig$neighbours,c)
  Jmat <- matrix(0, nrow = J, ncol = J)
  if(j > 1) {
    for(j in 1:J) {
      print(j)
      n1 = nb_lags[[1]][[j]]
      n2 = nb_lags[[2]][[j]]
      n1
      n2

      Jmat[j, n1] <- 1
      Jmat[j, n2] <- 1
    }
  }
  Jmat

  #' ////////////////////////////////////////////////////////////////////////////
  #' ============================================================================
  #' Run STAN
  #' ============================================================================
  #' #' /////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- run STAN\n")
  }

  # Note: Be REALLY careful about the data types here
  stan_data <- list(
    J = as.integer(J),  # integer
    Jmat = Jmat,        # matrix so you can do math on it
    N = as.integer(N),
    K = as.integer(K),
    X = X,              # matrix so you can do math
    y = y,
    S = S,              # matrix so you can do math on it
    n_strata = n_strata,
    max_in_strata = max_in_strata,
    S_condensed = S_condensed,
    stratum_id = stratum_id
  )

  # from here https://github.com/rok-cesnovar/misc/blob/master/democmdstanr/R/bernoulli.R
  # and here https://discourse.mc-stan.org/t/r-package-using-cmdstanr/32758
  stan_file <- system.file("stan", "SB_CondPoisson.stan",
                           package = "cityHeatHealth")

  mod <- cmdstanr::cmdstan_model(stan_file,
                                 cpp_options = list(stan_threads = TRUE))

  # (1) Variational
  out2_var <- mod$variational(
    data = stan_data,
    iter = 10000,
    refresh = 10,
    threads = 2
  )

  # (2) LaPLACE

  # MCMC -- Much longer, especially for things with lags
  #         wait ..... can  you do this just on the crossreduce though ...
  #         why do you need to do this on
  out2_mcmc <- mod$sample(
    data = stan_data,
    chains = 2,
    #iter_warmup = 3000,
    #iter_sampling = 5000,
    parallel_chains = 2,
    threads_per_chain = 1,
    refresh = 10,
    #adapt_delta = 0.8,
    #max_treedepth = 7 # .... ?
  )

  # or another version here for mcmc

  # or laplace i guess


  # ****************
  # **** COEF ******
  # ****************

  # and then finally the draws etc
  out2 <- out2_var
  out2 <- out2_mcmc
  out2_data <- posterior::as_draws_df(out2$draws())
  setDT(out2_data)
  beta_reg_all <- out2_data[, .SD, .SDcols = patterns("^beta_out")]
  beta_reg_all <- apply(beta_reg_all, 2, mean)
  beta_mat <- matrix(beta_reg_all, nrow = K, ncol = J, byrow = F)

  # with spatial smoothing
  t(beta_mat)

  # without spatial smoothing
  do.call(rbind, coef_list)


  # and using more threads? doesn't seem like comp is working hard rn
  # but great progress! almost there :)

  # and obv you have to compare to the non-SB data but this is pretty good

  # ****************
  # **** VCOV ******
  # ****************

  # so this gives you the coef
  # but then you also need to get VCOV
  # calc_vcov <- function(y, X, beta, stratum_vector)
  # calc_dispersion <- function(y, X, beta, stratum_vector)
  #
  # dispersion <- calc_dispersion(y = boston_deaths_tbl$daily_deaths,
  #                               X = cb,
  #                               stratum_vector = boston_deaths_tbl$strata,
  #                               beta = coef(m_sub))
  #
  # vcov_beta <- calc_vcov(y = boston_deaths_tbl$daily_deaths,
  #                        X = cb,
  #                        stratum_vector = boston_deaths_tbl$strata,
  #                        beta = coef(m_sub))
  #
  # # update with the dipserion parameter
  # t1 <- vcov_beta * dispersion


  # xcoef = sb_geo[[i]]$coef
  # xvcov = sb_geo[[i]]$vcov

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' STAGE 2 - GET FINAL ESTIMATES
  #'
  #' sb geo needs coef and vcov
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- apply estimates\n")
  }

  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure
  outcomes_col <- attributes(outcomes_tbl)$column_mapping$outcome

  blup_out   <- vector("list", n_geos);

  # loop through geos
  for(i in 1:n_geos) {

    #
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]

    # get x
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    # get centered basis
    blup_cp <- get_centered_cp(argvar = argvar_list[[i]],
                               xcoef = sb_geo[[i]]$coef,
                               xvcov = sb_geo[[i]]$vcov,
                               global_cen = global_cen,
                               cen = cen_list[[i]],
                               this_exp = this_exp,
                               x_b = x_b)

    ## make the out
    RRdf <- data.frame(
      geo_unit = this_geo,
      x = blup_cp$cp$predvar,
      RR = blup_cp$cp$allRRfit,
      RRlb = blup_cp$cp$allRRlow,
      RRub = blup_cp$cp$allRRhigh
    )
    names(RRdf)[2] <- exposure_col
    setDT(RRdf)

    ## attach outcomes vector
    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]
    outcomes_vec = single_outcomes_tbl[, get(outcomes_col)]
    stopifnot(sum(outcomes_vec) == sum(outc_list[[i]]))

    #
    blup_out[[i]] <- list(
      geo_unit = this_geo,
      basis_cen = blup_cp$basis_cen,
      exposure_col = exposure_col,
      this_exp = this_exp,
      cen = blup_cp$cp$cen,
      global_cen = global_cen,
      outcomes = outc_list[[i]],
      coef = sb_geo[[i]]$coef,
      vcov = sb_geo[[i]]$vcov,
      RRdf = RRdf
    )

  }

  # set names
  names(blup_out) <- unique_geos[, get(outcome_columns$geo_unit)]

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # yes, a list in a list, this will make sense later
  # aka, in the recursive call to this function that happens above
  # you could modify this to also output the mixmeta object
  # but not clear that you need that
  outlist = list(list(meta_fit = meta_fit, blup_out = blup_out,
                      grp_plt = grp_plt))
  names(outlist) = "_"
  class(outlist) = 'condPois_sb'

  return(outlist)

}

#'@export
#' print.condPois_sb
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_sb <- function(x) {
  cat("< an object of class `condPois_sb` >\n")
  invisible(x)
}

#'@export
#' print.condPois_sb_list
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
print.condPois_sb_list <- function(x) {
  cat("< an object of class `condPois_sb_list`:",
      paste(names(x), collapse = ",")," >\n")
  invisible(x)
}

#'@export
#' plot.condPois_sb
#'
#' @param x an object of class condPois_sb
#' @param geo_unit a geo_unit to investigate
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
plot.condPois_sb <- function(x, geo_unit,
                                 xlab = NULL, ylab = NULL, title = NULL) {

  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = geo_unit

  obj <- x$`_`$blup_out[[geo_unit]]

  if(is.null(xlab)) xlab = obj$exposure_col

  ggplot(obj$RRdf, aes(x = !!sym(obj$exposure_col), y = RR,
                       ymin = RRlb, ymax = RRub)) +
    geom_hline(yintercept = 1, linetype = '11') +
    theme_classic() +
    ggtitle(title) +
    geom_ribbon(fill = 'lightblue', alpha = 0.2) +
    geom_line() + xlab(xlab) + ylab(ylab)
}


#'@export
#' plot.condPois_sb_list
#'
#' @param x an object of class condPois_sb_list
#' @param geo_unit a geo_unit to investigate
#' @param xlab xlab override
#' @param ylab ylab override
#' @param title title override
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
plot.condPois_sb_list <- function(x, geo_unit,
                                      xlab = NULL, ylab = NULL, title = NULL) {

  if(is.null(ylab)) ylab = "RR"
  if(is.null(title)) title = geo_unit

  obj_l <- vector("list", length(names(x)))
  fct_lab <- x[[names(x)[1]]]$factor_col
  exp_lab <- x[[names(x)[1]]]$"_"$blup_out[[geo_unit]]$exposure_col

  for(i in seq_along(obj_l)) {
    yy <- x[[names(x)[i]]]$"_"$blup_out[[geo_unit]]$RRdf
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
    scale_color_viridis_d() +
    scale_fill_viridis_d()

}

