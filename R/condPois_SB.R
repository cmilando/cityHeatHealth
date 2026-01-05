#' Run a conditional poisson model with Spatial Bayesian inference
#'
#' https://mc-stan.org/docs/2_20/functions-reference/multinomial-distribution.html
#'
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes_tbl a data.table of outcomes, created by `make_outcome_table`
#' @param shp_sf a shapefile object describing the geo_units
#' @param stan_type a choice of either laplace or sample
#' @param stan_opts a list of parameters to send to stan
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
#' @importFrom spdep poly2nb
#' @import data.table
#'
#' @returns
#' @export
#'
#' @examples
condPois_sb <- function(exposure_matrix,
                        outcomes_tbl,
                        shp_sf,
                        stan_type = 'laplace',
                        use_spatial_model = 'bym2',
                        stan_opts = NULL,
                        global_cen = NULL,
                        argvar = NULL,
                        arglag = NULL,
                        maxlag = NULL,
                        min_n = NULL,
                        strata_min = 0,
                        verbose = 0) {

  # 1. Is the R package installed?
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "The 'cmdstanr' package is required but not installed.\n",
      "Install it with:\n",
      "  install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))",
      call. = FALSE
    )
  }

  # 2. Is CmdStan installed?
  ver <- cmdstanr::cmdstan_version(error_on_NA = FALSE)

  if (is.na(ver)) {
    stop(
      "CmdStan is not installed.\n",
      "Install it with:\n",
      "  cmdstanr::install_cmdstan()",
      call. = FALSE
    )
  }

  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcomes_tbl))

  #
  stopifnot(stan_type %in% c('laplace', 'mcmc'))
  cat(" STAN TYPE =", stan_type, '\n')

  #
  stopifnot(use_spatial_model %in% c('none', 'bym2', 'leroux'))
  cat(" SPATIAL MODEL =", use_spatial_model, '\n')

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
  orig_basis  <- vector("list", n_geos);
  orig_coef   <- vector("list", n_geos);
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

    local_cp <- local_cp$`_`$out[[1]]

    blup_cp <- get_centered_cp(argvar = local_cp$argvar,
                               xcoef = local_cp$coef,
                               xvcov = local_cp$vcov,
                               global_cen = global_cen,
                               cen = local_cp$cen,
                               this_exp = this_exp,
                               x_b = x_b)

    # get cb and outcomes lists
    orig_basis[[i]]  <- local_cp$orig_basis
    orig_coef[[i]]   <- local_cp$orig_coef
    coef_list[[i]]   <- local_cp$coef
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
  k_l <- lapply(orig_basis, \(x) ncol(as.matrix(x)))
  k_l <- unique(do.call(c, k_l))
  stopifnot(length(k_l) == 1)
  K = as.integer(k_l)

  # include the intercept
  X = array(dim = c(dim(orig_basis[[1]]), J))
  for(j in 1:J) X[,,j] = as.matrix(orig_basis[[j]])

  # outcome in J regions
  y = array(dim = c(n_l, J))
  for(j in 1:J) y[,j] = outc_list[[j]][, get(outcome_columns$outcome)]

  # create S matrix
  S <- getSmat(factor(outc_list[[1]]$strata))
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

  #' ////////////////////////////////////////////////////////////////////////////
  #' ============================================================================
  #' Get graph, has to make sure the shp is in the right order!
  #' ============================================================================
  #' #' /////////////////////////////////////////////////////////////////////////

  exp_geo_unit_col <- exposure_columns$geo_unit

  geo_units_in_order <- unique(exposure_matrix[, get(exp_geo_unit_col)])

  stopifnot(exp_geo_unit_col %in% names(shp_sf))

  rr <- sapply(geo_units_in_order, \(x) which(shp_sf[[exp_geo_unit_col]] == x))
  stopifnot(length(rr) == J)

  shp_sf_safe <- shp_sf[rr, ]

  nb_subset = spdep::poly2nb(shp_sf_safe);

  nbs = nb2graph(nb_subset);

  node1 = nbs$node1;
  node2 = nbs$node2;
  N_edges = nbs$N_edges;
  scaling_factor = scale_nb_components(nb_subset)[1];

  SW <- getSW(shp = shp_sf_safe, ni = 1, include_self = F)

  if(verbose > 0) {
    cat("\n")
  }

  #' ////////////////////////////////////////////////////////////////////////////
  #' ============================================================================
  #' Run STAN
  #' ============================================================================
  #' #' /////////////////////////////////////////////////////////////////////////

  if(verbose > 0) {
    cat("-- run STAN\n")
  }

  stan_data <- list(
    J = J,
    Jmat = SW,
    N = N,
    K = K,
    X = X,
    y = y,
    n_strata = n_strata,
    max_in_strata = max_in_strata,
    S_condensed = S_condensed,
    stratum_id = stratum_id,
    N_edges = N_edges,
    node1 = node1,
    node2 = node2,
    scaling_factor = scaling_factor
  )


  init_fun <- function() {
    list(
      beta = as.matrix(do.call(cbind, orig_coef)),
      bym2_sigma = 0.5,
      bym2_rho = 0.5,
      bym2_theta = rep(0, J),
      bym2_phi = rep(0, J),
      mu = rep(0, K),
      sigma = rep(1.0, K),
      q = 0.5,
      bstar = matrix(0, K, J)
    )
  }

  # from here: https://discourse.mc-stan.org/t/r-package-using-cmdstanr/32758
  # says you can use instantiate
  if(use_spatial_model == 'bym2') {

    stan_file <- system.file("stan", "bym2_condPois.stan",
                             package = "cityHeatHealth")

    mod <- cmdstanr::cmdstan_model(stan_file)

  }

  if(use_spatial_model == 'leroux') {

    stan_file <- system.file("stan", "leroux_condPois.stan",
                             package = "cityHeatHealth")

    mod <- cmdstanr::cmdstan_model(stan_file)
  }

  if(use_spatial_model == 'none') {

    stan_file <- system.file("stan", "condPois.stan",
                             package = "cityHeatHealth")

    mod <- cmdstanr::cmdstan_model(stan_file)

  }


  # ****************
  # **** RUN STAN ***
  # ****************
  out_df <- tryCatch({

    if(stan_type == 'mcmc') {

      default_stan_opts <- list(
        chains = 2,
        parallel_chains = 2,
        refresh = 100,
        init = init_fun
      )

      stan_opts <- modifyList(default_stan_opts, stan_opts %||% list())

      cat(" ...mcmc... \n")
      fit_mcmc <- do.call(
        mod$sample,
        c(list(data = stan_data), stan_opts)
      )

      cat(" ...mcmc draws... \n")
      mcmc_array <- fit_mcmc$draws()

      stan_summary <- fit_mcmc$summary()

      # outdf
      oo <- posterior::as_draws_df(mcmc_array)
    }


    if(stan_type == 'laplace') {

      default_stan_opts <- list(
        init = init_fun,
        jacobian = TRUE,
        iter = 1e4
      )

      stan_opts <- modifyList(default_stan_opts, stan_opts %||% list())

      cat(" ...laplace optimize... \n")
      fit_mode <- do.call(
        mod$optimize,
        c(list(data = stan_data), stan_opts)
      )

      cat(" ...laplace sample... \n")
      fit_laplace <- mod$laplace(data = stan_data,
                                 mode = fit_mode)

      cat(" ...laplace draws... \n")
      laplace_array <- fit_laplace$draws()

      stan_summary <- fit_laplace$summary()

      # outdf
      oo <- posterior::as_draws_df(laplace_array)
    }

    if(stan_type == 'variational') {

      default_stan_opts <- list(
        init = init_fun,
      )

      stan_opts <- modifyList(default_stan_opts, stan_opts %||% list())

      cat(" ...variational... \n")
      fit_var <- do.call(
        mod$variational,
        c(list(data = stan_data), stan_opts)
      )

      cat(" ...variational draws... \n")
      var_array <- fit_var$draws()

      stan_summary <- fit_var$summary()

      #out df
      oo <- posterior::as_draws_df(var_array)
    }

    if(stan_type == 'pathfinder') {

      default_stan_opts <- list(
        init = init_fun,
      )

      stan_opts <- modifyList(default_stan_opts, stan_opts %||% list())

      cat(" ...pathfinder... \n")
      fit_path <- do.call(
        mod$pathfinder,
        c(list(data = stan_data), stan_opts)
      )

      cat(" ...pathfinder draws... \n")
      path_array <- fit_path$draws()

      stan_summary <- fit_path$summary()

      #outdf
      oo <- posterior::as_draws_df(path_array)
    }

    oo
  },
  error = function(e) {
    # Exit the *entire* function immediately
    warning("some STAN error, exiting early with STAN objects so manual testing can occur")
    return(list(stan_data = stan_data, stan_opts = stan_opts))
  }
  )

  # ****************
  # **** COEF ******
  # ****************

  # and then finally the draws etc
  setDT(out_df)
  beta_reg_all <- out_df[, .SD, .SDcols = patterns("^beta")]
  beta_reg_all <- apply(beta_reg_all, 2, mean)
  beta_mat <- matrix(beta_reg_all, nrow = K, ncol = J, byrow = F)
  colnames(beta_mat) = geo_units_in_order

  # ***********************
  # **** VCOV and CR ******
  # ***********************

  # so this gives you the coef
  # but then you also need to get VCOV

  cr_coef <- vector("list", n_geos)
  cr_vcov <- vector("list", n_geos)
  cen_list_updated <- vector("list", n_geos)

  outcomes_col <- attributes(outcomes_tbl)$column_mapping$outcome

  geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit
  geo_unit_grp_col <- attributes(outcomes_tbl)$column_mapping$geo_unit_grp

  for(ni in 1:n_geos) {

    # get the name, which you know exists in both datasets
    this_geo <- unique_geos[ni, get(outcome_columns$geo_unit)]

    if(verbose > 1) {
      cat(this_geo, '\t')
    }

    # this cities exposure matrix
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_columns$exposure)]
    exp_mean = mean(this_exp)
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]

    # dispersion param
    dispersion <- calc_dispersion(y = single_outcomes_tbl[, get(outcomes_col)],
                                  X = orig_basis[[ni]],
                                  stratum_vector = single_outcomes_tbl$strata,
                                  beta = beta_mat[,ni])

    # VCOV
    vcov_beta <- calc_vcov(y = single_outcomes_tbl[, get(outcomes_col)],
                           X = orig_basis[[ni]],
                           stratum_vector = single_outcomes_tbl$strata,
                           beta = beta_mat[,ni])

    vcov_beta <- vcov_beta * dispersion

    # now -- get CROSSREDUCE Object !

    # the crossreduce coefficients are not affected by the centering point
    # but it does make a message if you dont put something there
    # so i center on the min to avoid that fate
    cb = orig_basis[[ni]]

    cp <- crosspred(basis = cb,
                    coef = beta_mat[,ni],
                    vcov = vcov_beta,
                    model.link = "log",
                    cen = exp_mean,
                    by = 0.1)

    if(!is.null(global_cen)) {
      cen = global_cen
    } else {
      cen = cp$predvar[which.min(cp$allRRfit)]
    }

    # now apply to cr and export
    cr <- crossreduce(basis = orig_basis[[ni]],
                      coef = beta_mat[,ni],
                      vcov = vcov_beta,
                      model.link = "log",
                      cen = cen,
                      by = 0.1)

    cr_coef[[ni]] <- coef(cr)
    cr_vcov[[ni]] <- vcov(cr)
    cen_list_updated[[ni]] = cen

  }

  if(verbose > 0) {
    cat("\n")
  }

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

  out   <- vector("list", n_geos);

  # loop through geos
  for(i in 1:n_geos) {

    #
    this_geo <- unique_geos[i, get(outcome_columns$geo_unit)]
    this_geo_grp <- unique_geos[i, get(outcome_columns$geo_unit_grp)]

    # get x
    rr <- exposure_matrix[, get(exp_geo_unit_col)] == this_geo
    single_exposure_matrix = exposure_matrix[rr, ,drop = FALSE]
    this_exp <- single_exposure_matrix[, get(exposure_col)]
    x_b <- c(floor(min(this_exp)), ceiling(max(this_exp)))

    # get centered basis
    blup_cp <- get_centered_cp(argvar = argvar_list[[i]],
                               xcoef = cr_coef[[i]],
                               xvcov = cr_vcov[[i]],
                               global_cen = global_cen,
                               cen = cen_list_updated[[i]],
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
    names(RRdf)[1] <- outcome_columns$geo_unit
    names(RRdf)[2] <- outcome_columns$geo_unit_grp
    names(RRdf)[3] <- exposure_col
    setDT(RRdf)

    ## attach outcomes vector
    rr <- outcomes_tbl[, get(out_geo_unit_col)] == this_geo
    single_outcomes_tbl = outcomes_tbl[rr, ,drop = FALSE]
    outcomes_vec = single_outcomes_tbl[, get(outcomes_col)]
    stopifnot(sum(outcomes_vec) == sum(outc_list[[i]][, get(outcomes_col)]))

    #
    out[[i]] <- list(
      geo_unit = this_geo,
      basis_cen = blup_cp$basis_cen,
      exposure_col = exposure_col,
      this_exp = this_exp,
      cen = blup_cp$cp$cen,
      global_cen = global_cen,
      outcomes = outc_list[[i]],
      coef = cr_coef[[i]],
      vcov = cr_vcov[[i]],
      RRdf = RRdf
    )

  }

  # set names
  names(out) <- unique_geos[, get(outcome_columns$geo_unit)]

  if(verbose > 0) {
    cat("\n")
  }

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' OUTPUT
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # yes, a list in a list, this will make sense later
  # aka, in the recursive call to this function that happens above
  # you could modify this to also output the mixmeta object
  # but not clear that you need that
  outlist = list(list(out = out, beta_mat = beta_mat,
                      stan_summary = stan_summary))
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
#' forest_plot.condPois_sb
#'
#' @param x an object of class condPois_sb
#' @param exposure_val exposure value at which to plot
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
forest_plot.condPois_sb <- function(x, exposure_val) {

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
#' spatial_plot.condPois_sb
#'
#' @param x an object of class condPois_sb
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param exposure_val exposure value at which to plot
#' @importFrom ggplot2 ggplot
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_sb <- function(x, shp, exposure_val,
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
#' spatial_plot.condPois_sb_list
#'
#' @param x an object of class condPois_sb_list
#' @param shp an sf shapefile with an appropriate column at which to join
#' @param exposure_val exposure value at which to plot
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @returns
#' @export
#'
#' @examples
spatial_plot.condPois_sb_list <- function(x, shp, exposure_val) {

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
#' forest_plot.condPois_sb_list
#'
#' @param x an object of class condPois_sb_list
#' @param exposure_val exposure value at which to plot
#' @returns
#' @export
#'
#' @examples
forest_plot.condPois_sb_list <- function(x, exposure_val) {

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
      scale_x_continuous(transform = 'log') +
      scale_color_viridis_d() +
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

#'@export
#' getRR.condPois_sb
#'
#' @param x
#' @importFrom data.table setDT
#' @returns
#' @export
#'
#' @examples
getRR.condPois_sb <- function(x) {
  oo <- do.call(rbind, lapply(x$`_`$out, \(obj) obj$RRdf))
  oo$model_class = class(x)
  setDT(oo)
  return(oo)
}

#'@export
#' getRR.condPois_sb_list
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
getRR.condPois_sb_list <- function(x) {

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
