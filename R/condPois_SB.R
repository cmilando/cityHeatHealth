#' Run a conditional poisson model with Spatial Bayesian inference
#'
#' https://mc-stan.org/docs/2_20/functions-reference/multinomial-distribution.html
#'
#' @import data.table
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom gnm gnm
#' @param exposure_matrix a matrix of exposures, with columns for lag, usually created by `make_exposure_matrix`
#' @param outcomes a data.table of outcomes, created by `make_outcome_table`
#' @param shp a shapefile describing the network
#' @param argvar a list containing the `argvar` components for the `crossbasis`
#' @param arglag a list containing the `arglag` components for the `crossbasis`
#' @param maxlag an integer of the maximum lag
#' @param min_n an integer describing the minimum number of cases for the entire run
#'
#' @returns
#' @export
#'
#' @examples
condPoisSB <- function(exposure_matrix, outcomes_tbl, shp,
                     argvar = NULL, arglag = NULL, maxlag = NULL,
                     min_n = 50) {

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' VALIDATIONS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  ## Check 1 -- that both inputs are the right class of variables
  stopifnot("exposure" %in% class(exposure_matrix))
  stopifnot("outcome" %in% class(outcomes_tbl))

  ## Check 2
  ## probably should make sure that exposure_matrix and outcome_tbl
  ## are the same size, at least
  exp_date_col <- attributes(exposure_matrix)$column_mapping$date
  outcome_date_col <- attributes(outcomes_tbl)$column_mapping$date

  exp_geo_unit_col <- attributes(exposure_matrix)$column_mapping$geo_unit
  outcome_geo_unit_col <- attributes(outcomes_tbl)$column_mapping$geo_unit

  setorderv(
    exposure_matrix,
    c(exp_date_col, exp_geo_unit_col)
  )

  setorderv(
    outcomes_tbl,
    c(outcome_date_col, outcome_geo_unit_col)
  )

  stopifnot(dim(exposure_matrix)[1] == dim(outcomes_tbl)[1])
  stopifnot(identical(exposure_matrix[, get(exp_date_col)],
                      outcomes_tbl[, get(outcome_date_col)]))

  # CHECK 4 geo_unit is the same for both"
  stopifnot(all(outcomes_tbl[, get(outcome_geo_unit_col)] %in%
                  exposure_matrix[, get(exp_geo_unit_col)]))

  # CHECK 5
  if("factor" %in% names(attributes(outcomes_tbl)$column_mapping)) {
    stop("if outcome has a factor, thats a problem")
  }

  ## CHECK 6 - minN
  outcome_col <- attributes(outcomes_tbl)$column_mapping$outcome
  stopifnot(sum(outcomes_tbl[, get(outcome_col)]) >= min_n)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' CREATE CROSSBASIS
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  stop("this has to be by region")

  exposure_col <- attributes(exposure_matrix)$column_mapping$exposure

  if(is.null(maxlag)) {
    maxlag = 5
  }

  if(is.null(argvar)) {
    x_knots = quantile(exposure_matrix[, get(exposure_col)], probs = c(0.5, 0.9))
    argvar <- list(fun = 'ns', knots = x_knots)
  }

  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = dlnm::logknots(maxlag, nk = 2))
  }

  ## get the columns you need
  xcols <- c(exposure_col, paste0('explag',1:maxlag))
  x_mat <- exposure_matrix[, ..xcols]

  ## if you are safe to proceed, make the x_mat
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  #' //////////////////////////////////////////////////////////////////////////
  #' ==========================================================================
  #' Prepare STAN inputs
  #' ==========================================================================
  #' //////////////////////////////////////////////////////////////////////////

  # n regions
  J = as.integer(length(data_l))

  # nrows
  N = as.integer(nrow(data_l[[1]]))

  # beta values, withouth the intercept
  K = as.integer(ncol(list_X[[1]]))

  # include the intercept
  X = array(dim = c(dim(list_X[[1]]), J))
  for(j in 1:J) X[,,j] = as.matrix(list_X[[j]])

  # outcome in two regions
  y = array(dim = c(nrow(list_X[[1]]), J))
  for(j in 1:J) y[,j] = data_l[[j]]$mort

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
  S <- getSmat(factor(data_l[[1]]$strata), include_self = T)
  head(S)

  # get strata vars
  n_strata <- as.integer(length(unique(data_l[[1]]$strata))) # 72, cool
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
  stratum_id = as.integer(factor(data_l[[1]]$strata))
  stratum_id

  ## *** THIS BECOMES THE J MATRIX
  # ok now do the same as include self but with J matrix
  # start simple and you can generalize later
  # this is an adjacency matrix that DOES NOT include itself
  shapefile_bcn <- read_sf("input/shapefile_bcn.shp")

  # Generate a list of spatial structure from the shapefile for use in WinBUGS.
  list_neig <- nb2listw(poly2nb(shapefile_bcn))

  neighbors <- lapply(list_neig$neighbours,c)
  Jmat <- matrix(0, nrow = J, ncol = J)
  if(j > 1) {
    for(j in 1:J) {
      Jmat[j, neighbors[[j]]] <- 1
    }
  }
  head(Jmat)
  ## ******



  #' ////////////////////////////////////////////////////////////////////////////
  #' ============================================================================
  #' Run STAN
  #' ============================================================================
  #' #' /////////////////////////////////////////////////////////////////////////

  grainsize = as.integer(4)

  stan_data <- list(
    J = J,
    Jmat = Jmat,
    N = N,
    K = K,
    X = X,
    y = y,
    S = S,
    n_strata = n_strata,
    max_in_strata = max_in_strata,
    S_condensed = S_condensed,
    stratum_id = stratum_id,
    grainsize = grainsize
  )

  # Set path to model
  # always compile with threads, even if you only use 1
  stan_model <- cmdstan_model("SB_CondPoisson.stan")

  # # probably takes ~ 6 hours or so
  out2 <- stan_model$sample(
    data = stan_data,
    chains = 2,
    iter_warmup = 3000,
    iter_sampling = 3000,
    parallel_chains = 2,
    threads_per_chain = 1,
    refresh = 10,
    #adapt_delta = 0.8,
    max_treedepth = 7 # .... ?
  )
  #
  # ##
  # draws1_array <- out1$draws()
  #
  # # Convert to data.frame (flattened, easier to use like extract())
  # draws1_df <- posterior::as_draws_df(draws1_array)
  # head(draws_df)

  ##
  # out2 <- stan_model$variational(
  #   data = stan_data,
  #   iter = 100000,
  #   threads = 4,
  #   refresh = 500
  # )

  ## if using GNM, you get COEF and VCOV as part of the model objects
  m_sub <- gnm(daily_deaths ~ cb,
               data = outcomes_tbl,
               family = quasipoisson,
               eliminate = factor(strata),
               subset = keep == 1)

  m_coef <- coef(m_sub)
  m_vcov <- vcov(m_sub)
  ## ******************************************

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = mean(exposure_matrix[, get(exposure_col)]),
                  by = 0.1)

  cen = cp$predvar[which.min(cp$allRRfit)]

  cp <- crosspred(cb,
                  coef = m_coef,
                  vcov = m_vcov,
                  model.link = "log",
                  cen = cen,
                  by = 0.1)

  return(cp)


}
