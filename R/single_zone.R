#' Run a conditional poisson model for a single zone

#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr join_by
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom gnm gnm
#' @importFrom lubridate month
#' @importFrom lubridate year
#' @importFrom lubridate wday
#' @param exposure_matrix
#' @param outcomes
#' @param argvar
#' @param arglag
#' @param maxlag
#'
#' @returns
#' @export
#'
#' @examples
single_zone <- function(exposure_matrix, outcomes,
                        argvar = NULL, arglag = NULL, maxlag = 5) {

  ## ******************
  ## TODO
  ## - I want to have defaults for argvar and arglag
  ## - tmaxC has to not be hard-coded, you'll figure this tou
  ## - TOWN20 has to be renamed, and flexible to what the strata are
  ## - and keep in mind you'll be working on different resolutions
  ## - give a warning that the strata are fixed at day, month, year
  ##   and area unit. this is the default unless you really want to
  ##   do things differently
  ## ******************

  if(is.null(argvar)) {
    x_knots = quantile(exposure_matrix$tmax_C, probs = c(0.5, 0.9))
    argvar <- list(fun = 'ns', knots = x_knots)
  }

  if(is.null(arglag)) {
    arglag <- list(fun = 'ns', knots = dlnm::logknots(maxlag, nk = 2))
  }

  ## keep only the columns you need
  x_mat <- exposure_matrix[, c('tmax_C', paste0('Templag',1:maxlag))]

  ## create the crossbasis
  cb <- crossbasis(x_mat, lag = maxlag, argvar = argvar, arglag = arglag)

  ## create the strata
  outcomes$dow   <- lubridate::wday(outcomes$date, label = T)
  outcomes$month <- lubridate::month(outcomes$date, label = T)
  outcomes$year  <- lubridate::year(outcomes$date)

  outcomes$strata <- paste0(outcomes$TOWN20, ":",
                            outcomes$year, ":",
                            outcomes$month, ":",
                            outcomes$dow)

  ## get rid of 0 strata
  ## you make sure there are no empty strata
  outcomes_agg <- outcomes %>%
    group_by(strata) %>%
    summarize(
      .groups = 'keep',
      total_daily_deaths = sum(daily_deaths)
    ) %>%
    mutate(keep = ifelse(total_daily_deaths > 0, 1, 0))

  outcomes_comb <- left_join(outcomes, outcomes_agg, by = join_by(strata))

  ## ******************************************
  ## if using GNM, you get COEF and VCOV as part of the model objects
  m_sub <- gnm(daily_deaths ~ cb,
               data = outcomes_comb,
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
                  cen = mean(exposure_matrix$tmax_C),
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

#' Calculate the dispersion parameter for a quasi-poisson model
#'
#' @param y
#' @param X
#' @param beta
#' @param stratum_vector
#'
#' @returns
#' @export
#'
#' @examples
calc_dispersion <- function(y, X, beta, stratum_vector) {

  # converted from STATA from armstrong 2014.

  # get dims of X
  N = nrow(X)
  K = ncol(X)

  # replace stratum_vector with counts
  n_strata <- length(unique(stratum_vector))
  stratum_vector <- as.numeric(factor(boston_deaths$strata, labels = 1:n_strata))

  # sum observed and predicted counts per stratum
  sum_y_stratum    = rep(0, n_strata)
  sum_pred_stratum = rep(0, n_strata)
  xBeta_out = exp(X %*% beta)
  for (n in 1:N) {
      s = stratum_vector[n]
      sum_y_stratum[s] = sum_y_stratum[s] + y[n]
      sum_pred_stratum[s] = sum_pred_stratum[s] + xBeta_out[n]
  }

  # rescale predictions so that the sum of the predicted in each group
  # equals the sum of the observed in each group
  pred_rescaled <- rep(NA, N)
  for (n in 1:N) {
    s = stratum_vector[n]
    pred_rescaled[n] = xBeta_out[n] * sum_y_stratum[s] / sum_pred_stratum[s]
  }
  any(is.na(pred_rescaled))

  # compute pearson chi-squared
  pearson_x2 = 0
  for (n in 1:N) {
      pearson_x2 = pearson_x2 + (y[n] - pred_rescaled[n])^2 / pred_rescaled[n]
  }
  pearson_x2

  # calculate dispersion
  df_resid = N - K - n_strata
  dispersion = pearson_x2 / df_resid

  return(dispersion)

}

#' Calculate the variance-covariance matrix of a poisson model
#'
#' https://link.springer.com/article/10.1007/s11222-005-4069-4
#' https://statomics.github.io/SGA2019/assets/poissonIRWLS-implemented.html#variance-covariance-matrix-of-the-model-parameters
#'
#' and extended to multi-nomial case using ChatGPT, and confirmed it works.
#'
#' @param X
#' @param beta
#' @param stratum_vector
#' @importFrom MASS ginv
#' @returns
#' @export
#'
#' @examples
calc_vcov <- function(y, X, beta, stratum_vector) {
  # X      : n x p design matrix (no stratum intercepts)
  # beta   : length-p coefficient vector
  # strata : length-n vector defining conditioning strata

  # --- basic checks ----------------------------------------------------------
  stopifnot(
    is.matrix(X),
    length(y) == nrow(X),
    length(beta) == ncol(X),
    length(stratum_vector) == nrow(X)
  )

  # Recode strata as consecutive integers 1, ..., S
  n_strata <- length(unique(stratum_vector))
  stratum_vector <- as.numeric(factor(stratum_vector, labels = 1:n_strata))

  # Initialize the Fisher information matrix for beta
  # Dimension: p x p, where p = number of regression coefficients
  # We will accumulate information across strata
  p <- ncol(X)
  I <- matrix(0, p, p)

  # Loop over strata
  # Each iteration adds the Fisher information contribution
  # from one multinomial likelihood (one conditional Poisson stratum)
  for (s in seq_len(n_strata)) {

    # Identify which observations belong to stratum s
    # These observations share a fixed total count
    idx <- which(stratum_vector == s)

    # Skip degenerate strata
    if (length(idx) <= 1) next

    # Extract the design matrix rows for stratum s
    # Dimension: n_s x p
    # drop = FALSE ensures Xs stays a matrix even if n_s = 1
    Xs <- X[idx, , drop = FALSE]

    # Compute the linear predictor for stratum s
    # eta_is = x_is^T beta
    eta <- as.vector(Xs %*% beta)

    # Convert linear predictors to unnormalized intensities
    # These are proportional to multinomial probabilities
    ps <- exp(eta)

    # Normalize to obtain multinomial probabilities
    # p_is = exp(eta_is) / sum_j exp(eta_js)
    # These probabilities sum to 1 within each stratum
    ps <- ps / sum(ps)

    # Compute the observed total count in stratum s
    # This is the conditioning value in the conditional Poisson
    Ns <- sum(y[idx])

    # Construct the multinomial covariance (weight) matrix
    # diag(ps) gives Var(Y_is | Ns)
    # tcrossprod(ps) = ps %*% t(ps) gives Cov(Y_is, Y_js | Ns)
    #
    # Ws = diag(ps) - ps ps^T
    #
    # This matrix:
    # - encodes negative correlation within strata
    # - has rank (n_s - 1)
    # - removes one degree of freedom due to conditioning
    Ws <- diag(ps) - tcrossprod(ps)

    # Add the Fisher information contribution from stratum s
    #
    # Multinomial Fisher information:
    # I_s = Ns * X_s^T Ws X_s
    #
    # Ns scales the information by the total count in the stratum
    I <- I + Ns * crossprod(Xs, Ws %*% Xs)
  }

  # Invert the total Fisher information matrix
  # This yields the variance–covariance matrix of beta
  # under the conditional Poisson (multinomial) likelihood
  return(ginv(I))
}


