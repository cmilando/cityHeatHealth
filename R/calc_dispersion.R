#' Calculate the dispersion parameter for a quasi-poisson model
#'
#' Converted from STATA from armstrong 2014 supplement with guidance from ChatGPT
#'
#' @param y a vector of outcomes
#' @param X a matrix of predictors, typically the crossbasis output
#' @param beta a vector of coefficients
#' @param stratum_vector a vector describing the stratum
#'
#' @returns
#' @export
#'
#' @examples
calc_dispersion <- function(y, X, beta, stratum_vector) {

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
