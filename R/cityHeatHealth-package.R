#' The 'cityHeatHealth' package.
#'
#' @description Tools for estimating heat-related health impacts using daily exposure and
#' outcome data at small spatial scales. The package supports one-stage
#' conditional Poisson models, two-stage meta-analytic designs, and spatial
#' Bayesian approaches that borrow strength across neighboring geographic
#' units when case counts are small. Methods include distributed lag
#' non-linear models and attributable number calculations, with workflows
#' designed for messy real-world epidemiologic data.
#'
#' @name cityHeatHealth-package
#' @aliases cityHeatHealth
#'
#' @import methods
#' @import Matrix
#' @import data.table
#'
#' @import ggplot2
#' @importFrom scales number
#'
#' @import sf
#' @import spdep
#'
#' @importFrom INLA inla.qinv
#' @importFrom patchwork wrap_plots
#' @importFrom MASS ginv
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom dlnm crossbasis
#' @importFrom dlnm crosspred
#' @importFrom dlnm crossreduce
#' @importFrom dlnm logknots
#' @importFrom dlnm onebasis
#' @importFrom gnm gnm
#' @importFrom lubridate days_in_month
#' @importFrom lubridate make_date
#' @importFrom mixmeta blup
#' @importFrom mixmeta mixmeta
#' @importFrom tidyr expand_grid
#' @importFrom zoo na.approx
#'
#' @references
#'
"_PACKAGE"
NULL
