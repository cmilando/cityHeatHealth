# Run a conditional poisson model with Spatial Bayesian inference

https://mc-stan.org/docs/2_20/functions-reference/multinomial-distribution.html

## Usage

``` r
condPois_sb(
  exposure_matrix,
  outcomes_tbl,
  shp_sf,
  stan_type = "laplace",
  use_spatial_model = "none",
  stan_opts = NULL,
  global_cen = NULL,
  argvar = NULL,
  arglag = NULL,
  maxlag = NULL,
  min_n = NULL,
  strata_min = 0,
  verbose = 0
)
```

## Arguments

- exposure_matrix:

  a matrix of exposures, with columns for lag, usually created by
  `make_exposure_matrix`

- outcomes_tbl:

  a data.table of outcomes, created by `make_outcome_table`

- shp_sf:

  a shapefile object describing the geo_units

- stan_type:

  a choice of either laplace or sample

- stan_opts:

  a list of parameters to send to stan

- global_cen:

  a global centering point

- argvar:

  a list containing the `argvar` components for the `crossbasis`

- arglag:

  a list containing the `arglag` components for the `crossbasis`

- maxlag:

  an integer of the maximum lag

- min_n:

  an integer describing the minimum number of cases for a single region

- strata_min:

  minimum number of cases per strata

- verbose:

  whether to print geo_units as they are run

## Examples
