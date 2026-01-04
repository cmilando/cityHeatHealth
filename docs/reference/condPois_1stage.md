# Run a conditional poisson model for a single geographic unit

Run a conditional poisson model for a single geographic unit

## Usage

``` r
condPois_1stage(
  exposure_matrix,
  outcomes_tbl,
  argvar = NULL,
  arglag = NULL,
  maxlag = NULL,
  min_n = NULL,
  strata_min = 0,
  global_cen = NULL,
  multi_zone = FALSE
)
```

## Arguments

- exposure_matrix:

  a matrix of exposures, with columns for lag, usually created by
  `make_exposure_matrix`

- outcomes_tbl:

  a data.table of outcomes, created by `make_outcome_table`

- argvar:

  a list containing the `argvar` components for the `crossbasis`

- arglag:

  a list containing the `arglag` components for the `crossbasis`

- maxlag:

  an integer of the maximum lag

- min_n:

  an integer describing the minimum number of cases for a single region

- strata_min:

  an integer describing the minimum number of cases for a single strata

- global_cen:

  global centering point

- multi_zone:

  are multiple strata being used.

## Examples
