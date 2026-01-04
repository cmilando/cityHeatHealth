# Calculate attributable number and attributable rates

Calculate attributable number and attributable rates

## Usage

``` r
calc_AN(
  model,
  outcomes_tbl,
  pop_data,
  agg_type,
  join_cols,
  nsim = 300,
  verbose = 0
)
```

## Arguments

- model:

  a model object of class `condPois_2stage`, `condPois_1stage`, or
  `condPois_bayes`

- outcomes_tbl:

  a table of outcomes, of class `outcomes`

- pop_data:

  population data

- agg_type:

  what is the spatial resolution you are aggregating to

- join_cols:

  how should you join population data to the outcome table

- nsim:

  number of simulations required for calculation of empirical CI
  (default = 300)

- verbose:

  0 = no printing, 1 = headers, 2 = detailed

## Examples
