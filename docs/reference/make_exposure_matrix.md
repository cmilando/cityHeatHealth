# Function to clean and prepare the exposure data matrix

Function to clean and prepare the exposure data matrix

## Usage

``` r
make_exposure_matrix(
  data,
  column_mapping,
  months_subset = 5:9,
  dt_by = "day",
  maxgap = 5,
  maxlag = 5,
  grp_level = FALSE
)
```

## Arguments

- data:

  a dataset of exposures

- column_mapping:

  a named list that indicates relevant columns in `data`. for the
  exposure data table, these need to be one of: c('date', "exposure",
  'geo_unit', 'geo_unit_grp')

- months_subset:

  the warm season months for this region, default is Northern
  Hemisphere's May through September (5 through 9)

- dt_by:

  is it daily data, or weekly or ...

- maxlag:

  the number of lags for the exposure variable (default is 5)

- grp_level:

  whether to summarize to the group level or not (default)

- maxgaps:

  the maximum allowable missing exposure data gap, to be passed to
  zoo::na.approx (default is 5)

## Value

a data.table of class("exposure")

## Examples
