# Function to create the outcome table

Function to create the outcome table

## Usage

``` r
make_outcome_table(
  data,
  column_mapping,
  warm_season_months = 5:9,
  collapse_to = NULL,
  grp_level = FALSE
)
```

## Arguments

- column_mapping:

  a named list that indicates relevant columns in `data`. for the
  exposure data table, these need to be one of: c('date',
  "outcome",'factor, 'geo_unit', 'geo_unit_grp')

- warm_season_months:

  the warm season months for this region, default is Northern
  Hemisphere's May through September (5 through 9)

- collapse_to:

  which factors to collapse across

- grp_level:

  whether to summarize to the group level or not (default)

## Examples
