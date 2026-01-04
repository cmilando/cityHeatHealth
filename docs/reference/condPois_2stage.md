# Run a conditional poisson model for across many geographic units

MixMeta nope its this:
https://www.sciencedirect.com/science/article/pii/S0140673614621140?via%3Dihub
avgtmax and range and random effect for county unclear, Gasp does both
here
https://github.com/gasparrini/RiskExtrapolation/blob/3caa5c46a748d6a789021554f435e53f9960d4c1/24_MetaRegression.R
and it goes counrty/city
https://github.com/gasparrini/MCC-SO2/blob/1a158b31bf83f3bbca3449b553a6570280306299/Rcode/06.nonlin.R#L47
nlmeta \<- mixmeta(nlcoefall, nlvcovall, random=~1\|country/city,
data=cities, control=list(showiter=T, igls.inititer=10)) basically, did
Gasp 2015 (temp_mean, temp_iqr, random effect for county) and included
the hierarchicacy structure for ids within countries
https://github.com/gasparrini/2015_gasparrini_EHP_Rcodedata
https://github.com/gasparrini/2019_sera_StatMed_Rcode/blob/master/03.MultilevelMA.R

## Usage

``` r
condPois_2stage(
  exposure_matrix,
  outcomes_tbl,
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

## Details

update as of 2025.12.21, its this
https://github.com/gasparrini/Extended2stage/blob/cf82fd2f0d616210e8aa6911d2fb49a732f65795/02.multilevel.R#L68

## Examples
