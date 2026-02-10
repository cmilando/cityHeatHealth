# QA_testing

``` r

library(cityHeatHealth)
```

``` r



# *****
#slow way to update
# unloadNamespace("cityHeatHealth")
# remotes::install_github("cmilando/cityHeatHealth")

# once you update your R version you can just do
#load locally with devtools::load_all()
#if you want to actually make changes to the files.
# *****

library(data.table)

data <- readRDS("../data-raw/weekly_mort_updated_nov17.RDS")
setDT(data)

# just get a specific subset
data <- subset(data, region %in% 'Southeast')

# Select only the needed columns
data <- data[, .(
  X_week_end_date,
  X_deaths,
  mean_temp_C,
  occurrence_county_code,
  region
)]

# Order by county code and week end date
setorder(data, occurrence_county_code, X_week_end_date)
```

Exposure data

``` r

exposure_columns <- list(
  "date" = "X_week_end_date",
  "exposure" = "mean_temp_C",
  "geo_unit" = "occurrence_county_code",
  "geo_unit_grp" = "region"
)

exp_data <- data[, .(
  X_week_end_date,
  mean_temp_C,
  occurrence_county_code,
  region
)]

exposure_matrix <- make_exposure_matrix(
  exp_data, exposure_columns, dt_by = 'week')

head(exposure_matrix)
#>    X_week_end_date mean_temp_C occurrence_county_code    region
#>             <IDat>       <num>                 <char>    <char>
#> 1:      2018-05-05    25.03714                  01073 Southeast
#> 2:      2018-05-12    28.61286                  01073 Southeast
#> 3:      2018-05-19    29.71857                  01073 Southeast
#> 4:      2018-05-26    28.69714                  01073 Southeast
#> 5:      2018-06-02    28.35000                  01073 Southeast
#> 6:      2018-06-09    30.47857                  01073 Southeast
#>                     strata     match_strata  explag1  explag2  explag3  explag4
#>                     <char>           <char>    <num>    <num>    <num>    <num>
#> 1: 01073:yr2018:mn05:dow07 01073:2018-05-05 20.02857 19.90857 21.19714 20.80714
#> 2: 01073:yr2018:mn05:dow07 01073:2018-05-12 25.03714 20.02857 19.90857 21.19714
#> 3: 01073:yr2018:mn05:dow07 01073:2018-05-19 28.61286 25.03714 20.02857 19.90857
#> 4: 01073:yr2018:mn05:dow07 01073:2018-05-26 29.71857 28.61286 25.03714 20.02857
#> 5: 01073:yr2018:mn06:dow07 01073:2018-06-02 28.69714 29.71857 28.61286 25.03714
#> 6: 01073:yr2018:mn06:dow07 01073:2018-06-09 28.35000 28.69714 29.71857 28.61286
#>     explag5
#>       <num>
#> 1: 19.44286
#> 2: 20.80714
#> 3: 21.19714
#> 4: 19.90857
#> 5: 20.02857
#> 6: 25.03714
```

Outcome

``` r

outcome_columns <- list(
  "date" = "X_week_end_date",
  "outcome" = "X_deaths",
  # "factor" = 'age_grp',
  # "factor" = 'sex',
  "geo_unit" = "occurrence_county_code",
  "geo_unit_grp" = "region"
)

out_data <- data[, .(
  X_week_end_date,
  X_deaths,
  occurrence_county_code,
  region
)]

outcomes_tbl <- make_outcome_table(out_data, outcome_columns)
#> Some NA values in outcome xgrid were removed
outcomes_tbl
#>       X_week_end_date occurrence_county_code    region X_deaths
#>                <IDat>                 <char>    <char>    <int>
#>    1:      2018-05-05                  01073 Southeast      205
#>    2:      2018-05-12                  01073 Southeast      186
#>    3:      2018-05-19                  01073 Southeast      174
#>    4:      2018-05-26                  01073 Southeast      176
#>    5:      2018-06-02                  01073 Southeast      167
#>   ---                                                          
#> 9297:      2024-08-31                  51810 Southeast       58
#> 9298:      2024-09-07                  51810 Southeast       72
#> 9299:      2024-09-14                  51810 Southeast       54
#> 9300:      2024-09-21                  51810 Southeast       58
#> 9301:      2024-09-28                  51810 Southeast       59
#>                        strata strata_total     match_strata
#>                        <char>        <num>           <char>
#>    1: 01073:yr2018:mn05:dow07          741 01073:2018-05-05
#>    2: 01073:yr2018:mn05:dow07          741 01073:2018-05-12
#>    3: 01073:yr2018:mn05:dow07          741 01073:2018-05-19
#>    4: 01073:yr2018:mn05:dow07          741 01073:2018-05-26
#>    5: 01073:yr2018:mn06:dow07          893 01073:2018-06-02
#>   ---                                                      
#> 9297: 51810:yr2024:mn08:dow07          312 51810:2024-08-31
#> 9298: 51810:yr2024:mn09:dow07          243 51810:2024-09-07
#> 9299: 51810:yr2024:mn09:dow07          243 51810:2024-09-14
#> 9300: 51810:yr2024:mn09:dow07          243 51810:2024-09-21
#> 9301: 51810:yr2024:mn09:dow07          243 51810:2024-09-28
```

``` r

model <- condPois_1stage(exposure_matrix, outcomes_tbl, 
                         multi_zone = T)
#> Warning in value[[3L]](cond): a check of making the centered basis for a
#> geo-unit 37021 did not pass. this likely means that the knots for the overall
#> basis are outside the range of exposures in this geographic unit. Consider
#> adjusting either the geo-units you are passing in, or the exposure variable
#> (e.g., switching from absolute to relative measures)
```

``` r

plot(model, title = 'Southeast')
```

![](QA_testing_files/figure-html/plot%20v1-1.png)

``` r

model2 <- condPois_2stage(exposure_matrix, outcomes_tbl,
                          verbose = 1)
#> -- validation passed
#> -- stage 1
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12009. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12031. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12057. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12073. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12099. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 12109. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 22051. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit 45019. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> 
#> -- mixmeta
#> formula: ~ 1 | region/occurrence_county_code 
#> -- stage 2
```

\`\`\`
