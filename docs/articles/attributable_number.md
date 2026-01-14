# Calculating attributable numbers and rates in \`cityHeatHealth\`

``` r

library(cityHeatHealth)
```

Estimating Attributable numbers (and rates of attributable numbers) are
a key way that we can translate relative risks into numbers that are
more tangible in public health settings. Below we provide easy
functionality to go from our model objects to estimates of attributable
numbers and rates.

### Setup

The first step of calculating attributable numbers is having a
population data estimate.

This varies a lot by place and dataset, so we don’t include
functionality for it (but an example of how this could be done can be
seen in
[`vignette("get_pop_estimates")`](http://chadmilando.com/cityHeatHealth/articles/get_pop_estimates.md)).

Assume you are starting with a dataset for the entire timeframe that
looks like this:

``` r

library(data.table)
data("ma_pop_data")
setDT(ma_pop_data)
ma_pop_data
#>               TOWN20 Female_0-17 Female_18-64 Female_65+ Male_0-17 Male_18-64
#>               <char>       <num>        <num>      <num>     <num>      <num>
#>   1:      BARNSTABLE        3899        15017       6014      4499      14035
#>   2:          BOURNE        1891         5751       3212      1489       5302
#>   3:        BREWSTER         634         2518       2007       833       2628
#>   4:         CHATHAM         163         1477       1759       480       1265
#>   5:          DENNIS         573         3792       3133       784       4101
#>  ---                                                                         
#> 347:   WEST BOYLSTON         619         2021       1107       604       2554
#> 348: WEST BROOKFIELD         343         1162        578       243       1002
#> 349:     WESTMINSTER         847         2371       1131       762       2028
#> 350:      WINCHENDON        1254         3318        711      1031       3134
#> 351:       WORCESTER       18779        67750      15995     21129      69365
#>      Male_65+
#>         <num>
#>   1:     5458
#>   2:     2810
#>   3:     1721
#>   4:     1463
#>   5:     2359
#>  ---         
#> 347:      790
#> 348:      495
#> 349:     1081
#> 350:      924
#> 351:    11173
```

Need to do some transformations:

- pivot longer
- variable clean

Note again, this processing will vary by application so this approach is
not prescriptive !

Pivot longer:

``` r

ma_pop_data_long <- melt(
  ma_pop_data,
  id.vars = "TOWN20",
  variable.name = "sex_age",
  value.name = "population"
)
```

Variable clean:

``` r

ma_pop_data_long$sex_age <- as.character(ma_pop_data_long$sex_age)
varnames <- strsplit(ma_pop_data_long$sex_age, "_", fixed = T)
varnames <- data.frame(do.call(rbind, varnames))
names(varnames) <- c('sex', 'age_grp')
rr <- which(varnames$sex == 'Female')
varnames$sex[rr] <- 'F'
rr <- which(varnames$sex == 'Male')
varnames$sex[rr] <- 'M'
ma_pop_data_long$sex = varnames$sex
ma_pop_data_long$age_grp = varnames$age_grp
ma_pop_data_long$sex_age <- NULL
ma_pop_data_long
#>                TOWN20 population    sex age_grp
#>                <char>      <num> <char>  <char>
#>    1:      BARNSTABLE       3899      F    0-17
#>    2:          BOURNE       1891      F    0-17
#>    3:        BREWSTER        634      F    0-17
#>    4:         CHATHAM        163      F    0-17
#>    5:          DENNIS        573      F    0-17
#>   ---                                          
#> 2102:   WEST BOYLSTON        790      M     65+
#> 2103: WEST BROOKFIELD        495      M     65+
#> 2104:     WESTMINSTER       1081      M     65+
#> 2105:      WINCHENDON        924      M     65+
#> 2106:       WORCESTER      11173      M     65+
```

We assume that these properties hold for the entire timeframe of our
analysis, but you could also make a version of this dataset with a
‘year’ column.

Now, quickly get a
[`condPois_1stage()`](http://chadmilando.com/cityHeatHealth/reference/condPois_1stage.md)
and
[`condPois_2stage()`](http://chadmilando.com/cityHeatHealth/reference/condPois_2stage.md)
objects to use in testing: exposures

``` r

library(data.table)
exposure_columns <- list(
  "date" = "date",
  "exposure" = "tmax_C",
  "geo_unit" = "TOWN20",
  "geo_unit_grp" = "COUNTY20"
)

ma_exposure_matrix <- make_exposure_matrix(
  subset(ma_exposure, COUNTY20 %in% c('MIDDLESEX', 'WORCESTER') &
           year(date) %in% 2012:2015), 
  exposure_columns)
#> Warning in make_exposure_matrix(subset(ma_exposure, COUNTY20 %in% c("MIDDLESEX", : check about any NA, some corrections for this later,
#>             but only in certain columns
```

outcomes

``` r

outcome_columns <- list(
  "date" = "date",
  "outcome" = "daily_deaths",
  "factor" = 'age_grp',
  "factor" = 'sex',
  "geo_unit" = "TOWN20",
  "geo_unit_grp" = "COUNTY20"
)

ma_outcomes_tbl <- make_outcome_table(
  subset(ma_deaths,COUNTY20 %in% c('MIDDLESEX', 'WORCESTER') &
           year(date) %in% 2012:2015), outcome_columns)
```

models

``` r

ma_model <- condPois_2stage(ma_exposure_matrix, ma_outcomes_tbl, verbose = 1, global_cen = 20)
#> -- validation passed
#> -- stage 1
#> -- mixmeta
#> -- stage 2
```

### Estimating the AN

Ok so now you pass in `population`.

So now estimate the AN as a full object

Remember that this needs to be compatible for:

- single zone

- ma model with ma_model\$`_`

- ma model with factor ma_model\$`0-17` \> I think you can handle this
  the same way you did before, with recursion

Now in this second step, you can choose the aggregation level that you
want results to.

In this block you need:

- what spatial resolution are you summarizing to: -\>\> ‘geo_unit’,
  ‘geo_unit_grp’, or ‘all’

- are you just getting the impacts that are \> then the centering point:
  -\>\> lets just assume yes for now, can always go back and change it

``` r

ma_AN <- calc_AN(ma_model, ma_outcomes_tbl, ma_pop_data_long,
                 agg_type = 'TOWN20', 
                 join_cols = 'TOWN20', 
                 nsim = 100,
                 verbose = 2)
#> -- validation passed
#> -- estimate in each geo_unit
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     105     110     
#> -- summarize by simulation
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     
ma_AN$`_`$rate_table
#>          TOWN20  COUNTY20 population above_MMT mean_annual_attr_rate_est
#>          <char>    <char>      <num>    <lgcl>                     <num>
#>   1:      ACTON MIDDLESEX      23864      TRUE                 35.147083
#>   2:      ACTON MIDDLESEX      23864     FALSE                 -2.137110
#>   3:  ARLINGTON MIDDLESEX      45906      TRUE                 35.632597
#>   4:  ARLINGTON MIDDLESEX      45906     FALSE                 -2.458829
#>   5: ASHBURNHAM WORCESTER       6337      TRUE                 24.262269
#>  ---                                                                    
#> 224: WINCHESTER MIDDLESEX      22809     FALSE                 -3.112806
#> 225:     WOBURN MIDDLESEX      40992      TRUE                 34.088969
#> 226:     WOBURN MIDDLESEX      40992     FALSE                 -1.613120
#> 227:  WORCESTER WORCESTER     204191      TRUE                 31.074460
#> 228:  WORCESTER WORCESTER     204191     FALSE                 -2.681925
#>      mean_annual_attr_rate_lb mean_annual_attr_rate_ub
#>                         <num>                    <num>
#>   1:                27.747339                42.456996
#>   2:                -3.337663                -1.083483
#>   3:                28.014720                41.725483
#>   4:                -3.358902                -1.590887
#>   5:                17.394864                30.135514
#>  ---                                                  
#> 224:                -4.236540                -2.113749
#> 225:                28.538648                39.654689
#> 226:                -2.336127                -0.798479
#> 227:                26.747567                35.914842
#> 228:                -3.802537                -1.792960
ma_AN$`_`$number_table
#>          TOWN20  COUNTY20 population above_MMT mean_annual_attr_num_est
#>          <char>    <char>      <num>    <lgcl>                    <num>
#>   1:      ACTON MIDDLESEX      23864      TRUE                  8.38750
#>   2:      ACTON MIDDLESEX      23864     FALSE                 -0.51000
#>   3:  ARLINGTON MIDDLESEX      45906      TRUE                 16.35750
#>   4:  ARLINGTON MIDDLESEX      45906     FALSE                 -1.12875
#>   5: ASHBURNHAM WORCESTER       6337      TRUE                  1.53750
#>  ---                                                                   
#> 224: WINCHESTER MIDDLESEX      22809     FALSE                 -0.71000
#> 225:     WOBURN MIDDLESEX      40992      TRUE                 13.97375
#> 226:     WOBURN MIDDLESEX      40992     FALSE                 -0.66125
#> 227:  WORCESTER WORCESTER     204191      TRUE                 63.45125
#> 228:  WORCESTER WORCESTER     204191     FALSE                 -5.47625
#>      mean_annual_attr_num_lb mean_annual_attr_num_ub
#>                        <num>                   <num>
#>   1:               6.6216250              10.1319375
#>   2:              -0.7965000              -0.2585625
#>   3:              12.8604375              19.1545000
#>   4:              -1.5419375              -0.7303125
#>   5:               1.1023125               1.9096875
#>  ---                                                
#> 224:              -0.9663125              -0.4821250
#> 225:              11.6985625              16.2552500
#> 226:              -0.9576250              -0.3273125
#> 227:              54.6161250              73.3348750
#> 228:              -7.7644375              -3.6610625
```

you can change `agg_type` to be a different spatial resolution – either
whatever the group variable was or “all”

``` r

ma_AN <- calc_AN(ma_model, ma_outcomes_tbl, ma_pop_data_long,
                 agg_type = 'COUNTY20', 
                 join_cols = 'TOWN20', 
                 nsim = 100,
                 verbose = 2)
#> -- validation passed
#> -- estimate in each geo_unit
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     105     110     
#> -- summarize by simulation
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     
ma_AN$`_`$rate_table
#>     COUNTY20 population above_MMT mean_annual_attr_rate_est
#>       <char>      <num>    <lgcl>                     <num>
#> 1: MIDDLESEX    1623109      TRUE                 34.534957
#> 2: MIDDLESEX    1623109     FALSE                 -2.103140
#> 3: WORCESTER     858898      TRUE                 30.078804
#> 4: WORCESTER     858898     FALSE                 -2.199912
#>    mean_annual_attr_rate_lb mean_annual_attr_rate_ub
#>                       <num>                    <num>
#> 1:                33.473387                35.811705
#> 2:                -2.231200                -1.953061
#> 3:                28.579740                31.536261
#> 4:                -2.470106                -1.964167
ma_AN$`_`$number_table
#>     COUNTY20 population above_MMT mean_annual_attr_num_est
#>       <char>      <num>    <lgcl>                    <num>
#> 1: MIDDLESEX    1623109      TRUE                560.54000
#> 2: MIDDLESEX    1623109     FALSE                -34.13625
#> 3: WORCESTER     858898      TRUE                258.34625
#> 4: WORCESTER     858898     FALSE                -18.89500
#>    mean_annual_attr_num_lb mean_annual_attr_num_ub
#>                      <num>                   <num>
#> 1:               543.30956               581.26300
#> 2:               -36.21481               -31.70031
#> 3:               245.47081               270.86431
#> 4:               -21.21569               -16.87019
```

See that the numbers are roughly the same for Suffolk county ? They
won’t be exactly the same because of how the averaging works.

Some plot functions exist:

``` r

plot(ma_AN, table_type = 'rate', above_MMT = T)
```

![](attributable_number_files/figure-html/grpsum3-1.png)

``` r

plot(ma_AN, table_type = 'rate', above_MMT = F)
```

![](attributable_number_files/figure-html/grpsum3-2.png)

### Estimating the AN - single

check of single

``` r


# run the model
m2 <- condPois_1stage(exposure_matrix = ma_exposure_matrix, 
                  outcomes_tbl = ma_outcomes_tbl, 
                  multi_zone = TRUE)
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ACTON. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ARLINGTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ASHBURNHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ASHBY. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ASHLAND. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ATHOL. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit AUBURN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit AYER. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BARRE. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BEDFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BELMONT. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BERLIN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BILLERICA. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BLACKSTONE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BOLTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BOXBOROUGH. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BOYLSTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BROOKFIELD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit BURLINGTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CAMBRIDGE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CARLISLE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CHARLTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CHELMSFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CLINTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit CONCORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit DOUGLAS. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit DRACUT. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit DUDLEY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit DUNSTABLE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit EAST BROOKFIELD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit EVERETT. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit FITCHBURG. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit FRAMINGHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit GARDNER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit GRAFTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit GROTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HARDWICK. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HARVARD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HOLDEN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HOLLISTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HOPEDALE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HOPKINTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HUBBARDSTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit HUDSON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LANCASTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LEICESTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LEOMINSTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LEXINGTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LINCOLN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LITTLETON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LOWELL. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit LUNENBURG. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MALDEN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MARLBOROUGH. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MAYNARD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MEDFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MELROSE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MENDON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MILFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MILLBURY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit MILLVILLE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NATICK. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NEW BRAINTREE. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NEWTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NORTH BROOKFIELD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NORTH READING. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NORTHBOROUGH. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit NORTHBRIDGE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit OAKHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit OXFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit PAXTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit PEPPERELL. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit PETERSHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit PHILLIPSTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit PRINCETON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit READING. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit ROYALSTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit RUTLAND. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SHERBORN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SHIRLEY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SHREWSBURY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SOMERVILLE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SOUTHBOROUGH. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SOUTHBRIDGE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SPENCER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit STERLING. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit STONEHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit STOW. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit STURBRIDGE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SUDBURY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit SUTTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit TEMPLETON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit TEWKSBURY. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit TOWNSEND. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit TYNGSBOROUGH. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit UPTON. This means your zones are across too large of an area, or there
#> are differences in exposures so much that the bases are quite different. Try
#> limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit UXBRIDGE. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WAKEFIELD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WALTHAM. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WARREN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WATERTOWN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WAYLAND. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WEBSTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WEST BOYLSTON. This means your zones are across too large of an area,
#> or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WEST BROOKFIELD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WESTBOROUGH. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WESTFORD. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WESTMINSTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WESTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WILMINGTON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WINCHENDON. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WINCHESTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WOBURN. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.
#> Warning in condPois_1stage(exposure_matrix = ma_exposure_matrix, outcomes_tbl =
#> ma_outcomes_tbl, : Centering point is outside the range of exposures in
#> geo-unit WORCESTER. This means your zones are across too large of an area, or
#> there are differences in exposures so much that the bases are quite different.
#> Try limiting the geo-units passed in to those that are more similar, manually
#> setting a centering point that you know each geo-unit has, or changing your
#> exposure variable.

ma_AN_s1 <- calc_AN(m2, ma_outcomes_tbl, ma_pop_data_long,
                 agg_type = 'COUNTY20', 
                 join_cols = 'TOWN20', 
                 nsim = 100,
                 verbose = 2)
#> -- validation passed
#> -- estimate in each geo_unit
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     105     110     
#> -- summarize by simulation
#> 5    10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100     

ma_AN_s1$`_`$rate_table
#>     COUNTY20 population above_MMT mean_annual_attr_rate_est
#>       <char>      <num>    <lgcl>                     <num>
#> 1: MIDDLESEX    1623109      TRUE                  83.09755
#> 2: MIDDLESEX    1623109     FALSE                   0.00000
#> 3: WORCESTER     858898      TRUE                  80.60896
#> 4: WORCESTER     858898     FALSE                   0.00000
#>    mean_annual_attr_rate_lb mean_annual_attr_rate_ub
#>                       <num>                    <num>
#> 1:                 82.14084                 84.07091
#> 2:                  0.00000                  0.00000
#> 3:                 78.76969                 81.88347
#> 4:                  0.00000                  0.00000
plot(ma_AN_s1, "num", above_MMT = T)
```

![](attributable_number_files/figure-html/checkSingle1-1.png)

### Estimating the AN - with factors

In the case where you have factors, you can easily extend this

``` r

ma_outcomes_tbl_fct <- make_outcome_table(
  subset(ma_deaths,COUNTY20 %in% c('MIDDLESEX', 'WORCESTER') &
           year(date) %in% 2012:2015), 
  outcome_columns, collapse_to = 'age_grp')

ma_model_fct <- condPois_2stage(ma_exposure_matrix, ma_outcomes_tbl_fct, verbose = 1)
#> < age_grp : 0-17 >
#> -- validation passed
#> -- stage 1
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BILLERICA. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BLACKSTONE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BOLTON. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BOXBOROUGH. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit CARLISLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit CONCORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit DRACUT. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FITCHBURG. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FRAMINGHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit GARDNER. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit GROTON. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HARDWICK. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOLDEN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOPKINTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LANCASTER. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LEXINGTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LITTLETON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LOWELL. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MARLBOROUGH. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MAYNARD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MILFORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit NORTHBRIDGE. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit READING. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SHREWSBURY. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOMERVILLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOUTHBOROUGH. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOUTHBRIDGE. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STERLING. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STONEHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STURBRIDGE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit TOWNSEND. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit UXBRIDGE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WAKEFIELD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WALTHAM. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WEBSTER. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WESTFORD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WESTON. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WINCHENDON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WINCHESTER. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WOBURN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> -- mixmeta
#> -- stage 2
#> < age_grp : 18-64 >
#> -- validation passed
#> -- stage 1
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BARRE. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BILLERICA. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit CARLISLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit CONCORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit DRACUT. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FITCHBURG. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FRAMINGHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit GARDNER. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOPEDALE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOPKINTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LEXINGTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LOWELL. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MARLBOROUGH. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MILFORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MILLVILLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit OAKHAM. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit READING. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SHREWSBURY. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOMERVILLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOUTHBRIDGE. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STERLING. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STONEHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STURBRIDGE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WAKEFIELD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WALTHAM. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WESTFORD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WINCHESTER. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WOBURN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> -- mixmeta
#> -- stage 2
#> < age_grp : 65+ >
#> -- validation passed
#> -- stage 1
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BERLIN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BILLERICA. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit BLACKSTONE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit CONCORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit DRACUT. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit DUNSTABLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FITCHBURG. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit FRAMINGHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit GARDNER. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit GROTON. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOLDEN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit HOPKINTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LEXINGTON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LITTLETON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit LOWELL. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MARLBOROUGH. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MAYNARD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit MILFORD. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit NORTHBRIDGE. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit READING. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SHREWSBURY. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOMERVILLE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOUTHBOROUGH. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit SOUTHBRIDGE. This means your zones are across too large
#> of an area, or there are differences in exposures so much that the bases are
#> quite different. Try limiting the geo-units passed in to those that are more
#> similar, manually setting a centering point that you know each geo-unit has, or
#> changing your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STONEHAM. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit STURBRIDGE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit UXBRIDGE. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WAKEFIELD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WALTHAM. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WEBSTER. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WESTFORD. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WESTON. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WINCHENDON. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WINCHESTER. This means your zones are across too large of
#> an area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> Warning in condPois_1stage(exposure_matrix = single_exposure_matrix,
#> outcomes_tbl = single_outcomes_tbl, : Centering point is outside the range of
#> exposures in geo-unit WOBURN. This means your zones are across too large of an
#> area, or there are differences in exposures so much that the bases are quite
#> different. Try limiting the geo-units passed in to those that are more similar,
#> manually setting a centering point that you know each geo-unit has, or changing
#> your exposure variable.
#> -- mixmeta
#> -- stage 2

ma_AN_fct <- calc_AN(ma_model_fct, ma_outcomes_tbl_fct,
                     ma_pop_data_long,
                 agg_type = 'COUNTY20', 
                 join_cols = 'TOWN20', 
                 nsim = 100,
                 verbose = 1)
#> < age_grp : 0-17 >
#> -- validation passed
#> -- estimate in each geo_unit
#> -- summarize by simulation
#> < age_grp : 18-64 >
#> -- validation passed
#> -- estimate in each geo_unit
#> -- summarize by simulation
#> < age_grp : 65+ >
#> -- validation passed
#> -- estimate in each geo_unit
#> -- summarize by simulation

plot(ma_AN_fct, "num", above_MMT = T)
```

![](attributable_number_files/figure-html/fctrun-1.png)

These results are fictional of course but show what kind of outputs can
be made easily.

``` r

spatial_plot(ma_AN_fct, shp = ma_counties, table_type = "num", above_MMT = T)
```

![](attributable_number_files/figure-html/multi_plot3d-1.png)
