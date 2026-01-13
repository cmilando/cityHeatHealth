# group-check

``` r

library(cityHeatHealth)
```

``` r

data("ma_exposure")
data("ma_deaths")

# create exposure matrix
exposure_columns <- list(
  "date" = "date",
  "exposure" = "tmax_C",
  "geo_unit" = "TOWN20",
  "geo_unit_grp" = "COUNTY20"
)

TOWNLIST <- c('CHELSEA', 'EVERETT', 'REVERE', 'MALDEN')

exposure <- subset(ma_exposure, TOWN20 %in%  TOWNLIST)
summary(exposure)
#>       date                tmax_C           TOWN20            COUNTY20        
#>  Min.   :2010-01-01   Min.   :-17.128   Length:15731       Length:15731      
#>  1st Qu.:2012-09-28   1st Qu.:  5.814   Class :character   Class :character  
#>  Median :2015-06-29   Median : 15.670   Mode  :character   Mode  :character  
#>  Mean   :2015-07-01   Mean   : 14.968                                        
#>  3rd Qu.:2018-04-03   3rd Qu.: 24.552                                        
#>  Max.   :2020-12-31   Max.   : 38.557                                        
#>                       NA's   :148

exposure_mat <- make_exposure_matrix(exposure, exposure_columns, 
                                     grp_level = T)
#> Warning in make_exposure_matrix(exposure, exposure_columns, grp_level = T): check about any NA, some corrections for this later,
#>             but only in certain columns

exposure_mat
#>        COUNTY20       date   tmax_C  explag1  explag2  explag3  explag4
#>          <char>     <IDat>    <num>    <num>    <num>    <num>    <num>
#>    1: MIDDLESEX 2010-05-01 22.73500 15.94895  8.32810 10.25350 15.61525
#>    2: MIDDLESEX 2010-05-02 29.20320 22.73500 15.94895  8.32810 10.25350
#>    3: MIDDLESEX 2010-05-03 30.83745 29.20320 22.73500 15.94895  8.32810
#>    4: MIDDLESEX 2010-05-04 25.41450 30.83745 29.20320 22.73500 15.94895
#>    5: MIDDLESEX 2010-05-05 24.12005 25.41450 30.83745 29.20320 22.73500
#>   ---                                                                  
#> 3362:   SUFFOLK 2020-09-26 25.15465 24.78400 24.66000 18.28090 17.02280
#> 3363:   SUFFOLK 2020-09-27 24.72725 25.15465 24.78400 24.66000 18.28090
#> 3364:   SUFFOLK 2020-09-28 25.19490 24.72725 25.15465 24.78400 24.66000
#> 3365:   SUFFOLK 2020-09-29 26.22930 25.19490 24.72725 25.15465 24.78400
#> 3366:   SUFFOLK 2020-09-30 25.48680 26.22930 25.19490 24.72725 25.15465
#>        explag5
#>          <num>
#>    1: 13.13640
#>    2: 15.61525
#>    3: 10.25350
#>    4:  8.32810
#>    5: 15.94895
#>   ---         
#> 3362: 15.72065
#> 3363: 17.02280
#> 3364: 18.28090
#> 3365: 24.66000
#> 3366: 24.78400

# create outcome table
outcome_columns <- list(
  "date" = "date",
  "outcome" = "daily_deaths",
  "factor" = 'age_grp',
  "factor" = 'sex',
  "geo_unit" = "TOWN20",
  "geo_unit_grp" = "COUNTY20"
)
deaths   <- subset(ma_deaths, TOWN20 %in%  TOWNLIST)
deaths_tbl <- make_outcome_table(deaths,  outcome_columns, 
                                 grp_level = T)
#> Warning in make_outcome_table(deaths, outcome_columns, grp_level = T): make type checks  (e.g., so Date == Date),
#>          for some reason this doesn't work in some cases? but ok in others?
deaths_tbl
#>             date  COUNTY20 daily_deaths spatial_grp                    strata
#>           <IDat>    <char>        <int>      <char>                    <char>
#>    1: 2010-05-01 MIDDLESEX          395         ALL MIDDLESEX:yr2010:mn5:dow7
#>    2: 2010-05-02 MIDDLESEX          406         ALL MIDDLESEX:yr2010:mn5:dow1
#>    3: 2010-05-03 MIDDLESEX          411         ALL MIDDLESEX:yr2010:mn5:dow2
#>    4: 2010-05-04 MIDDLESEX          411         ALL MIDDLESEX:yr2010:mn5:dow3
#>    5: 2010-05-05 MIDDLESEX          427         ALL MIDDLESEX:yr2010:mn5:dow4
#>   ---                                                                        
#> 3362: 2020-09-26   SUFFOLK          421         ALL   SUFFOLK:yr2020:mn9:dow7
#> 3363: 2020-09-27   SUFFOLK          433         ALL   SUFFOLK:yr2020:mn9:dow1
#> 3364: 2020-09-28   SUFFOLK          467         ALL   SUFFOLK:yr2020:mn9:dow2
#> 3365: 2020-09-29   SUFFOLK          459         ALL   SUFFOLK:yr2020:mn9:dow3
#> 3366: 2020-09-30   SUFFOLK          437         ALL   SUFFOLK:yr2020:mn9:dow4
#>       strata_total
#>              <num>
#>    1:         1987
#>    2:         2018
#>    3:         1951
#>    4:         1576
#>    5:         1615
#>   ---             
#> 3362:         1728
#> 3363:         1674
#> 3364:         1704
#> 3365:         2144
#> 3366:         2094

#
stopifnot(identical(deaths_tbl$date, exposure_mat$date ))
```
