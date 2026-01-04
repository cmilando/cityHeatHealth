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
#>          <char>     <Date>    <num>    <num>    <num>    <num>    <num>
#>    1: MIDDLESEX 2010-05-01 22.73500 15.94895  8.32810 10.25350 15.61525
#>    2:   SUFFOLK 2010-05-01 22.97060 14.95275  7.70270 10.44610 15.55035
#>    3: MIDDLESEX 2010-05-02 29.20320 22.73500 15.94895  8.32810 10.25350
#>    4:   SUFFOLK 2010-05-02 27.38805 22.97060 14.95275  7.70270 10.44610
#>    5: MIDDLESEX 2010-05-03 30.83745 29.20320 22.73500 15.94895  8.32810
#>   ---                                                                  
#> 3362:   SUFFOLK 2020-09-28 25.19490 24.72725 25.15465 24.78400 24.66000
#> 3363: MIDDLESEX 2020-09-29 24.69490 24.29570 23.52120 26.10840 23.63025
#> 3364:   SUFFOLK 2020-09-29 26.22930 25.19490 24.72725 25.15465 24.78400
#> 3365: MIDDLESEX 2020-09-30 22.34675 24.69490 24.29570 23.52120 26.10840
#> 3366:   SUFFOLK 2020-09-30 25.48680 26.22930 25.19490 24.72725 25.15465
#>        explag5
#>          <num>
#>    1: 13.13640
#>    2: 16.71785
#>    3: 15.61525
#>    4: 15.55035
#>    5: 10.25350
#>   ---         
#> 3362: 18.28090
#> 3363: 23.57620
#> 3364: 24.66000
#> 3365: 23.63025
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
#>           <Date>    <char>        <int>      <char>                    <char>
#>    1: 2010-05-01 MIDDLESEX          395         ALL MIDDLESEX:yr2010:mn5:dow7
#>    2: 2010-05-01   SUFFOLK          329         ALL   SUFFOLK:yr2010:mn5:dow7
#>    3: 2010-05-02 MIDDLESEX          406         ALL MIDDLESEX:yr2010:mn5:dow1
#>    4: 2010-05-02   SUFFOLK          336         ALL   SUFFOLK:yr2010:mn5:dow1
#>    5: 2010-05-03 MIDDLESEX          411         ALL MIDDLESEX:yr2010:mn5:dow2
#>   ---                                                                        
#> 3362: 2020-09-28   SUFFOLK          467         ALL   SUFFOLK:yr2020:mn9:dow2
#> 3363: 2020-09-29 MIDDLESEX          336         ALL MIDDLESEX:yr2020:mn9:dow3
#> 3364: 2020-09-29   SUFFOLK          459         ALL   SUFFOLK:yr2020:mn9:dow3
#> 3365: 2020-09-30 MIDDLESEX          320         ALL MIDDLESEX:yr2020:mn9:dow4
#> 3366: 2020-09-30   SUFFOLK          437         ALL   SUFFOLK:yr2020:mn9:dow4
#>       strata_total
#>              <num>
#>    1:         1987
#>    2:         1684
#>    3:         2018
#>    4:         1675
#>    5:         1951
#>   ---             
#> 3362:         1704
#> 3363:         1538
#> 3364:         2144
#> 3365:         1561
#> 3366:         2094

#
stopifnot(identical(deaths_tbl$date, exposure_mat$date ))
```
