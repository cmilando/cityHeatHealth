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
head(exposure)
#>              date  tmax_C  TOWN20 COUNTY20
#> 224432 2010-01-01 -0.7221 CHELSEA  SUFFOLK
#> 224433 2010-01-02  1.2943 CHELSEA  SUFFOLK
#> 224434 2010-01-03 -1.6351 CHELSEA  SUFFOLK
#> 224435 2010-01-04 -0.6948 CHELSEA  SUFFOLK
#> 224436 2010-01-05  0.6686 CHELSEA  SUFFOLK
#> 224437 2010-01-06  0.9271 CHELSEA  SUFFOLK

# so this has strata inncorrect because
# it should only be by town if keep_unit = True
exposure_mat <- make_exposure_matrix(exposure, exposure_columns, 
                                     grp_level = T)
#> Warning in make_exposure_matrix(exposure, exposure_columns, grp_level = T): check about any NA, some corrections for this later,
#>             but only in certain columns

exposure_mat
#>        COUNTY20 spatial_grp       date                    strata   tmax_C
#>          <char>      <char>     <IDat>                    <char>    <num>
#>    1: MIDDLESEX         ALL 2010-05-02 MIDDLESEX:yr2010:mn5:dow1 29.20320
#>    2: MIDDLESEX         ALL 2010-05-09 MIDDLESEX:yr2010:mn5:dow1 16.72915
#>    3: MIDDLESEX         ALL 2010-05-16 MIDDLESEX:yr2010:mn5:dow1 19.63385
#>    4: MIDDLESEX         ALL 2010-05-23 MIDDLESEX:yr2010:mn5:dow1 24.13790
#>    5: MIDDLESEX         ALL 2010-05-30 MIDDLESEX:yr2010:mn5:dow1 24.61435
#>   ---                                                                    
#> 3362:   SUFFOLK         ALL 2020-09-25   SUFFOLK:yr2020:mn9:dow6 24.78400
#> 3363:   SUFFOLK         ALL 2020-09-05   SUFFOLK:yr2020:mn9:dow7 27.85625
#> 3364:   SUFFOLK         ALL 2020-09-12   SUFFOLK:yr2020:mn9:dow7 20.56520
#> 3365:   SUFFOLK         ALL 2020-09-19   SUFFOLK:yr2020:mn9:dow7 17.62225
#> 3366:   SUFFOLK         ALL 2020-09-26   SUFFOLK:yr2020:mn9:dow7 25.15465
#>        explag1  explag2  explag3  explag4  explag5
#>          <num>    <num>    <num>    <num>    <num>
#>    1: 22.73500 15.94895  8.32810 10.25350 15.61525
#>    2: 20.17370 23.08275 25.47170 24.12005 25.41450
#>    3: 23.74405 18.69100  7.95205 13.63525 12.07127
#>    4: 28.12350 27.39070 13.89735 17.25370 23.89320
#>    5: 25.24900 24.75090 34.15800 31.06935 26.86420
#>   ---                                             
#> 3362: 24.66000 18.28090 17.02280 15.72065 16.50645
#> 3363: 26.16430 21.98300 24.33585 23.54830 22.96530
#> 3364: 27.49125 27.96035 29.19290 27.01440 25.66710
#> 3365: 23.84470 22.32660 17.85965 22.99420 23.41645
#> 3366: 24.78400 24.66000 18.28090 17.02280 15.72065

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

# seems like this isn't working because its
deaths_tbl <- make_outcome_table(deaths,  outcome_columns, 
                                 grp_level = T)
deaths_tbl
#>             date  COUNTY20 daily_deaths spatial_grp                    strata
#>           <IDat>    <char>        <int>      <char>                    <char>
#>    1: 2010-05-02 MIDDLESEX          406         ALL MIDDLESEX:yr2010:mn5:dow1
#>    2: 2010-05-09 MIDDLESEX          399         ALL MIDDLESEX:yr2010:mn5:dow1
#>    3: 2010-05-16 MIDDLESEX          368         ALL MIDDLESEX:yr2010:mn5:dow1
#>    4: 2010-05-23 MIDDLESEX          398         ALL MIDDLESEX:yr2010:mn5:dow1
#>    5: 2010-05-30 MIDDLESEX          447         ALL MIDDLESEX:yr2010:mn5:dow1
#>   ---                                                                        
#> 3362: 2020-09-25   SUFFOLK          421         ALL   SUFFOLK:yr2020:mn9:dow6
#> 3363: 2020-09-05   SUFFOLK          458         ALL   SUFFOLK:yr2020:mn9:dow7
#> 3364: 2020-09-12   SUFFOLK          435         ALL   SUFFOLK:yr2020:mn9:dow7
#> 3365: 2020-09-19   SUFFOLK          414         ALL   SUFFOLK:yr2020:mn9:dow7
#> 3366: 2020-09-26   SUFFOLK          421         ALL   SUFFOLK:yr2020:mn9:dow7
#>       strata_total
#>              <num>
#>    1:         2018
#>    2:         2018
#>    3:         2018
#>    4:         2018
#>    5:         2018
#>   ---             
#> 3362:         1746
#> 3363:         1728
#> 3364:         1728
#> 3365:         1728
#> 3366:         1728
unique(deaths_tbl$strata)
#>   [1] "MIDDLESEX:yr2010:mn5:dow1" "MIDDLESEX:yr2010:mn5:dow2"
#>   [3] "MIDDLESEX:yr2010:mn5:dow3" "MIDDLESEX:yr2010:mn5:dow4"
#>   [5] "MIDDLESEX:yr2010:mn5:dow5" "MIDDLESEX:yr2010:mn5:dow6"
#>   [7] "MIDDLESEX:yr2010:mn5:dow7" "MIDDLESEX:yr2010:mn6:dow1"
#>   [9] "MIDDLESEX:yr2010:mn6:dow2" "MIDDLESEX:yr2010:mn6:dow3"
#>  [11] "MIDDLESEX:yr2010:mn6:dow4" "MIDDLESEX:yr2010:mn6:dow5"
#>  [13] "MIDDLESEX:yr2010:mn6:dow6" "MIDDLESEX:yr2010:mn6:dow7"
#>  [15] "MIDDLESEX:yr2010:mn7:dow1" "MIDDLESEX:yr2010:mn7:dow2"
#>  [17] "MIDDLESEX:yr2010:mn7:dow3" "MIDDLESEX:yr2010:mn7:dow4"
#>  [19] "MIDDLESEX:yr2010:mn7:dow5" "MIDDLESEX:yr2010:mn7:dow6"
#>  [21] "MIDDLESEX:yr2010:mn7:dow7" "MIDDLESEX:yr2010:mn8:dow1"
#>  [23] "MIDDLESEX:yr2010:mn8:dow2" "MIDDLESEX:yr2010:mn8:dow3"
#>  [25] "MIDDLESEX:yr2010:mn8:dow4" "MIDDLESEX:yr2010:mn8:dow5"
#>  [27] "MIDDLESEX:yr2010:mn8:dow6" "MIDDLESEX:yr2010:mn8:dow7"
#>  [29] "MIDDLESEX:yr2010:mn9:dow1" "MIDDLESEX:yr2010:mn9:dow2"
#>  [31] "MIDDLESEX:yr2010:mn9:dow3" "MIDDLESEX:yr2010:mn9:dow4"
#>  [33] "MIDDLESEX:yr2010:mn9:dow5" "MIDDLESEX:yr2010:mn9:dow6"
#>  [35] "MIDDLESEX:yr2010:mn9:dow7" "MIDDLESEX:yr2011:mn5:dow1"
#>  [37] "MIDDLESEX:yr2011:mn5:dow2" "MIDDLESEX:yr2011:mn5:dow3"
#>  [39] "MIDDLESEX:yr2011:mn5:dow4" "MIDDLESEX:yr2011:mn5:dow5"
#>  [41] "MIDDLESEX:yr2011:mn5:dow6" "MIDDLESEX:yr2011:mn5:dow7"
#>  [43] "MIDDLESEX:yr2011:mn6:dow1" "MIDDLESEX:yr2011:mn6:dow2"
#>  [45] "MIDDLESEX:yr2011:mn6:dow3" "MIDDLESEX:yr2011:mn6:dow4"
#>  [47] "MIDDLESEX:yr2011:mn6:dow5" "MIDDLESEX:yr2011:mn6:dow6"
#>  [49] "MIDDLESEX:yr2011:mn6:dow7" "MIDDLESEX:yr2011:mn7:dow1"
#>  [51] "MIDDLESEX:yr2011:mn7:dow2" "MIDDLESEX:yr2011:mn7:dow3"
#>  [53] "MIDDLESEX:yr2011:mn7:dow4" "MIDDLESEX:yr2011:mn7:dow5"
#>  [55] "MIDDLESEX:yr2011:mn7:dow6" "MIDDLESEX:yr2011:mn7:dow7"
#>  [57] "MIDDLESEX:yr2011:mn8:dow1" "MIDDLESEX:yr2011:mn8:dow2"
#>  [59] "MIDDLESEX:yr2011:mn8:dow3" "MIDDLESEX:yr2011:mn8:dow4"
#>  [61] "MIDDLESEX:yr2011:mn8:dow5" "MIDDLESEX:yr2011:mn8:dow6"
#>  [63] "MIDDLESEX:yr2011:mn8:dow7" "MIDDLESEX:yr2011:mn9:dow1"
#>  [65] "MIDDLESEX:yr2011:mn9:dow2" "MIDDLESEX:yr2011:mn9:dow3"
#>  [67] "MIDDLESEX:yr2011:mn9:dow4" "MIDDLESEX:yr2011:mn9:dow5"
#>  [69] "MIDDLESEX:yr2011:mn9:dow6" "MIDDLESEX:yr2011:mn9:dow7"
#>  [71] "MIDDLESEX:yr2012:mn5:dow1" "MIDDLESEX:yr2012:mn5:dow2"
#>  [73] "MIDDLESEX:yr2012:mn5:dow3" "MIDDLESEX:yr2012:mn5:dow4"
#>  [75] "MIDDLESEX:yr2012:mn5:dow5" "MIDDLESEX:yr2012:mn5:dow6"
#>  [77] "MIDDLESEX:yr2012:mn5:dow7" "MIDDLESEX:yr2012:mn6:dow1"
#>  [79] "MIDDLESEX:yr2012:mn6:dow2" "MIDDLESEX:yr2012:mn6:dow3"
#>  [81] "MIDDLESEX:yr2012:mn6:dow4" "MIDDLESEX:yr2012:mn6:dow5"
#>  [83] "MIDDLESEX:yr2012:mn6:dow6" "MIDDLESEX:yr2012:mn6:dow7"
#>  [85] "MIDDLESEX:yr2012:mn7:dow1" "MIDDLESEX:yr2012:mn7:dow2"
#>  [87] "MIDDLESEX:yr2012:mn7:dow3" "MIDDLESEX:yr2012:mn7:dow4"
#>  [89] "MIDDLESEX:yr2012:mn7:dow5" "MIDDLESEX:yr2012:mn7:dow6"
#>  [91] "MIDDLESEX:yr2012:mn7:dow7" "MIDDLESEX:yr2012:mn8:dow1"
#>  [93] "MIDDLESEX:yr2012:mn8:dow2" "MIDDLESEX:yr2012:mn8:dow3"
#>  [95] "MIDDLESEX:yr2012:mn8:dow4" "MIDDLESEX:yr2012:mn8:dow5"
#>  [97] "MIDDLESEX:yr2012:mn8:dow6" "MIDDLESEX:yr2012:mn8:dow7"
#>  [99] "MIDDLESEX:yr2012:mn9:dow1" "MIDDLESEX:yr2012:mn9:dow2"
#> [101] "MIDDLESEX:yr2012:mn9:dow3" "MIDDLESEX:yr2012:mn9:dow4"
#> [103] "MIDDLESEX:yr2012:mn9:dow5" "MIDDLESEX:yr2012:mn9:dow6"
#> [105] "MIDDLESEX:yr2012:mn9:dow7" "MIDDLESEX:yr2013:mn5:dow1"
#> [107] "MIDDLESEX:yr2013:mn5:dow2" "MIDDLESEX:yr2013:mn5:dow3"
#> [109] "MIDDLESEX:yr2013:mn5:dow4" "MIDDLESEX:yr2013:mn5:dow5"
#> [111] "MIDDLESEX:yr2013:mn5:dow6" "MIDDLESEX:yr2013:mn5:dow7"
#> [113] "MIDDLESEX:yr2013:mn6:dow1" "MIDDLESEX:yr2013:mn6:dow2"
#> [115] "MIDDLESEX:yr2013:mn6:dow3" "MIDDLESEX:yr2013:mn6:dow4"
#> [117] "MIDDLESEX:yr2013:mn6:dow5" "MIDDLESEX:yr2013:mn6:dow6"
#> [119] "MIDDLESEX:yr2013:mn6:dow7" "MIDDLESEX:yr2013:mn7:dow1"
#> [121] "MIDDLESEX:yr2013:mn7:dow2" "MIDDLESEX:yr2013:mn7:dow3"
#> [123] "MIDDLESEX:yr2013:mn7:dow4" "MIDDLESEX:yr2013:mn7:dow5"
#> [125] "MIDDLESEX:yr2013:mn7:dow6" "MIDDLESEX:yr2013:mn7:dow7"
#> [127] "MIDDLESEX:yr2013:mn8:dow1" "MIDDLESEX:yr2013:mn8:dow2"
#> [129] "MIDDLESEX:yr2013:mn8:dow3" "MIDDLESEX:yr2013:mn8:dow4"
#> [131] "MIDDLESEX:yr2013:mn8:dow5" "MIDDLESEX:yr2013:mn8:dow6"
#> [133] "MIDDLESEX:yr2013:mn8:dow7" "MIDDLESEX:yr2013:mn9:dow1"
#> [135] "MIDDLESEX:yr2013:mn9:dow2" "MIDDLESEX:yr2013:mn9:dow3"
#> [137] "MIDDLESEX:yr2013:mn9:dow4" "MIDDLESEX:yr2013:mn9:dow5"
#> [139] "MIDDLESEX:yr2013:mn9:dow6" "MIDDLESEX:yr2013:mn9:dow7"
#> [141] "MIDDLESEX:yr2014:mn5:dow1" "MIDDLESEX:yr2014:mn5:dow2"
#> [143] "MIDDLESEX:yr2014:mn5:dow3" "MIDDLESEX:yr2014:mn5:dow4"
#> [145] "MIDDLESEX:yr2014:mn5:dow5" "MIDDLESEX:yr2014:mn5:dow6"
#> [147] "MIDDLESEX:yr2014:mn5:dow7" "MIDDLESEX:yr2014:mn6:dow1"
#> [149] "MIDDLESEX:yr2014:mn6:dow2" "MIDDLESEX:yr2014:mn6:dow3"
#> [151] "MIDDLESEX:yr2014:mn6:dow4" "MIDDLESEX:yr2014:mn6:dow5"
#> [153] "MIDDLESEX:yr2014:mn6:dow6" "MIDDLESEX:yr2014:mn6:dow7"
#> [155] "MIDDLESEX:yr2014:mn7:dow1" "MIDDLESEX:yr2014:mn7:dow2"
#> [157] "MIDDLESEX:yr2014:mn7:dow3" "MIDDLESEX:yr2014:mn7:dow4"
#> [159] "MIDDLESEX:yr2014:mn7:dow5" "MIDDLESEX:yr2014:mn7:dow6"
#> [161] "MIDDLESEX:yr2014:mn7:dow7" "MIDDLESEX:yr2014:mn8:dow1"
#> [163] "MIDDLESEX:yr2014:mn8:dow2" "MIDDLESEX:yr2014:mn8:dow3"
#> [165] "MIDDLESEX:yr2014:mn8:dow4" "MIDDLESEX:yr2014:mn8:dow5"
#> [167] "MIDDLESEX:yr2014:mn8:dow6" "MIDDLESEX:yr2014:mn8:dow7"
#> [169] "MIDDLESEX:yr2014:mn9:dow1" "MIDDLESEX:yr2014:mn9:dow2"
#> [171] "MIDDLESEX:yr2014:mn9:dow3" "MIDDLESEX:yr2014:mn9:dow4"
#> [173] "MIDDLESEX:yr2014:mn9:dow5" "MIDDLESEX:yr2014:mn9:dow6"
#> [175] "MIDDLESEX:yr2014:mn9:dow7" "MIDDLESEX:yr2015:mn5:dow1"
#> [177] "MIDDLESEX:yr2015:mn5:dow2" "MIDDLESEX:yr2015:mn5:dow3"
#> [179] "MIDDLESEX:yr2015:mn5:dow4" "MIDDLESEX:yr2015:mn5:dow5"
#> [181] "MIDDLESEX:yr2015:mn5:dow6" "MIDDLESEX:yr2015:mn5:dow7"
#> [183] "MIDDLESEX:yr2015:mn6:dow1" "MIDDLESEX:yr2015:mn6:dow2"
#> [185] "MIDDLESEX:yr2015:mn6:dow3" "MIDDLESEX:yr2015:mn6:dow4"
#> [187] "MIDDLESEX:yr2015:mn6:dow5" "MIDDLESEX:yr2015:mn6:dow6"
#> [189] "MIDDLESEX:yr2015:mn6:dow7" "MIDDLESEX:yr2015:mn7:dow1"
#> [191] "MIDDLESEX:yr2015:mn7:dow2" "MIDDLESEX:yr2015:mn7:dow3"
#> [193] "MIDDLESEX:yr2015:mn7:dow4" "MIDDLESEX:yr2015:mn7:dow5"
#> [195] "MIDDLESEX:yr2015:mn7:dow6" "MIDDLESEX:yr2015:mn7:dow7"
#> [197] "MIDDLESEX:yr2015:mn8:dow1" "MIDDLESEX:yr2015:mn8:dow2"
#> [199] "MIDDLESEX:yr2015:mn8:dow3" "MIDDLESEX:yr2015:mn8:dow4"
#> [201] "MIDDLESEX:yr2015:mn8:dow5" "MIDDLESEX:yr2015:mn8:dow6"
#> [203] "MIDDLESEX:yr2015:mn8:dow7" "MIDDLESEX:yr2015:mn9:dow1"
#> [205] "MIDDLESEX:yr2015:mn9:dow2" "MIDDLESEX:yr2015:mn9:dow3"
#> [207] "MIDDLESEX:yr2015:mn9:dow4" "MIDDLESEX:yr2015:mn9:dow5"
#> [209] "MIDDLESEX:yr2015:mn9:dow6" "MIDDLESEX:yr2015:mn9:dow7"
#> [211] "MIDDLESEX:yr2016:mn5:dow1" "MIDDLESEX:yr2016:mn5:dow2"
#> [213] "MIDDLESEX:yr2016:mn5:dow3" "MIDDLESEX:yr2016:mn5:dow4"
#> [215] "MIDDLESEX:yr2016:mn5:dow5" "MIDDLESEX:yr2016:mn5:dow6"
#> [217] "MIDDLESEX:yr2016:mn5:dow7" "MIDDLESEX:yr2016:mn6:dow1"
#> [219] "MIDDLESEX:yr2016:mn6:dow2" "MIDDLESEX:yr2016:mn6:dow3"
#> [221] "MIDDLESEX:yr2016:mn6:dow4" "MIDDLESEX:yr2016:mn6:dow5"
#> [223] "MIDDLESEX:yr2016:mn6:dow6" "MIDDLESEX:yr2016:mn6:dow7"
#> [225] "MIDDLESEX:yr2016:mn7:dow1" "MIDDLESEX:yr2016:mn7:dow2"
#> [227] "MIDDLESEX:yr2016:mn7:dow3" "MIDDLESEX:yr2016:mn7:dow4"
#> [229] "MIDDLESEX:yr2016:mn7:dow5" "MIDDLESEX:yr2016:mn7:dow6"
#> [231] "MIDDLESEX:yr2016:mn7:dow7" "MIDDLESEX:yr2016:mn8:dow1"
#> [233] "MIDDLESEX:yr2016:mn8:dow2" "MIDDLESEX:yr2016:mn8:dow3"
#> [235] "MIDDLESEX:yr2016:mn8:dow4" "MIDDLESEX:yr2016:mn8:dow5"
#> [237] "MIDDLESEX:yr2016:mn8:dow6" "MIDDLESEX:yr2016:mn8:dow7"
#> [239] "MIDDLESEX:yr2016:mn9:dow1" "MIDDLESEX:yr2016:mn9:dow2"
#> [241] "MIDDLESEX:yr2016:mn9:dow3" "MIDDLESEX:yr2016:mn9:dow4"
#> [243] "MIDDLESEX:yr2016:mn9:dow5" "MIDDLESEX:yr2016:mn9:dow6"
#> [245] "MIDDLESEX:yr2016:mn9:dow7" "MIDDLESEX:yr2017:mn5:dow1"
#> [247] "MIDDLESEX:yr2017:mn5:dow2" "MIDDLESEX:yr2017:mn5:dow3"
#> [249] "MIDDLESEX:yr2017:mn5:dow4" "MIDDLESEX:yr2017:mn5:dow5"
#> [251] "MIDDLESEX:yr2017:mn5:dow6" "MIDDLESEX:yr2017:mn5:dow7"
#> [253] "MIDDLESEX:yr2017:mn6:dow1" "MIDDLESEX:yr2017:mn6:dow2"
#> [255] "MIDDLESEX:yr2017:mn6:dow3" "MIDDLESEX:yr2017:mn6:dow4"
#> [257] "MIDDLESEX:yr2017:mn6:dow5" "MIDDLESEX:yr2017:mn6:dow6"
#> [259] "MIDDLESEX:yr2017:mn6:dow7" "MIDDLESEX:yr2017:mn7:dow1"
#> [261] "MIDDLESEX:yr2017:mn7:dow2" "MIDDLESEX:yr2017:mn7:dow3"
#> [263] "MIDDLESEX:yr2017:mn7:dow4" "MIDDLESEX:yr2017:mn7:dow5"
#> [265] "MIDDLESEX:yr2017:mn7:dow6" "MIDDLESEX:yr2017:mn7:dow7"
#> [267] "MIDDLESEX:yr2017:mn8:dow1" "MIDDLESEX:yr2017:mn8:dow2"
#> [269] "MIDDLESEX:yr2017:mn8:dow3" "MIDDLESEX:yr2017:mn8:dow4"
#> [271] "MIDDLESEX:yr2017:mn8:dow5" "MIDDLESEX:yr2017:mn8:dow6"
#> [273] "MIDDLESEX:yr2017:mn8:dow7" "MIDDLESEX:yr2017:mn9:dow1"
#> [275] "MIDDLESEX:yr2017:mn9:dow2" "MIDDLESEX:yr2017:mn9:dow3"
#> [277] "MIDDLESEX:yr2017:mn9:dow4" "MIDDLESEX:yr2017:mn9:dow5"
#> [279] "MIDDLESEX:yr2017:mn9:dow6" "MIDDLESEX:yr2017:mn9:dow7"
#> [281] "MIDDLESEX:yr2018:mn5:dow1" "MIDDLESEX:yr2018:mn5:dow2"
#> [283] "MIDDLESEX:yr2018:mn5:dow3" "MIDDLESEX:yr2018:mn5:dow4"
#> [285] "MIDDLESEX:yr2018:mn5:dow5" "MIDDLESEX:yr2018:mn5:dow6"
#> [287] "MIDDLESEX:yr2018:mn5:dow7" "MIDDLESEX:yr2018:mn6:dow1"
#> [289] "MIDDLESEX:yr2018:mn6:dow2" "MIDDLESEX:yr2018:mn6:dow3"
#> [291] "MIDDLESEX:yr2018:mn6:dow4" "MIDDLESEX:yr2018:mn6:dow5"
#> [293] "MIDDLESEX:yr2018:mn6:dow6" "MIDDLESEX:yr2018:mn6:dow7"
#> [295] "MIDDLESEX:yr2018:mn7:dow1" "MIDDLESEX:yr2018:mn7:dow2"
#> [297] "MIDDLESEX:yr2018:mn7:dow3" "MIDDLESEX:yr2018:mn7:dow4"
#> [299] "MIDDLESEX:yr2018:mn7:dow5" "MIDDLESEX:yr2018:mn7:dow6"
#> [301] "MIDDLESEX:yr2018:mn7:dow7" "MIDDLESEX:yr2018:mn8:dow1"
#> [303] "MIDDLESEX:yr2018:mn8:dow2" "MIDDLESEX:yr2018:mn8:dow3"
#> [305] "MIDDLESEX:yr2018:mn8:dow4" "MIDDLESEX:yr2018:mn8:dow5"
#> [307] "MIDDLESEX:yr2018:mn8:dow6" "MIDDLESEX:yr2018:mn8:dow7"
#> [309] "MIDDLESEX:yr2018:mn9:dow1" "MIDDLESEX:yr2018:mn9:dow2"
#> [311] "MIDDLESEX:yr2018:mn9:dow3" "MIDDLESEX:yr2018:mn9:dow4"
#> [313] "MIDDLESEX:yr2018:mn9:dow5" "MIDDLESEX:yr2018:mn9:dow6"
#> [315] "MIDDLESEX:yr2018:mn9:dow7" "MIDDLESEX:yr2019:mn5:dow1"
#> [317] "MIDDLESEX:yr2019:mn5:dow2" "MIDDLESEX:yr2019:mn5:dow3"
#> [319] "MIDDLESEX:yr2019:mn5:dow4" "MIDDLESEX:yr2019:mn5:dow5"
#> [321] "MIDDLESEX:yr2019:mn5:dow6" "MIDDLESEX:yr2019:mn5:dow7"
#> [323] "MIDDLESEX:yr2019:mn6:dow1" "MIDDLESEX:yr2019:mn6:dow2"
#> [325] "MIDDLESEX:yr2019:mn6:dow3" "MIDDLESEX:yr2019:mn6:dow4"
#> [327] "MIDDLESEX:yr2019:mn6:dow5" "MIDDLESEX:yr2019:mn6:dow6"
#> [329] "MIDDLESEX:yr2019:mn6:dow7" "MIDDLESEX:yr2019:mn7:dow1"
#> [331] "MIDDLESEX:yr2019:mn7:dow2" "MIDDLESEX:yr2019:mn7:dow3"
#> [333] "MIDDLESEX:yr2019:mn7:dow4" "MIDDLESEX:yr2019:mn7:dow5"
#> [335] "MIDDLESEX:yr2019:mn7:dow6" "MIDDLESEX:yr2019:mn7:dow7"
#> [337] "MIDDLESEX:yr2019:mn8:dow1" "MIDDLESEX:yr2019:mn8:dow2"
#> [339] "MIDDLESEX:yr2019:mn8:dow3" "MIDDLESEX:yr2019:mn8:dow4"
#> [341] "MIDDLESEX:yr2019:mn8:dow5" "MIDDLESEX:yr2019:mn8:dow6"
#> [343] "MIDDLESEX:yr2019:mn8:dow7" "MIDDLESEX:yr2019:mn9:dow1"
#> [345] "MIDDLESEX:yr2019:mn9:dow2" "MIDDLESEX:yr2019:mn9:dow3"
#> [347] "MIDDLESEX:yr2019:mn9:dow4" "MIDDLESEX:yr2019:mn9:dow5"
#> [349] "MIDDLESEX:yr2019:mn9:dow6" "MIDDLESEX:yr2019:mn9:dow7"
#> [351] "MIDDLESEX:yr2020:mn5:dow1" "MIDDLESEX:yr2020:mn5:dow2"
#> [353] "MIDDLESEX:yr2020:mn5:dow3" "MIDDLESEX:yr2020:mn5:dow4"
#> [355] "MIDDLESEX:yr2020:mn5:dow5" "MIDDLESEX:yr2020:mn5:dow6"
#> [357] "MIDDLESEX:yr2020:mn5:dow7" "MIDDLESEX:yr2020:mn6:dow1"
#> [359] "MIDDLESEX:yr2020:mn6:dow2" "MIDDLESEX:yr2020:mn6:dow3"
#> [361] "MIDDLESEX:yr2020:mn6:dow4" "MIDDLESEX:yr2020:mn6:dow5"
#> [363] "MIDDLESEX:yr2020:mn6:dow6" "MIDDLESEX:yr2020:mn6:dow7"
#> [365] "MIDDLESEX:yr2020:mn7:dow1" "MIDDLESEX:yr2020:mn7:dow2"
#> [367] "MIDDLESEX:yr2020:mn7:dow3" "MIDDLESEX:yr2020:mn7:dow4"
#> [369] "MIDDLESEX:yr2020:mn7:dow5" "MIDDLESEX:yr2020:mn7:dow6"
#> [371] "MIDDLESEX:yr2020:mn7:dow7" "MIDDLESEX:yr2020:mn8:dow1"
#> [373] "MIDDLESEX:yr2020:mn8:dow2" "MIDDLESEX:yr2020:mn8:dow3"
#> [375] "MIDDLESEX:yr2020:mn8:dow4" "MIDDLESEX:yr2020:mn8:dow5"
#> [377] "MIDDLESEX:yr2020:mn8:dow6" "MIDDLESEX:yr2020:mn8:dow7"
#> [379] "MIDDLESEX:yr2020:mn9:dow1" "MIDDLESEX:yr2020:mn9:dow2"
#> [381] "MIDDLESEX:yr2020:mn9:dow3" "MIDDLESEX:yr2020:mn9:dow4"
#> [383] "MIDDLESEX:yr2020:mn9:dow5" "MIDDLESEX:yr2020:mn9:dow6"
#> [385] "MIDDLESEX:yr2020:mn9:dow7" "SUFFOLK:yr2010:mn5:dow1"  
#> [387] "SUFFOLK:yr2010:mn5:dow2"   "SUFFOLK:yr2010:mn5:dow3"  
#> [389] "SUFFOLK:yr2010:mn5:dow4"   "SUFFOLK:yr2010:mn5:dow5"  
#> [391] "SUFFOLK:yr2010:mn5:dow6"   "SUFFOLK:yr2010:mn5:dow7"  
#> [393] "SUFFOLK:yr2010:mn6:dow1"   "SUFFOLK:yr2010:mn6:dow2"  
#> [395] "SUFFOLK:yr2010:mn6:dow3"   "SUFFOLK:yr2010:mn6:dow4"  
#> [397] "SUFFOLK:yr2010:mn6:dow5"   "SUFFOLK:yr2010:mn6:dow6"  
#> [399] "SUFFOLK:yr2010:mn6:dow7"   "SUFFOLK:yr2010:mn7:dow1"  
#> [401] "SUFFOLK:yr2010:mn7:dow2"   "SUFFOLK:yr2010:mn7:dow3"  
#> [403] "SUFFOLK:yr2010:mn7:dow4"   "SUFFOLK:yr2010:mn7:dow5"  
#> [405] "SUFFOLK:yr2010:mn7:dow6"   "SUFFOLK:yr2010:mn7:dow7"  
#> [407] "SUFFOLK:yr2010:mn8:dow1"   "SUFFOLK:yr2010:mn8:dow2"  
#> [409] "SUFFOLK:yr2010:mn8:dow3"   "SUFFOLK:yr2010:mn8:dow4"  
#> [411] "SUFFOLK:yr2010:mn8:dow5"   "SUFFOLK:yr2010:mn8:dow6"  
#> [413] "SUFFOLK:yr2010:mn8:dow7"   "SUFFOLK:yr2010:mn9:dow1"  
#> [415] "SUFFOLK:yr2010:mn9:dow2"   "SUFFOLK:yr2010:mn9:dow3"  
#> [417] "SUFFOLK:yr2010:mn9:dow4"   "SUFFOLK:yr2010:mn9:dow5"  
#> [419] "SUFFOLK:yr2010:mn9:dow6"   "SUFFOLK:yr2010:mn9:dow7"  
#> [421] "SUFFOLK:yr2011:mn5:dow1"   "SUFFOLK:yr2011:mn5:dow2"  
#> [423] "SUFFOLK:yr2011:mn5:dow3"   "SUFFOLK:yr2011:mn5:dow4"  
#> [425] "SUFFOLK:yr2011:mn5:dow5"   "SUFFOLK:yr2011:mn5:dow6"  
#> [427] "SUFFOLK:yr2011:mn5:dow7"   "SUFFOLK:yr2011:mn6:dow1"  
#> [429] "SUFFOLK:yr2011:mn6:dow2"   "SUFFOLK:yr2011:mn6:dow3"  
#> [431] "SUFFOLK:yr2011:mn6:dow4"   "SUFFOLK:yr2011:mn6:dow5"  
#> [433] "SUFFOLK:yr2011:mn6:dow6"   "SUFFOLK:yr2011:mn6:dow7"  
#> [435] "SUFFOLK:yr2011:mn7:dow1"   "SUFFOLK:yr2011:mn7:dow2"  
#> [437] "SUFFOLK:yr2011:mn7:dow3"   "SUFFOLK:yr2011:mn7:dow4"  
#> [439] "SUFFOLK:yr2011:mn7:dow5"   "SUFFOLK:yr2011:mn7:dow6"  
#> [441] "SUFFOLK:yr2011:mn7:dow7"   "SUFFOLK:yr2011:mn8:dow1"  
#> [443] "SUFFOLK:yr2011:mn8:dow2"   "SUFFOLK:yr2011:mn8:dow3"  
#> [445] "SUFFOLK:yr2011:mn8:dow4"   "SUFFOLK:yr2011:mn8:dow5"  
#> [447] "SUFFOLK:yr2011:mn8:dow6"   "SUFFOLK:yr2011:mn8:dow7"  
#> [449] "SUFFOLK:yr2011:mn9:dow1"   "SUFFOLK:yr2011:mn9:dow2"  
#> [451] "SUFFOLK:yr2011:mn9:dow3"   "SUFFOLK:yr2011:mn9:dow4"  
#> [453] "SUFFOLK:yr2011:mn9:dow5"   "SUFFOLK:yr2011:mn9:dow6"  
#> [455] "SUFFOLK:yr2011:mn9:dow7"   "SUFFOLK:yr2012:mn5:dow1"  
#> [457] "SUFFOLK:yr2012:mn5:dow2"   "SUFFOLK:yr2012:mn5:dow3"  
#> [459] "SUFFOLK:yr2012:mn5:dow4"   "SUFFOLK:yr2012:mn5:dow5"  
#> [461] "SUFFOLK:yr2012:mn5:dow6"   "SUFFOLK:yr2012:mn5:dow7"  
#> [463] "SUFFOLK:yr2012:mn6:dow1"   "SUFFOLK:yr2012:mn6:dow2"  
#> [465] "SUFFOLK:yr2012:mn6:dow3"   "SUFFOLK:yr2012:mn6:dow4"  
#> [467] "SUFFOLK:yr2012:mn6:dow5"   "SUFFOLK:yr2012:mn6:dow6"  
#> [469] "SUFFOLK:yr2012:mn6:dow7"   "SUFFOLK:yr2012:mn7:dow1"  
#> [471] "SUFFOLK:yr2012:mn7:dow2"   "SUFFOLK:yr2012:mn7:dow3"  
#> [473] "SUFFOLK:yr2012:mn7:dow4"   "SUFFOLK:yr2012:mn7:dow5"  
#> [475] "SUFFOLK:yr2012:mn7:dow6"   "SUFFOLK:yr2012:mn7:dow7"  
#> [477] "SUFFOLK:yr2012:mn8:dow1"   "SUFFOLK:yr2012:mn8:dow2"  
#> [479] "SUFFOLK:yr2012:mn8:dow3"   "SUFFOLK:yr2012:mn8:dow4"  
#> [481] "SUFFOLK:yr2012:mn8:dow5"   "SUFFOLK:yr2012:mn8:dow6"  
#> [483] "SUFFOLK:yr2012:mn8:dow7"   "SUFFOLK:yr2012:mn9:dow1"  
#> [485] "SUFFOLK:yr2012:mn9:dow2"   "SUFFOLK:yr2012:mn9:dow3"  
#> [487] "SUFFOLK:yr2012:mn9:dow4"   "SUFFOLK:yr2012:mn9:dow5"  
#> [489] "SUFFOLK:yr2012:mn9:dow6"   "SUFFOLK:yr2012:mn9:dow7"  
#> [491] "SUFFOLK:yr2013:mn5:dow1"   "SUFFOLK:yr2013:mn5:dow2"  
#> [493] "SUFFOLK:yr2013:mn5:dow3"   "SUFFOLK:yr2013:mn5:dow4"  
#> [495] "SUFFOLK:yr2013:mn5:dow5"   "SUFFOLK:yr2013:mn5:dow6"  
#> [497] "SUFFOLK:yr2013:mn5:dow7"   "SUFFOLK:yr2013:mn6:dow1"  
#> [499] "SUFFOLK:yr2013:mn6:dow2"   "SUFFOLK:yr2013:mn6:dow3"  
#> [501] "SUFFOLK:yr2013:mn6:dow4"   "SUFFOLK:yr2013:mn6:dow5"  
#> [503] "SUFFOLK:yr2013:mn6:dow6"   "SUFFOLK:yr2013:mn6:dow7"  
#> [505] "SUFFOLK:yr2013:mn7:dow1"   "SUFFOLK:yr2013:mn7:dow2"  
#> [507] "SUFFOLK:yr2013:mn7:dow3"   "SUFFOLK:yr2013:mn7:dow4"  
#> [509] "SUFFOLK:yr2013:mn7:dow5"   "SUFFOLK:yr2013:mn7:dow6"  
#> [511] "SUFFOLK:yr2013:mn7:dow7"   "SUFFOLK:yr2013:mn8:dow1"  
#> [513] "SUFFOLK:yr2013:mn8:dow2"   "SUFFOLK:yr2013:mn8:dow3"  
#> [515] "SUFFOLK:yr2013:mn8:dow4"   "SUFFOLK:yr2013:mn8:dow5"  
#> [517] "SUFFOLK:yr2013:mn8:dow6"   "SUFFOLK:yr2013:mn8:dow7"  
#> [519] "SUFFOLK:yr2013:mn9:dow1"   "SUFFOLK:yr2013:mn9:dow2"  
#> [521] "SUFFOLK:yr2013:mn9:dow3"   "SUFFOLK:yr2013:mn9:dow4"  
#> [523] "SUFFOLK:yr2013:mn9:dow5"   "SUFFOLK:yr2013:mn9:dow6"  
#> [525] "SUFFOLK:yr2013:mn9:dow7"   "SUFFOLK:yr2014:mn5:dow1"  
#> [527] "SUFFOLK:yr2014:mn5:dow2"   "SUFFOLK:yr2014:mn5:dow3"  
#> [529] "SUFFOLK:yr2014:mn5:dow4"   "SUFFOLK:yr2014:mn5:dow5"  
#> [531] "SUFFOLK:yr2014:mn5:dow6"   "SUFFOLK:yr2014:mn5:dow7"  
#> [533] "SUFFOLK:yr2014:mn6:dow1"   "SUFFOLK:yr2014:mn6:dow2"  
#> [535] "SUFFOLK:yr2014:mn6:dow3"   "SUFFOLK:yr2014:mn6:dow4"  
#> [537] "SUFFOLK:yr2014:mn6:dow5"   "SUFFOLK:yr2014:mn6:dow6"  
#> [539] "SUFFOLK:yr2014:mn6:dow7"   "SUFFOLK:yr2014:mn7:dow1"  
#> [541] "SUFFOLK:yr2014:mn7:dow2"   "SUFFOLK:yr2014:mn7:dow3"  
#> [543] "SUFFOLK:yr2014:mn7:dow4"   "SUFFOLK:yr2014:mn7:dow5"  
#> [545] "SUFFOLK:yr2014:mn7:dow6"   "SUFFOLK:yr2014:mn7:dow7"  
#> [547] "SUFFOLK:yr2014:mn8:dow1"   "SUFFOLK:yr2014:mn8:dow2"  
#> [549] "SUFFOLK:yr2014:mn8:dow3"   "SUFFOLK:yr2014:mn8:dow4"  
#> [551] "SUFFOLK:yr2014:mn8:dow5"   "SUFFOLK:yr2014:mn8:dow6"  
#> [553] "SUFFOLK:yr2014:mn8:dow7"   "SUFFOLK:yr2014:mn9:dow1"  
#> [555] "SUFFOLK:yr2014:mn9:dow2"   "SUFFOLK:yr2014:mn9:dow3"  
#> [557] "SUFFOLK:yr2014:mn9:dow4"   "SUFFOLK:yr2014:mn9:dow5"  
#> [559] "SUFFOLK:yr2014:mn9:dow6"   "SUFFOLK:yr2014:mn9:dow7"  
#> [561] "SUFFOLK:yr2015:mn5:dow1"   "SUFFOLK:yr2015:mn5:dow2"  
#> [563] "SUFFOLK:yr2015:mn5:dow3"   "SUFFOLK:yr2015:mn5:dow4"  
#> [565] "SUFFOLK:yr2015:mn5:dow5"   "SUFFOLK:yr2015:mn5:dow6"  
#> [567] "SUFFOLK:yr2015:mn5:dow7"   "SUFFOLK:yr2015:mn6:dow1"  
#> [569] "SUFFOLK:yr2015:mn6:dow2"   "SUFFOLK:yr2015:mn6:dow3"  
#> [571] "SUFFOLK:yr2015:mn6:dow4"   "SUFFOLK:yr2015:mn6:dow5"  
#> [573] "SUFFOLK:yr2015:mn6:dow6"   "SUFFOLK:yr2015:mn6:dow7"  
#> [575] "SUFFOLK:yr2015:mn7:dow1"   "SUFFOLK:yr2015:mn7:dow2"  
#> [577] "SUFFOLK:yr2015:mn7:dow3"   "SUFFOLK:yr2015:mn7:dow4"  
#> [579] "SUFFOLK:yr2015:mn7:dow5"   "SUFFOLK:yr2015:mn7:dow6"  
#> [581] "SUFFOLK:yr2015:mn7:dow7"   "SUFFOLK:yr2015:mn8:dow1"  
#> [583] "SUFFOLK:yr2015:mn8:dow2"   "SUFFOLK:yr2015:mn8:dow3"  
#> [585] "SUFFOLK:yr2015:mn8:dow4"   "SUFFOLK:yr2015:mn8:dow5"  
#> [587] "SUFFOLK:yr2015:mn8:dow6"   "SUFFOLK:yr2015:mn8:dow7"  
#> [589] "SUFFOLK:yr2015:mn9:dow1"   "SUFFOLK:yr2015:mn9:dow2"  
#> [591] "SUFFOLK:yr2015:mn9:dow3"   "SUFFOLK:yr2015:mn9:dow4"  
#> [593] "SUFFOLK:yr2015:mn9:dow5"   "SUFFOLK:yr2015:mn9:dow6"  
#> [595] "SUFFOLK:yr2015:mn9:dow7"   "SUFFOLK:yr2016:mn5:dow1"  
#> [597] "SUFFOLK:yr2016:mn5:dow2"   "SUFFOLK:yr2016:mn5:dow3"  
#> [599] "SUFFOLK:yr2016:mn5:dow4"   "SUFFOLK:yr2016:mn5:dow5"  
#> [601] "SUFFOLK:yr2016:mn5:dow6"   "SUFFOLK:yr2016:mn5:dow7"  
#> [603] "SUFFOLK:yr2016:mn6:dow1"   "SUFFOLK:yr2016:mn6:dow2"  
#> [605] "SUFFOLK:yr2016:mn6:dow3"   "SUFFOLK:yr2016:mn6:dow4"  
#> [607] "SUFFOLK:yr2016:mn6:dow5"   "SUFFOLK:yr2016:mn6:dow6"  
#> [609] "SUFFOLK:yr2016:mn6:dow7"   "SUFFOLK:yr2016:mn7:dow1"  
#> [611] "SUFFOLK:yr2016:mn7:dow2"   "SUFFOLK:yr2016:mn7:dow3"  
#> [613] "SUFFOLK:yr2016:mn7:dow4"   "SUFFOLK:yr2016:mn7:dow5"  
#> [615] "SUFFOLK:yr2016:mn7:dow6"   "SUFFOLK:yr2016:mn7:dow7"  
#> [617] "SUFFOLK:yr2016:mn8:dow1"   "SUFFOLK:yr2016:mn8:dow2"  
#> [619] "SUFFOLK:yr2016:mn8:dow3"   "SUFFOLK:yr2016:mn8:dow4"  
#> [621] "SUFFOLK:yr2016:mn8:dow5"   "SUFFOLK:yr2016:mn8:dow6"  
#> [623] "SUFFOLK:yr2016:mn8:dow7"   "SUFFOLK:yr2016:mn9:dow1"  
#> [625] "SUFFOLK:yr2016:mn9:dow2"   "SUFFOLK:yr2016:mn9:dow3"  
#> [627] "SUFFOLK:yr2016:mn9:dow4"   "SUFFOLK:yr2016:mn9:dow5"  
#> [629] "SUFFOLK:yr2016:mn9:dow6"   "SUFFOLK:yr2016:mn9:dow7"  
#> [631] "SUFFOLK:yr2017:mn5:dow1"   "SUFFOLK:yr2017:mn5:dow2"  
#> [633] "SUFFOLK:yr2017:mn5:dow3"   "SUFFOLK:yr2017:mn5:dow4"  
#> [635] "SUFFOLK:yr2017:mn5:dow5"   "SUFFOLK:yr2017:mn5:dow6"  
#> [637] "SUFFOLK:yr2017:mn5:dow7"   "SUFFOLK:yr2017:mn6:dow1"  
#> [639] "SUFFOLK:yr2017:mn6:dow2"   "SUFFOLK:yr2017:mn6:dow3"  
#> [641] "SUFFOLK:yr2017:mn6:dow4"   "SUFFOLK:yr2017:mn6:dow5"  
#> [643] "SUFFOLK:yr2017:mn6:dow6"   "SUFFOLK:yr2017:mn6:dow7"  
#> [645] "SUFFOLK:yr2017:mn7:dow1"   "SUFFOLK:yr2017:mn7:dow2"  
#> [647] "SUFFOLK:yr2017:mn7:dow3"   "SUFFOLK:yr2017:mn7:dow4"  
#> [649] "SUFFOLK:yr2017:mn7:dow5"   "SUFFOLK:yr2017:mn7:dow6"  
#> [651] "SUFFOLK:yr2017:mn7:dow7"   "SUFFOLK:yr2017:mn8:dow1"  
#> [653] "SUFFOLK:yr2017:mn8:dow2"   "SUFFOLK:yr2017:mn8:dow3"  
#> [655] "SUFFOLK:yr2017:mn8:dow4"   "SUFFOLK:yr2017:mn8:dow5"  
#> [657] "SUFFOLK:yr2017:mn8:dow6"   "SUFFOLK:yr2017:mn8:dow7"  
#> [659] "SUFFOLK:yr2017:mn9:dow1"   "SUFFOLK:yr2017:mn9:dow2"  
#> [661] "SUFFOLK:yr2017:mn9:dow3"   "SUFFOLK:yr2017:mn9:dow4"  
#> [663] "SUFFOLK:yr2017:mn9:dow5"   "SUFFOLK:yr2017:mn9:dow6"  
#> [665] "SUFFOLK:yr2017:mn9:dow7"   "SUFFOLK:yr2018:mn5:dow1"  
#> [667] "SUFFOLK:yr2018:mn5:dow2"   "SUFFOLK:yr2018:mn5:dow3"  
#> [669] "SUFFOLK:yr2018:mn5:dow4"   "SUFFOLK:yr2018:mn5:dow5"  
#> [671] "SUFFOLK:yr2018:mn5:dow6"   "SUFFOLK:yr2018:mn5:dow7"  
#> [673] "SUFFOLK:yr2018:mn6:dow1"   "SUFFOLK:yr2018:mn6:dow2"  
#> [675] "SUFFOLK:yr2018:mn6:dow3"   "SUFFOLK:yr2018:mn6:dow4"  
#> [677] "SUFFOLK:yr2018:mn6:dow5"   "SUFFOLK:yr2018:mn6:dow6"  
#> [679] "SUFFOLK:yr2018:mn6:dow7"   "SUFFOLK:yr2018:mn7:dow1"  
#> [681] "SUFFOLK:yr2018:mn7:dow2"   "SUFFOLK:yr2018:mn7:dow3"  
#> [683] "SUFFOLK:yr2018:mn7:dow4"   "SUFFOLK:yr2018:mn7:dow5"  
#> [685] "SUFFOLK:yr2018:mn7:dow6"   "SUFFOLK:yr2018:mn7:dow7"  
#> [687] "SUFFOLK:yr2018:mn8:dow1"   "SUFFOLK:yr2018:mn8:dow2"  
#> [689] "SUFFOLK:yr2018:mn8:dow3"   "SUFFOLK:yr2018:mn8:dow4"  
#> [691] "SUFFOLK:yr2018:mn8:dow5"   "SUFFOLK:yr2018:mn8:dow6"  
#> [693] "SUFFOLK:yr2018:mn8:dow7"   "SUFFOLK:yr2018:mn9:dow1"  
#> [695] "SUFFOLK:yr2018:mn9:dow2"   "SUFFOLK:yr2018:mn9:dow3"  
#> [697] "SUFFOLK:yr2018:mn9:dow4"   "SUFFOLK:yr2018:mn9:dow5"  
#> [699] "SUFFOLK:yr2018:mn9:dow6"   "SUFFOLK:yr2018:mn9:dow7"  
#> [701] "SUFFOLK:yr2019:mn5:dow1"   "SUFFOLK:yr2019:mn5:dow2"  
#> [703] "SUFFOLK:yr2019:mn5:dow3"   "SUFFOLK:yr2019:mn5:dow4"  
#> [705] "SUFFOLK:yr2019:mn5:dow5"   "SUFFOLK:yr2019:mn5:dow6"  
#> [707] "SUFFOLK:yr2019:mn5:dow7"   "SUFFOLK:yr2019:mn6:dow1"  
#> [709] "SUFFOLK:yr2019:mn6:dow2"   "SUFFOLK:yr2019:mn6:dow3"  
#> [711] "SUFFOLK:yr2019:mn6:dow4"   "SUFFOLK:yr2019:mn6:dow5"  
#> [713] "SUFFOLK:yr2019:mn6:dow6"   "SUFFOLK:yr2019:mn6:dow7"  
#> [715] "SUFFOLK:yr2019:mn7:dow1"   "SUFFOLK:yr2019:mn7:dow2"  
#> [717] "SUFFOLK:yr2019:mn7:dow3"   "SUFFOLK:yr2019:mn7:dow4"  
#> [719] "SUFFOLK:yr2019:mn7:dow5"   "SUFFOLK:yr2019:mn7:dow6"  
#> [721] "SUFFOLK:yr2019:mn7:dow7"   "SUFFOLK:yr2019:mn8:dow1"  
#> [723] "SUFFOLK:yr2019:mn8:dow2"   "SUFFOLK:yr2019:mn8:dow3"  
#> [725] "SUFFOLK:yr2019:mn8:dow4"   "SUFFOLK:yr2019:mn8:dow5"  
#> [727] "SUFFOLK:yr2019:mn8:dow6"   "SUFFOLK:yr2019:mn8:dow7"  
#> [729] "SUFFOLK:yr2019:mn9:dow1"   "SUFFOLK:yr2019:mn9:dow2"  
#> [731] "SUFFOLK:yr2019:mn9:dow3"   "SUFFOLK:yr2019:mn9:dow4"  
#> [733] "SUFFOLK:yr2019:mn9:dow5"   "SUFFOLK:yr2019:mn9:dow6"  
#> [735] "SUFFOLK:yr2019:mn9:dow7"   "SUFFOLK:yr2020:mn5:dow1"  
#> [737] "SUFFOLK:yr2020:mn5:dow2"   "SUFFOLK:yr2020:mn5:dow3"  
#> [739] "SUFFOLK:yr2020:mn5:dow4"   "SUFFOLK:yr2020:mn5:dow5"  
#> [741] "SUFFOLK:yr2020:mn5:dow6"   "SUFFOLK:yr2020:mn5:dow7"  
#> [743] "SUFFOLK:yr2020:mn6:dow1"   "SUFFOLK:yr2020:mn6:dow2"  
#> [745] "SUFFOLK:yr2020:mn6:dow3"   "SUFFOLK:yr2020:mn6:dow4"  
#> [747] "SUFFOLK:yr2020:mn6:dow5"   "SUFFOLK:yr2020:mn6:dow6"  
#> [749] "SUFFOLK:yr2020:mn6:dow7"   "SUFFOLK:yr2020:mn7:dow1"  
#> [751] "SUFFOLK:yr2020:mn7:dow2"   "SUFFOLK:yr2020:mn7:dow3"  
#> [753] "SUFFOLK:yr2020:mn7:dow4"   "SUFFOLK:yr2020:mn7:dow5"  
#> [755] "SUFFOLK:yr2020:mn7:dow6"   "SUFFOLK:yr2020:mn7:dow7"  
#> [757] "SUFFOLK:yr2020:mn8:dow1"   "SUFFOLK:yr2020:mn8:dow2"  
#> [759] "SUFFOLK:yr2020:mn8:dow3"   "SUFFOLK:yr2020:mn8:dow4"  
#> [761] "SUFFOLK:yr2020:mn8:dow5"   "SUFFOLK:yr2020:mn8:dow6"  
#> [763] "SUFFOLK:yr2020:mn8:dow7"   "SUFFOLK:yr2020:mn9:dow1"  
#> [765] "SUFFOLK:yr2020:mn9:dow2"   "SUFFOLK:yr2020:mn9:dow3"  
#> [767] "SUFFOLK:yr2020:mn9:dow4"   "SUFFOLK:yr2020:mn9:dow5"  
#> [769] "SUFFOLK:yr2020:mn9:dow6"   "SUFFOLK:yr2020:mn9:dow7"

#
stopifnot(identical(deaths_tbl$date, exposure_mat$date ))
```
