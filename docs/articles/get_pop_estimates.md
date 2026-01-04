# Obtaining population estimates from \`tidycensus\`

``` r

library(cityHeatHealth)
```

A common thing we have to do is obtain population estimates for
`geo_unit`s in our analysis.

We are choosing *not* to include that as part of package functions, but
it is a common task, so here is some sample code that can get you
started.

``` r


# 2) 2020 ACS 5-year (estimates; B01003_001 = total population)
library(tidycensus)

vars_acs <- c(
  # ---- Male age groups (B01001_003–B01001_025) ----
  male_under5        = "B01001_003",
  male_5_9           = "B01001_004",
  male_10_14         = "B01001_005",
  male_15_17         = "B01001_006",
  male_18_19         = "B01001_007",
  male_20            = "B01001_008",
  male_21            = "B01001_009",
  male_22_24         = "B01001_010",
  male_25_29         = "B01001_011",
  male_30_34         = "B01001_012",
  male_35_39         = "B01001_013",
  male_40_44         = "B01001_014",
  male_45_49         = "B01001_015",
  male_50_54         = "B01001_016",
  male_55_59         = "B01001_017",
  male_60_61         = "B01001_018",
  male_62_64         = "B01001_019",
  male_65_66         = "B01001_020",
  male_67_69         = "B01001_021",
  male_70_74         = "B01001_022",
  male_75_79         = "B01001_023",
  male_80_84         = "B01001_024",
  male_85_over       = "B01001_025",

  # ---- Female age groups (B01001_027–B01001_049) ----
  female_under5      = "B01001_027",
  female_5_9         = "B01001_028",
  female_10_14       = "B01001_029",
  female_15_17       = "B01001_030",
  female_18_19       = "B01001_031",
  female_20          = "B01001_032",
  female_21          = "B01001_033",
  female_22_24       = "B01001_034",
  female_25_29       = "B01001_035",
  female_30_34       = "B01001_036",
  female_35_39       = "B01001_037",
  female_40_44       = "B01001_038",
  female_45_49       = "B01001_039",
  female_50_54       = "B01001_040",
  female_55_59       = "B01001_041",
  female_60_61       = "B01001_042",
  female_62_64       = "B01001_043",
  female_65_66       = "B01001_044",
  female_67_69       = "B01001_045",
  female_70_74       = "B01001_046",
  female_75_79       = "B01001_047",
  female_80_84       = "B01001_048",
  female_85_over     = "B01001_049"

) 
```

Cities and towns in Massachusetts are under `county subdivision`

``` r

ma_cities <- get_acs(
  geography = "county subdivision",
  variables = vars_acs,
  state = "MA",
  year = 2022,
  survey = "acs5",
  geometry = FALSE
)
#> Getting data from the 2018-2022 5-year ACS
```

Now some aggregation

``` r

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
ma_cities_long <- ma_cities %>%
  select(GEOID, NAME, variable, estimate) %>%
  mutate(
    sex = ifelse(grepl("^male_", variable), "Male", "Female"),
    age_group = case_when(
      variable %in% c(
        "male_under5", "male_5_9", "male_10_14", "male_15_17",
        "female_under5", "female_5_9", "female_10_14", "female_15_17"
      ) ~ "0-17",

      variable %in% c(
        "male_18_19", "male_20", "male_21", "male_22_24",
        "male_25_29", "male_30_34", "male_35_39", "male_40_44",
        "male_45_49", "male_50_54", "male_55_59", "male_60_61",
        "male_62_64",
        "female_18_19", "female_20", "female_21", "female_22_24",
        "female_25_29", "female_30_34", "female_35_39", "female_40_44",
        "female_45_49", "female_50_54", "female_55_59", "female_60_61",
        "female_62_64"
      ) ~ "18-64",

      TRUE ~ "65+"
    )
  )

ma_pop_town_sex_age <- ma_cities_long %>%
  group_by(GEOID, NAME, sex, age_group) %>%
  summarize(
    population = sum(estimate, na.rm = TRUE),
    .groups = "drop"
  )
```

And convert to wide formatting

``` r

library(tidyr)
ma_pop_wide <- ma_pop_town_sex_age %>%
  unite(sex_age, sex, age_group) %>%
  pivot_wider(
    names_from = sex_age,
    values_from = population
  ) %>%
  select(-GEOID)

# a little cleaning
rr <- grepl("County subdivisions not defined", ma_pop_wide$NAME)
ma_pop_wide <- ma_pop_wide[!rr, ]

# and now clean names
ma_pop_wide$NAME <- toupper(ma_pop_wide$NAME)
nn <- do.call(rbind, strsplit(ma_pop_wide$NAME, ",", fixed = T))
n1 <- nn[, 1]
n2 <- gsub(" TOWN", "", n1)
n3 <- gsub(" CITY", "", n2)
ma_pop_wide$NAME <- n3

head(ma_pop_wide)
#> # A tibble: 6 × 7
#>   NAME       `Female_0-17` `Female_18-64` `Female_65+` `Male_0-17` `Male_18-64`
#>   <chr>              <dbl>          <dbl>        <dbl>       <dbl>        <dbl>
#> 1 BARNSTABLE          3899          15017         6014        4499        14035
#> 2 BOURNE              1891           5751         3212        1489         5302
#> 3 BREWSTER             634           2518         2007         833         2628
#> 4 CHATHAM              163           1477         1759         480         1265
#> 5 DENNIS               573           3792         3133         784         4101
#> 6 EASTHAM              413           1460         1263         265         1382
#> # ℹ 1 more variable: `Male_65+` <dbl>
```
