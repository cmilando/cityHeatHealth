# QA_testing

``` r

library(knitr)
library(kableExtra)

df <- data.frame(
  Analysis_Option = c(
    "Single ERF",
    "Single ERF",
    "Single ERF",
    "Single ERF",
    "",
    "Multiple ERFs (rf = TOWN / ZCTA)",
    "Multiple ERFs (rf = TOWN / ZCTA)",
    "Multiple ERFs (rf = TOWN / ZCTA)",
    "Multiple ERFs (rf = TOWN / ZCTA)"
  ),
  Exposure = c(
    "for 1 ZCTA",
    "for 1 TOWN\nMaintain ZCTA exposures and outcomes\ngrp_level = TRUE\nkeep_unit_exposures = TRUE",
    "Average to TOWN\ngrp_level = TRUE\nkeep_unit_exposures = FALSE",
    "",
    "",
    "For ZCTAs",
    "For TOWNS\nMaintain ZCTA exposures and outcomes\ngrp_level = TRUE\nkeep_unit_exposures = TRUE",
    "Average to TOWN\ngrp_level = TRUE\nkeep_unit_exposures = FALSE",
    ""
  ),
  Outcome = c(
    "Subset to ZCTA",
    "Subset to TOWN",
    "Subset to TOWN\ngrp_level = TRUE\nkeep_unit_outcomes = FALSE",
    "",
    "",
    "",
    "Subset to TOWN\ngrp_level = TRUE\nkeep_unit_outcomes = TRUE",
    "Subset to TOWN\ngrp_level = TRUE\nkeep_unit_outcomes = FALSE",
    ""
  ),
  Model = c(
    "condPois_1stage",
    "condPois_1stage",
    "condPois_1stage",
    "",
    "",
    "condPois_2stage",
    "condPois_2stage\n(rf = TOWN)",
    "condPois_2stage\n(rf = TOWN)",
    ""
  )
)

kable(df, align = "l", escape = FALSE) %>%
  kable_styling(full_width = FALSE, position = "left")
```

| Analysis_Option | Exposure | Outcome | Model |
|:---|:---|:---|:---|
| Single ERF | for 1 ZCTA | Subset to ZCTA | condPois_1stage |
| Single ERF | for 1 TOWN Maintain ZCTA exposures and outcomes grp_level = TRUE keep_unit_exposures = TR | E \|Subset to TOWN | \|condPois_1stage |
| Single ERF | Average to TOWN grp_level = TRUE keep_unit_exposures = FALSE | \|Subset to TOWN grp_level = TRUE keep_unit_outcomes = FA | SE \|condPois_1stage |
| Single ERF |  |  |  |
|  |  |  |  |
| Multiple ERFs (rf = TOWN / ZCTA) | For ZCTAs |  | condPois_2stage |
| Multiple ERFs (rf = TOWN / ZCTA) | For TOWNS Maintain ZCTA exposures and outcomes grp_level = TRUE keep_unit_exposures = TRU | \|Subset to TOWN grp_level = TRUE keep_unit_outcomes = T | UE \|condPois_2stage (rf = |
| Multiple ERFs (rf = TOWN / ZCTA) | Average to TOWN grp_level = TRUE keep_unit_exposures = FALSE | \|Subset to TOWN grp_level = TRUE keep_unit_outcomes = FA | SE \|condPois_2stage (rf = T |
| Multiple ERFs (rf = TOWN / ZCTA) |  |  |  |
