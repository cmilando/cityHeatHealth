# Massachusetts town-level population by sex and age group

A town-level population dataset for Massachusetts derived from the
American Community Survey (ACS) 5-year estimates. Populations are
aggregated by sex (Male, Female) and broad age groups (0–17, 18–64, 65
and older) for all Massachusetts cities and towns.

## Usage

``` r
ma_pop_data
```

## Format

A data frame with one row per town and the following variables:

- TOWN20:

  Character. Name of the Massachusetts town

- Male_0_17:

  Numeric. Estimated population of males aged 0–17 years.

- Male_18_64:

  Numeric. Estimated population of males aged 18–64 years.

- Male_65_plus:

  Numeric. Estimated population of males aged 65 years and older.

- Female_0_17:

  Numeric. Estimated population of females aged 0–17 years.

- Female_18_64:

  Numeric. Estimated population of females aged 18–64 years.

- Female_65_plus:

  Numeric. Estimated population of females aged 65 years and older.

## Source

U.S. Census Bureau, American Community Survey 5-Year Estimates.

## Details

Geographic units correspond to Census *county subdivisions (COUSUB)*,
which align with cities and towns in New England states. This dataset
includes all 351 Massachusetts cities and towns.

Population estimates are based on ACS table `B01001` (Sex by Age) and
were aggregated from Census age bins into three broad age groups
commonly used in public health analyses. Estimates represent pooled
5-year ACS values and are subject to sampling error.

This dataset is suitable for use as denominators in rate calculations,
age-stratified analyses, and small-area public health modeling in
Massachusetts.

## References

U.S. Census Bureau. American Community Survey 5-Year Data (ACS5). Sex by
Age, Table B01001.

## See also

[`get_acs`](https://walker-data.com/tidycensus/reference/get_acs.html)

## Examples

``` r
# View population totals for a single town
subset(ma_pop_data, TOWN20 == "BOSTON")
#> # A tibble: 1 × 7
#>   TOWN20 `Female_0-17` `Female_18-64` `Female_65+` `Male_0-17` `Male_18-64`
#>   <chr>          <dbl>          <dbl>        <dbl>       <dbl>        <dbl>
#> 1 BOSTON         50490         247351        47313       53784       232284
#> # ℹ 1 more variable: `Male_65+` <dbl>
```
