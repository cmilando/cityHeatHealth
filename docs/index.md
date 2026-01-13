# cityHeatHealth

The package `cityHeatHealth` makes it simple to estimate heat-health
impacts at small spatial scales. Starting from a messy exposure and
outcome dataset, we can quickly estimate heat-health impacts.

## Installing the package

    remotes::install_github("cmilando/cityHeatHealth")

## Usage

This package can be used in three main ways:

| 1-stage design | 2-stage design | Spatial Bayes |
|----|----|----|
| A 1-stage conditional poisson model when estimating a *single set of beta coefficients* for heat-health impacts across single or multiple zones: [`vignette("one_stage_demo")`](http://chadmilando.com/cityHeatHealth/articles/one_stage_demo.md) | A 2-stage design is used when estimating heat-health impacts across many zones, but where *individual zone models* are desired: [`vignette("two_stage_demo")`](http://chadmilando.com/cityHeatHealth/articles/two_stage_demo.md) | If numbers are very small in the 2-stage design, spatial bayesian methods can be used to tighten confidence intervals: [`vignette("bayesian_demo")`](http://chadmilando.com/cityHeatHealth/articles/bayesian_demo.md) |

In implementations, an attributable number calculation is applied to
model outputs, see
[`vignette("attributable_number")`](http://chadmilando.com/cityHeatHealth/articles/attributable_number.md).

## Starting a new analysis

To start a new analysis, you will need the following **4** datasets:

| Exposure | Outcomes | Populations | Spatial |
|----|----|----|----|
| Exposures at the daily scale for each `geo_unit` | Health outcomes at the daily scale for each `geo_unit` | Population data for each subdivision of the health outcome data that you want results for | A map showing how the various `geo_unit`s are neighbors |

This package comes pre-loaded with sample datasets of each type (e.g.,
`ma_exposure`, `ma_deaths`, `ma_pop_data`, and `ma_towns` respectively)
so each of these methods can be explored.

## Additional vignettes

We also provide several additional vignettes for common questions – see
“More articles””
