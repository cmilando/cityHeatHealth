The package `cityHeatHealth` makes it simple to estimate heat-health impacts at small spatial scales. Starting from a messy exposure and outcome dataset, we can quickly estimate heat-health impacts.

This package can be used in three main ways:

| 1-stage design | 2-stage design | Spatial Bayes |
|----------------|----------------|---------------|
| A 1-stage conditional poisson model when estimating a *single set of beta coefficients* for heat-health impacts across single or multiple zones | A 2-stage design is used when estimating heat-health impacts across many zones, but where *individual zone models* are desired  | If numbers are very small in the 2-stage design, spatial bayesian methods can be used to tighten confidence intervals |

Essentially, the spatial bayes (SB) method is a bridge between and 1 stage and a 2 stage
if you want more granularity than a 1 stage, but the data aren't strong enough to do a 2 stage even with blups, you can do an SB model, which is sort of like many overlapping 1 stage models of groups of spatial units fused together.

The vignettes for these are:

* `vignette("one_stage_demo")`
* `vignette("two_stage_demo")`
* `vignette("bayesian_demo")`

In implementations, an attributable number calculation is applied to model outcomes to get heat-attributable outcomes. To avoid repetition, we provide a separate vignette that details the attributable number module: `vignette("attributable_number")`. 

## Starting a new analysis

To start a new analysis, you will need the following **4** datasets:

| Exposure | Outcomes | Populations | Spatial |
|----------|----------|----------|----------|
| Exposures at the daily scale for each `geo_unit` | Health outcomes at the daily scale for each `geo_unit` | Population data for each subdivision of the health outcome data that you want results for | A map showing how the various `geo_unit`s are neighbors |

This package comes pre-loaded with sample datasets of each type (e.g., `ma_exposure`, `ma_deaths`, `ma_pop_data`, and `ma_towns` respectively) so each of these methods can be explored. See the above vignettes for starter code for each type of analysis.

## Additional vignettes

We also provide several additional vignettes for common questions:

* `vignette("comparison")`: a comparison between our methods and outputs from other tutorials
* `vignette("dlnm_sensitivity")`: shows how to assess dlnm sensitivity
* `vignette("get_pop_estimates")`: sample code for how to get population estimates from `tidycensus`
* vignette("validation")`: showing that the spatial bayesian methods work when case numbers become very small
