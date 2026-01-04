# Calculate the dispersion parameter for a quasi-poisson model

Converted from STATA from armstrong 2014 supplement with LLM guidance.
Minimal error checking of these inputs at present.

## Usage

``` r
calc_dispersion(y, X, beta, stratum_vector)
```

## Arguments

- y:

  a vector of outcomes

- X:

  a matrix of predictors, typically the crossbasis output

- beta:

  a vector of coefficients

- stratum_vector:

  a vector describing the stratum
