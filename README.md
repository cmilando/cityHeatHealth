Some things that will eventually go in the paper

The reason we are doing a conditional poisson (over time-series analysis) is because it 
has some built-in properties that can help with small numbers 
(i.e., dropping low or empty strata).

The reason we are using data.table instead of tidyverse is so that joins
can be done more efficiently on smaller computers (which may be the case
for our hypothesized user)

