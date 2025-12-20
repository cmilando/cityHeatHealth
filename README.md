Some things that will eventually go in the paper:

(1) The reason we are doing a conditional poisson (over time-series analysis) is because it 
has some built-in properties that can help with small numbers 
(i.e., dropping low or empty strata).

(2) The reason we are using data.table instead of tidyverse is so that joins
can be done more efficiently on smaller computers (which may be the case
for our hypothesized user)

(3) We are specifically designing this for daily data. 
As we've seen, other time frames introduce other problems that are out of scope

(4) the motivation here is for small-area studies, hence the specific setup for grp and grp-level. as per the paper these are both (i) hard to do, and (ii) essential for public health officials

(5) There as several places where numeric cut-offs are necessary: MinN, 

(6) How to deal with differing underlying population sizes? scale by the log-offset? But should it matter since you are doing within year-month-dow strata? No, the strata always have to be within populations that are the same, else the math gets messed up. I dont want to include an adjustment for log-population so i can use strata that combine multiple populations, although we could somewhat easily. Im just not sure how this would work
in stan if its expecting the Y to be a number. I suppose in a future iteration you could build in scaling and un-scaling by some nominal amount.

(7) for both exposures and outcomes, `make_xgrid` is a key function that creates the
skeleton of what you expect.

(8) throughout we also want to keep the size of various datasets small.

(9) and you need to remember its crossreduce not crosspred at the end of the day

(10) you could make this slightly more size resilient by not expanding outcomes until the ineer step, but for now this is hopefully ok
