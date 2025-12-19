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
