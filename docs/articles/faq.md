# faq

``` r

library(cityHeatHealth)
```

Some things that will eventually go in the paper:

1.  The reason we are doing a conditional poisson (over time-series
    analysis) is because it has some built-in properties that can help
    with small numbers (i.e., dropping low or empty strata). and we
    potentially trust the structure here more than the splines for time
    approach. As such there as several places where numeric cut-offs are
    necessary: MinN, strata_total. tHese are sensitivitity points

2.  The reason we are using data.table instead of tidyverse is so that
    joins can be done more efficiently on smaller computers (which may
    be the case for our hypothesized user)

3.  We are specifically designing this for daily data. I think it could
    likely be used with other types of data though, might be worth
    testing to see where it breaks. As we’ve seen, other time frames
    introduce other problems that are out of scope.

(3.5) we are also testing this for non-fatal

4.  the motivation here is for small-area studies, hence the specific
    setup for grp and grp-level. as per the paper these are both (i)
    hard to do, and (ii) essential for public health officials

5.  How to deal with differing underlying population sizes? scale by the
    log-offset? But should it matter since you are doing within
    year-month-dow strata? No, the strata always have to be within
    populations that are the same, else the math gets messed up. I dont
    want to include an adjustment for log-population so i can use strata
    that combine multiple populations, although we could somewhat
    easily. Im just not sure how this would work in stan if its
    expecting the Y to be a number. I suppose in a future iteration you
    could build in scaling and un-scaling by some nominal amount.

6.  for both exposures and outcomes, `make_xgrid` is a key function that
    creates the skeleton of what you expect.

7.  throughout we also want to keep the size of various datasets small.

8.  and you need to remember its crossreduce not crosspred at the end of
    the day

9.  you could make this slightly more size resilient by not expanding
    outcomes until the inner step, but for now this is hopefully ok with
    data.table and pointers.

10. why aren’t you doing more with future? im not sure these tasks are
    best suited, its a lot of quickly running tasks with lots of i/o and
    things to be passed around. none of these calculations actually take
    that long. I suppose you could look into it for some of the smaller
    tasks. but you’d have to make copies of some of the large objects,
    which could be bad? unless they aren’t actually that big. Also the
    speed gains are pretty small, this really doesn’t take that much
    time at all, so your total time savings will be small. The other
    thing about not doing it in future\_ is that it becomes (a) harder
    to deploy (b) harder to ensure it works on mac and pc and (c) much
    harder to debug. would be worth it if the time-savings were really a
    lot but we are talking about a few minutes here
