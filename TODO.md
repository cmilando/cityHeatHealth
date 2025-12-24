
## In progress: 

* Bayesian code

>> need to update to be the number of unique networks, because right now there are too many degree of freedom

--> what if you did not just n = 1 neighbors. whats about n = 2 neighbors with
a weighting factor on neighbors that are 2 away? I think thins should work. but you'll have to rework the STAN code so there are fewer 

so you'll have to have windows
and update this   // *****
  matrix[K, J] beta = rep_matrix(mu, J) + beta_star;
  // *****
  and update the y[ ] to be compared to multiple versions
  this is the way that that comes back around
and the bigger smoothness, the fewer spatial networks 

--> right because again the problem of this is it doesn't actually reduce the number of the problem. so the level of smoothing yes is determined by the user by saying how many levels of neighbors to use but 

--> ah ha, the problem there is (just like EpiEstim) you won't actually get betas
for everyplace --> but wait, you actually can if you do the poor-man's recursion. so the point would be you can have a 2-degree neighbor beta expressed in each of the J places, and then you can do just like in otherwise the ys, which would then give you region-specific betas (much like daily R(t)s). 

--> and maybe the two stage can be a weighted average?

--> probably might make sense to do look at how other people do spatial-temporal Poisson rather than doing it yourself. 



* spatial plot

I think you should make a basic one then a have additional per class
so function(x, ...) usemethod("spatial_plot") or something like that

>> spatial_plot(calcAN, shp, "rate")


* Attributable number for `condPois_1stage



* Fix the 2_stage BLUP thing ...


