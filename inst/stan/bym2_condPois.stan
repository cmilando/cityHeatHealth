data {
    // First for the spatial component you'll need something similar
  int<lower=1> J; // regions, and each of the params will need a J

  // Basic params
  int<lower=1> N;               // Number of rows of data
  int<lower=1> K;               // Number of model coefficients, NO intercept
  array[N, K, J] real X;        // matrix of model parameters, NO intercept
  array[N, J] int<lower=0> y;   // outcomes

  // Now for the strata,  you need some way to signal which rows
  // are in the strata and subset to those
  // -- might be easy just to do a matrix multiplication but with 1s and 0s
  // NOTE IMPORTANT: THESE MUST BE THE SAME FOR EACH STRATA !!!
  int<lower=1> n_strata;
  int<lower=1> max_in_strata;
  array[n_strata, max_in_strata] int S_condensed;
  array[N] int<lower=0> stratum_id;

  // BYM2 adjacency info
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=J> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=J> node2; // and node1[i] < node2[i]

  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}

transformed data {


   // -------------------------------
   // https://mc-stan.org/docs/2_21/stan-users-guide/ragged-data-structs-section.html
   // first get all the strata, which include some spurious 0s
   // now, subset to just the ones that are not 0
   // updated to count backwards from the right side
   // -------------------------------

  array[n_strata] int strata_len;

  for(n in 1:n_strata) {

       int k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[n, k] > 0) k_not_zero += 1;
       }
      strata_len[n] = k_not_zero;
  }


  // -------------------------------
  // Re-mapping X into arrays, this actually makes things
  // much faster
  // -------------------------------
  array[J] matrix[N, K] Xj;

  for (j in 1:J) {
    for (n in 1:N) {
      for (k in 1:K) {
        Xj[j][n, k] = X[n, k, j];
      }
    }
  }

}

parameters {

  matrix[K, J] beta;  // attribute effects

  // BYM2 components
  // from https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/bym2_predictor_plus_offset.stan
  real<lower=0> bym2_sigma; // overall standard deviation
  real<lower=0, upper=1> bym2_rho; // proportion unstructured vs. spatially structured variance

  vector[J] bym2_theta; // heterogeneous effects
  vector[J] bym2_phi; // spatial effects


}

transformed parameters {

  // variance of each component should be approximately equal to 1

  vector[J] convolved_re;
  convolved_re = sqrt(1 - bym2_rho) * bym2_theta + sqrt(bym2_rho / scaling_factor) * bym2_phi;

  vector[J] u;
  u = convolved_re * bym2_sigma;

}

model {
  // -------------------------------
  // PRIORS
  // -------------------------------

  // Flatten beta (matrix[K, J]) to a vector and assign standard normal prior
  // This is a typical weakly informative prior for regression coefficients
  to_vector(beta) ~ std_normal();


  // See BYM2 ref above
  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(bym2_phi[node1] - bym2_phi[node2]);

  // soft sum-to-zero constraint on phi)
  sum(bym2_phi) ~ normal(0, 0.001 * J); // equivalent to mean(phi) ~ normal(0,0.001)


  bym2_theta ~ normal(0.0, 1.0);
  bym2_sigma ~ normal(0, 5);
  bym2_rho ~ beta(0.5, 0.5);


  // -------------------------------
  // LIKELIHOOD (softmax-style)
  // -------------------------------

  // Loop over strata
  for (i in 1:n_strata) {

    // idx = indices of rows belonging to this stratum
    // Only include non-zero entries from S_condensed
    array[strata_len[i]] int idx = S_condensed[i, 1:strata_len[i]];

    // Loop over regions / categories
    for (j in 1:J) {

      // Only evaluate likelihood if there are any counts in this stratum-region
      if (sum(y[idx, j]) > 0) {

        // eta_ij will store the linear predictor for each row in the stratum for this region
        vector[strata_len[i]] eta_ij;

        // Loop over rows in this stratum
        for (k in 1:strata_len[i]) {
          int n = idx[k];  // actual row index in the dataset

          // Compute linear predictor: dot product of covariates for this observation and region
          // + spatial effect for this region
          // you can remove u[j] to get normal estimates
          eta_ij[k] = dot_product(Xj[j][n], beta[, j]) + u[j];
        }

        // Add contribution to target log probability
        // multinomial_logit_lpmf automatically applies the softmax internally
        // Ensures numerical stability
        // y[idx, j] = observed counts for this stratum & region
        // eta_ij = logit-scale linear predictors for those counts
        target += multinomial_logit_lpmf(y[idx, j] | eta_ij);
      } // end if any counts in this stratum-region
    } // end loop over regions
  } // end loop over strata
}



