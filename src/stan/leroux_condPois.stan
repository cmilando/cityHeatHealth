data {
  // First for the spatial component you'll need something similar
  int<lower=1> J; // regions, and each of the params will need a J
  matrix[J, J] Jmat; // Has to be a matrix so you can do math on it

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

}

transformed data {

  //// ********************************************************
  // PRECOMPUTE THIS BLOCK IN TRANSFORMED DATA
  // you could even pre-compute this in transformed data
  // next get n_a, which is just the sum of this row of Jmat
  // HMM DOES THIS INCLUDE ITSELF? Or is this what the 1 is for?
  // You could always add 1 if so
  row_vector[J] n_a;

  for(j in 1:J) {
    n_a[j] = sum(Jmat[j, ]);
  }
  //// ********************************************************


   //// ********************************************************
   // PRECOMPUTE THIS BLOCK IN TRANSFORMED DATA
   // will have to do something like this:
   // https://mc-stan.org/docs/2_21/stan-users-guide/ragged-data-structs-section.html
   // first get all the strata, which include some spurious 0s
   // now, subset to just the ones that are not 0
   // updated to count backwards from the right side

  array[n_strata] int strata_len;

  for(n in 1:n_strata) {

       int k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[n, k] > 0) k_not_zero += 1;
       }
      strata_len[n] = k_not_zero;
  }
  //// ********************************************************
  array[n_strata] int dummy;
  for(n in 1:n_strata) {
    dummy[n] = n;
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
  // Ok so instead, there is a single centering value for each beta
  //  and a single centering value for sigma
  vector[K] mu;
  vector<lower=0>[K] sigma; // lower limit is 1e-6 so it can NEver be 0

  // the Leroux value, just one
  real<lower=0,upper=1> q;

  // and then you need Beta* which has dimension K by J
  matrix[K, J] beta_star;  // attribute effects
}

model {

  // --------------------------------------------------------------------------
  // PRIORS
  // --------------------------------------------------------------------------
  mu ~ std_normal();
  q ~ std_normal();
  sigma ~ std_normal();

  // --------------------------------------------------------------------------
  // GET BETA STAR
  // Note: I tried and tried to get this to work with a z variable for standard
  //       but it just did not work, I think because of the sum
  //       like the star_mean depends on itself.
  // --------------------------------------------------------------------------
  // ** J is region
  // ** K is beta cofficient


  // furhter vectorized
  // So FIRST, get the sum of BETA STAR neighors that aren't the current one
  // this gets the dot product
  // beta_star needs to be WITHIN k because you are going within each coefficient
  // but across space (so across the j dimension), and then ` transpose
  // hmm - this is updated each time, which is probably not correct
  // Jmat[, j] is  x J
  // beta_star is K x J so transpose is J x K
  // to product is 1 x K;
  // it likes it in vectors so tranpose again? this might be too expensive
  // this works because beta_star has initial values because its a parameter
  matrix[K, J] beta_star_sum = beta_star * Jmat;

  // in each J a single value for the denominator applied to
  // beta star mean and st_dev
  row_vector[J] beta_star_denom = 1 - q + q * n_a;

  // Now contruct the beta star mean and sigma
  // remember to square root denom in sigma
  // since in the paper its for the variance
  // make sure to use element division and multiplication
  // broadcast as a row_vector
  matrix[K,J] star_mean = rep_matrix(q ./ beta_star_denom, K) .* beta_star_sum;

  // ok so sigma is a K vector, so assume is K x 1
  // and beta*_denom is a J row_vector, so assume its 1 x J
  // sqrt inside rep so it does fewer calculations
  matrix[K, J] star_sd = rep_matrix(sigma, J) ./ rep_matrix(sqrt(beta_star_denom), K);

  // and Now vectorized to get priors for beta_star
  to_vector(beta_star) ~ normal(to_vector(star_mean), to_vector(star_sd));

  // --------------------------------------------------------------------------
  // and get the target
  // --------------------------------------------------------------------------
  // FIRST GET BETA
  // *****
  matrix[K, J] beta = rep_matrix(mu, J) + beta_star;
  // *****

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
          eta_ij[k] = dot_product(Xj[j][n], beta[, j]);
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

generated quantities {

  // (1) make a new BETA by randomly sampling from mu and beta_star
  matrix[K, J] beta_out;

  for(k in 1:K) {
    for(j in 1:J) {
      beta_out[k,j] = mu[k] + beta_star[k,j]; // probably some additional variance here
    }
  }
}
