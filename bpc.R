bpc <- function(pgs_liab, K, prior, r2l) {
  ## Bayesian polygenic score Probability Conversion (BPC) approach.
  ## See README.md for a tutorial on how to apply this function.
  
  ## pgs_liab - PGS (or vector of PGS values) on the liability scale.
  ## K        - Disease prevalence in the population.
  ## prior    - Sample prevalence in the (hypothetical) testing sample.
  ##            (prior disease probability)
  ## r2l      - R^2 of PGS on the liability scale.

  t <- -qnorm(K, mean = 0, sd = 1) # disease threshold
  z <- dnorm(t) # height of the normal distribution at t
  i1 <- z / K; i1 # mean liability of A1 (eg Falconer and Mackay) 
  k1 <- i1 * (i1 - t); 1 - k1 # reduction in variance in A1
  mean_pgs_case <- i1 * r2l # use the known relationship between liab and pgs
  var_pgs_case <- r2l -k1 * r2l * r2l # Thalis' rule
  i0 <- -z / (1 - K) ;i0 # mean liability of A0
  k0 <- i0 * (i0 - t) ; 1 - k0 # reduction in variance in A0
  mean_pgs_control <- i0 * r2l # use the known relationship between liab and pgs
  var_pgs_control <- r2l -k0*r2l*r2l # Thalis' rule
  
  # height of the normal distribution at t
  d_case <- dnorm(pgs_liab, mean = mean_pgs_case, sd = sqrt(var_pgs_case)) 
  d_control <- dnorm(pgs_liab, mean = mean_pgs_control, sd = sqrt(var_pgs_control))
  
  pred_prob <- {prior * d_case} / {prior * d_case + (1-prior) * d_control} 
    
  return(pred_prob)
}