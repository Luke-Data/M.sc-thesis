# AIM: two-step estimate of local linear model with optimized bandwidth

*Workflow:*

  1) Estimate initial coefficient function with local linear 
  
  2) Choose optimal bandwidth through MSE minimization criterion
  
  3) Select one singular coefficient function to improove estimate
  
  4) Compute another function estimate with new bandwidth through two_step estimate
  
  (New bandwidth calculate with MSE bandwidth optimization)
