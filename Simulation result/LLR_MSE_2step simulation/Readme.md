# AIM: two-step estimate of local linear model with optimized bandwidth

*Workflow:*

  1) Choose optimal bandwidth through MSE minimization criterion

  2) Estimate initial coefficient function with local linear 
  
  3) Select one singular coefficient function to estimate again with diffrent bandwidth
  
  4) Compute another function estimate with new bandwidth through two_step estimate
  
  (New bandwidth calculate with MSE bandwidth optimization)
