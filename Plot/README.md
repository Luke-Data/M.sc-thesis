# DPG

  set.seed(1234)
  n = 300
  
  X1 <- rnorm(n,mean = 2, sd = 1.5)
  X2 <- rbeta(n, shape1 = .5, shape2 = .5)
  X3 <- rgamma(n, shape = .5, rate = 1)
  Z1 <- runif(n, min = 0, max = 1)
  
  Y <- X1*f1(Z1) + X2*f2(Z1)

# starting value 

lambda grid --> seq(0.1,200,length.out = 20)

# PLOT

plotting some coefficient estimate function by diffrent lambda of variable X3.

RED LINE --> unpenalized function estimate
BLUE LINE --> increase penalty coefficient
