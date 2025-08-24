# DPG 1°

  set.seed(1234)
  n = 300
  
  X1 <- rnorm(n,mean = 2, sd = 1.5)
  
  X2 <- rbeta(n, shape1 = .5, shape2 = .5)
  
  X3 <- rgamma(n, shape = .5, rate = 1)
  
  Z1 <- runif(n, min = 0, max = 1)
  
  Y <- X1 * f1(Z1) + X2 * f2(Z1)

# Lambda value 

lambda grid --> seq(0.1,200,length.out = 20)

# DPG 2°

set.seed(123)
n = 300

mu <- c(0,0)

S <- matrix(c(1,0.5^5, .5^5,1),nrow = 2,ncol = 2)

X.dat <- mvrnorm(n, mu = mu, Sigma = S)

X3 <- rbeta(n,.5,.2)

Z <- runif(n, min = 0, max = 1) # Z \in [0,1]

Z <- sort(Z)

f1 <- function(z) 2*sin(2*pi*z)

f2 <- function(z) 4*z*(1-z)

f3 <- function(z) cos(pi*z)


Y <- X.dat[,1]*f1(Z) + X.dat[,2]*f2(Z) + 3*rnorm(n)


# PLOT

plotting some coefficient estimate function by diffrent lambda of variable X3. resulting plot is shrinkage of function f(X3)

RED LINE --> unpenalized function estimate

BLUE LINE --> increase penalty coefficient
