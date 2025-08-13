Comparison between estimate coefficient function using custom LLR with optimal bandwidth (using MSE criterion) vs ceooficent function obtained with mgcv package

# DGP 

library(huge)

n <- 500
p <- 4
set.seed(123456)

df <- huge.generator(n, p, graph = 'random')

colnames(df$data) <- c('X1','X2','X3','Z')

data <- as.data.frame(df$data)

# True coeff function

f1 <- function(z) exp(z)-z^2

f2 <- function(z) sin(0.5*pi*z)+4*z

f3 <- function(z) cos(pi*z)

Y <- f1(Z)*data$X1 + f2(Z)*data$X2 + f3(Z)*data$X3 + rnorm(n,0,1)
data <- cbind(Y,data)
