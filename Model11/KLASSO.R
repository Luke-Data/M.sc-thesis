# KLASSO

# b.tilde --> initialized mtrix of NW coeff, you must run at least one NW 
# penalty --> choosen for shrinkage; lenght(penalty) == number active var
# Y --> response var
# Z --> smooth var (single)
# X.matrix --> active variable

# DPG 

set.seed(1234)
n = 300

X1 <- rnorm(n,mean = 2, sd = 1.5)
X2 <- rbeta(n, shape1 = .5, shape2 = .5)
X3 <- rgamma(n, shape = .5, rate = 1)
Z <- runif(n, min = 0, max = 1) # Z \in [0,1]

Y <- X1*f1(Z) + X2*f2(Z) 

data <- cbind(Y,X1,X2,X3,Z)
data <- as.data.frame(data)

# set X.data active variable

X.data <- as.data.frame(data[, !names(data) %in%  c("Y","Z")])

# bw ottimale tramite cv

library(np)

bw <- npregbw(xdat = data.frame(data$Z),   # solo Z è variabile smoothing
              ydat = data$Y,
              zdat = X.data, 
              regtype = "lc",  
              bwmethod = "cv.ls")

bw.opt <- bw$bw

# starting coefficient function

b.tilde <- NW(data$Y,data$Z,X.data,h=bw.opt)

# creating a lambda penalty grid value

grid_lambda <- seq(0.01,1,length.out = 50)


KLASSO <- function(b.tilde, penalty, Z, Y, n.iter=100, X.matrix, tol = 1e-06){
  
  X.matrix <- as.matrix(X.matrix)
  
  n <- nrow(X.matrix)
  d <- ncol(X.matrix)
  
  for (i in 1:n.iter){
    L2 <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
    t.diag <- penalty / (L2 + 1e-6)
    D <- diag(t.diag)
    b_new <- matrix(NA, nrow = n, ncol = d)
    for (j in 1:NROW(Z)){
      w <- kernel(Z, Z[j], h = .4)
      W <- diag(w)
      XXW <- solve(t(X.matrix) %*% W %*% X.matrix + D)
      XYW <- t(X.matrix)%*% W %*% Y
      b_new[j,] <- XXW %*% XYW
    }
    if (sum((b_new - b.tilde)^2) < tol) {
      break
    }
    b.tilde <- b_new
  }
  return(b.tilde)
}


RSS.lambda <- function(Y, X.data, Z, b.lambda) {
  n <- length(Y)
  rss <- 0
  X.data <- as.matrix(X.data)
  for (t in 1:n) {
    weights <- kernel(Z, Z[t], .4)  
    b_t <- b.lambda[t, ]            
    residuals <- Y - X.data %*% b_t # matrice x vettore
    rss <- rss + sum((residuals^2) * weights)
  }
  return(rss / (n^2))
}


BIC.lambda <- function(RSS,b.tilde,h){
  
  n = nrow(b.tilde)
  
  df.check <- apply(b.tilde, 2, function(col) sum(col^2))
  df <- ncol(b.tilde)
  
  for (i in 1:length(df)){
    if (abs(df.check[i]) <= 0.09){
      df = df - 1
    }
  }
  return(log(RSS) + df*(log(n*h)/n*h))
}


grid_lambda <- seq(0.1,2,length.out = 20)

BIC <- rep(NA,length(grid_lambda))
RSS <- rep(NA,length(grid_lambda))
P <- matrix(data=0, ncol = ncol(X.data), nrow = length(grid_lambda))



for (i in 1:length(grid_lambda)){
  
  norma <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
  
  P[i,] <- grid_lambda[i]*(sqrt(nrow(X.data)))/norma
  
  b.tilde.temp <- KLASSO(b.tilde = b.tilde, penalty = P[i,], Y = data$Y, Z = data$Z, n.iter=100, X.matrix=X.data,tol=1e-06)
  
  RSS[i] <- RSS.lambda(data$Y,X.data,data$Z,b.tilde.temp)
  
  BIC[i] <- BIC.lambda(RSS[i], b.tilde.temp, .4)
  
}

par(mfrow=c(1,2))

plot(data$Z,b.tilde[,3], type = 'l', ylim = c(-.1,.1), xlim = c(-0.05,1))
abline(h=0)

plot(data$Z,b.tilde.temp[,3], type = 'l', ylim = c(-.1,.1),xlim = c(-0.05,1))
abline(h=0)










