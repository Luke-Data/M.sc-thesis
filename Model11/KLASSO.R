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

b.tilde <- NW(data$Z,data,h=bw.opt)

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

plot(data$Z,b.tilde[,3], ylim = c(-.1,.1), xlim = c(-0.05,1))
abline(h=0)

plot(data$Z,b.tilde.temp[,3], ylim = c(-.1,.1),xlim = c(-0.05,1))
abline(h=0)













# DGP and visualization ----

set.seed(1234)
n = 400

X1 <- rnorm(n,mean = 2, sd = 1.5)
X2 <- rbeta(n, shape1 = .5, shape2 = .5)
X3 <- rgamma(n, shape = .5, rate = 1)
Z <- runif(n, min = 0, max = 1) # Z \in [0,1]
Z <- sort(Z)

Y <- X1*f1(Z) + X2*f2(Z) + rnorm(n)

data <- cbind(Y,X1,X2,X3,Z)
data <- as.data.frame(data)

# set X.data active variable

X.data <- as.data.frame(data[, !names(data) %in%  c("Y","Z")])

# bw ottimale tramite cv ----

library(np)

bw <- npregbw(
  xdat = data.frame(data$Z),   # Z come smoothing var
  ydat = data$Y,
  zdat = X.data,               # regressori "lineari"
  regtype = "lc",              # local constant = NW
  bwmethod = "cv.ls"           # cross-validation least squares
)

bw.opt <- bw$bw
print(bw.opt)


# starting coefficient function ----

b.tilde <- NW(data$Z,data,h=bw.opt)

par(mfrow=c(1,1))

plot(Z,f1(Z),col='indianred2', type='l',lwd=3, main=expression('True F & ' * beta[1](z)))
lines(Z,b.tilde[,1], col='skyblue', lwd=3)
legend("topright", 
       legend = c("Funzione Vera (f1)", "Funzione Stimata (b.tilde)"),
       col = c('indianred2', 'skyblue'), 
       lty = c(1, 1))

grid_lambda <- seq(0.1,100,length.out = 30) # creating a lambda penalty grid value

KLASSO <- function(b.tilde, penalty, Z, Y, n.iter=100, X.matrix, tol = 1e-06,banda = .4){
  
  X.matrix <- as.matrix(X.matrix)
  
  n <- nrow(X.matrix)
  d <- ncol(X.matrix)
  
  for (i in 1:n.iter){
    L2 <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
    t.diag <- penalty / (L2 + 1e-6)
    D <- diag(t.diag)
    b_new <- matrix(NA, nrow = n, ncol = d)
    for (j in 1:NROW(Z)){
      w <- kernel(Z, Z[j], h = bw.opt)
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


RSS.lambda <- function(Y, X.data, Z, b.lambda,h) {
  n <- length(Y)
  rss <- 0
  X.data <- as.matrix(X.data)
  for (t in 1:n) {
    weights <- kernel(Z, Z[t], h)  
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
  return(log(RSS) + df*(log(n*h)/(n*h)))
}

BIC <- rep(NA,length(grid_lambda))
RSS <- rep(NA,length(grid_lambda))
P <- matrix(data=0, ncol = ncol(X.data), nrow = length(grid_lambda))
b.tilde.temp <- vector('list', length=length(grid_lambda))


for (i in 1:length(grid_lambda)){
  
  norma <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
  
  P[i,] <- grid_lambda[i]*(sqrt(nrow(X.data)))/norma
  
  b.tilde.temp[[i]] <- KLASSO(b.tilde = b.tilde, penalty = P[i,], Y = data$Y, Z = data$Z, n.iter=100, X.matrix=X.data,tol=1e-06,banda=bw.opt)
  
  RSS[i] <- RSS.lambda(data$Y,X.data,data$Z,b.tilde.temp[[i]], bw.opt)
  
  BIC[i] <- BIC.lambda(RSS[i], b.tilde.temp[[i]], bw.opt)
  
}

par(mfrow=c(1,2))

plot(Z,b.tilde[,3], type='l', ylim = c(-1,1), xlim = c(-0.05,1))
abline(h=0)

plot(Z,b.tilde.temp[[30]][,3], type = 'l', ylim = c(-1,1),xlim = c(-0.05,1))
abline(h=0)











# 2° DGP and visualization ----

rm(list = ls())
library(MASS)

kernel <- function(u,U,h=.5) (1/(sqrt(2*pi)))*exp(-(0.5)*((u-U)/h)^2)/h

NW <- function(Z, data, h){
  Z <- data$Z
  data <- as.data.frame(data[, !names(data) %in%  c("Z")])
  coeff <- matrix(NA, ncol = ncol(data)-1, nrow = nrow(data))
  for (i in 1:nrow(data)){
    z0 <- Z[i]
    w <- kernel(Z,z0,h = h)
    coeff[i,] <- lm(Y ~ . -1, data = data, weights = w)$coefficients
  }
  return(coeff)
} 

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

data <- cbind(Y,X.dat[,1],X.dat[,2],X3,Z)
data <- as.data.frame(data)

panel.d <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.05)) 
  lines(density(x), col="#6F7598", lwd = 2) 
}

pairs(data, panel = panel.smooth, diag.panel = panel.d, col="#6da7a7")

# set X.data active variable

X.data <- as.data.frame(data[, !names(data) %in%  c("Y","Z")])

# bw ottimale tramite cv ----

library(np)

bw <- npscoefbw(Y ~ V2 + V3 | Z, data = data, bwmethod = "cv.ls")  # cross-validation
bw.opt <- bw$bw

# starting coefficient function ----

b.tilde <- NW(data$Z,data,h=bw.opt)

par(mfrow=c(1,1))

plot(Z,f1(Z),col='indianred2', type='l',lwd=3, main=expression('True F & ' * beta[1](z)))
lines(Z,b.tilde[,1], col='skyblue', lwd=3)
legend("topright", 
       legend = c("Funzione Vera (f1)", "Funzione Stimata (b.tilde)"),
       col = c('indianred2', 'skyblue'), 
       lty = c(1, 1))

grid_lambda <- seq(0.1,400,length.out = 30) # creating a lambda penalty grid value

KLASSO <- function(b.tilde, penalty, Z, Y, n.iter=100, X.matrix, tol = 1e-06,banda = .4){
  
  X.matrix <- as.matrix(X.matrix)
  
  n <- nrow(X.matrix)
  d <- ncol(X.matrix)
  
  for (i in 1:n.iter){
    L2 <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
    t.diag <- penalty / (L2 + 1e-6)
    D <- diag(t.diag)
    b_new <- matrix(NA, nrow = n, ncol = d)
    for (j in 1:NROW(Z)){
      w <- kernel(Z, Z[j], h = bw.opt)
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


RSS.lambda <- function(Y, X.data, Z, b.lambda,h) {
  n <- length(Y)
  rss <- 0
  X.data <- as.matrix(X.data)
  for (t in 1:n) {
    weights <- kernel(Z, Z[t], h)  
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
  return(log(RSS) + df*(log(n*h)/(n*h)))
}

BIC <- rep(NA,length(grid_lambda))
RSS <- rep(NA,length(grid_lambda))
P <- matrix(data=0, ncol = ncol(X.data), nrow = length(grid_lambda))
b.tilde.temp <- vector('list', length=length(grid_lambda))


for (i in 1:length(grid_lambda)){
  
  norma <- apply(b.tilde, 2, function(col) sqrt(sum(col^2)))
  
  P[i,] <- grid_lambda[i]*(sqrt(nrow(X.data)))/norma
  
  b.tilde.temp[[i]] <- KLASSO(b.tilde = b.tilde, penalty = P[i,], Y = data$Y, Z = data$Z, n.iter=100, X.matrix=X.data,tol=1e-06,banda=bw.opt)
  
  RSS[i] <- RSS.lambda(data$Y,X.data,data$Z,b.tilde.temp[[i]], bw.opt)
  
  BIC[i] <- BIC.lambda(RSS[i], b.tilde.temp[[i]], bw.opt)
  
}

par(mfrow=c(1,2))

plot(Z,b.tilde[,3], type='l', ylim = c(-1,1), xlim = c(-0.05,1))
abline(h=0)

plot(Z,b.tilde.temp[[15]][,3], type = 'l', ylim = c(-1,1),xlim = c(-0.05,1))
abline(h=0)
















