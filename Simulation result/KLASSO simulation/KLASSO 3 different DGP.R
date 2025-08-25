# KLSSO Simulation 

# DPG 1° ----

set.seed(123)
n = 400

X1 <- rnorm(n,mean = 2, sd = 1.5)
X2 <- rbeta(n, shape1 = .5, shape2 = .5)
X3 <- rgamma(n, shape = .5, rate = 1)
Z <- runif(n, min = 0, max = 1) # Z \in [0,1]
Z <- sort(Z)

f1 <- function(z) sin(z)
f2 <- function(z) 4*z*(1-z)
f3 <- function(z) cos(pi*z)

Y <- X1*f1(Z) + X2*f2(Z) + rnorm(n)

data <- cbind(Y,X1,X2,X3,Z)
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

bw <- npscoefbw(Y ~ X1 + X2 | Z, data = data, bwmethod = "cv.ls")  # cross-validation
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

grid_lambda <- seq(0.1,200,length.out = 30) # creating a lambda penalty grid value

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


library(ggplot2)

# Supponiamo che b.tilde.temp contenga le stime per ogni lambda
# Li combino in un unico data.frame lungo (long format)
plot_data <- data.frame()

for (i in seq_along(grid_lambda)) {
  df_tmp <- data.frame(
    Z = Z,
    f1 = b.tilde.temp[[i]][,1],
    f2 = b.tilde.temp[[i]][,2],
    f3 = b.tilde.temp[[i]][,3],
    Lambda = grid_lambda[i]
  )
  plot_data <- rbind(plot_data, df_tmp)
}

# porto in formato lungo per ggplot
library(tidyr)
plot_long <- pivot_longer(plot_data, cols = c("f1","f2","f3"),
                          names_to = "Function", values_to = "Value")

# True functions
truth <- data.frame(
  Z = Z,
  f1 = f1(Z),
  f2 = f2(Z),
  f3 = f3(Z)
) |> pivot_longer(cols=c("f1","f2","f3"),
                  names_to="Function", values_to="Value")

# Plot progressivo
ggplot() +
  geom_line(data = plot_long,
            aes(x = Z, y = Value, group = Lambda, color = Lambda),
            alpha = 0.7) +
  scale_color_viridis_c(option = "magma") +
  geom_line(data = truth,
            aes(x = Z, y = Value),
            color = "red", size = 1.2) +
  facet_wrap(~Function, scales = "free_y") +
  theme_minimal() +
  labs(y = "f(Z)", color = "Lambda",
       title = "Stime progressive delle coefficient functions")













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


library(ggplot2)

# Supponiamo che b.tilde.temp contenga le stime per ogni lambda
# Li combino in un unico data.frame lungo (long format)
plot_data <- data.frame()

for (i in seq_along(grid_lambda)) {
  df_tmp <- data.frame(
    Z = Z,
    f1 = b.tilde.temp[[i]][,1],
    f2 = b.tilde.temp[[i]][,2],
    f3 = b.tilde.temp[[i]][,3],
    Lambda = grid_lambda[i]
  )
  plot_data <- rbind(plot_data, df_tmp)
}

# porto in formato lungo per ggplot
library(tidyr)
plot_long <- pivot_longer(plot_data, cols = c("f1","f2","f3"),
                          names_to = "Function", values_to = "Value")

# True functions
truth <- data.frame(
  Z = Z,
  f1 = f1(Z),
  f2 = f2(Z),
  f3 = f3(Z)
) |> pivot_longer(cols=c("f1","f2","f3"),
                  names_to="Function", values_to="Value")

# Plot progressivo
ggplot() +
  geom_line(data = plot_long,
            aes(x = Z, y = Value, group = Lambda, color = Lambda),
            alpha = 0.7) +
  scale_color_viridis_c() +
  geom_line(data = truth,
            aes(x = Z, y = Value),
            color = "red", size = 1.2) +
  facet_wrap(~Function, scales = "free_y") +
  theme_minimal() +
  labs(y = "f(Z)", color = "Lambda",
       title = "Stime progressive delle coefficient functions")










# DGP 3° and visualization huge.generator() ----

rm(list = ls())
library(huge)
library(ggplot2)
library(dplyr)

n <- 200
p <- 5
set.seed(123)

df <- huge.generator(n, p, graph = 'cluster', u = .1, vis = T, g=2)
colnames(df$data) <- c('X1','X2','X3','X4','Z')
data <- as.data.frame(df$data)
Z <- data$Z

f1 <- function(z) 2*sin(2*pi*z)
f2 <- function(z) 4*z*(1-z)
f3 <- function(z) cos(pi*z)
f4 <- function(z) z^2

Y <-  f1(Z)*data$X1 + f2(Z)*data$X2 + f3(Z)*data$X3 + rnorm(n,0,1)
data <- cbind(Y,data)

data <- data[order(data$Z), ]

panel.d <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.05)) 
  lines(density(x), col="#6F7598", lwd = 2) 
}

pairs(data, panel = panel.smooth, diag.panel = panel.d, col="#6da7a7")

# NW vcm ----

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

# set X.data active variable

X.data <- as.data.frame(data[, !names(data) %in%  c("Y","Z")])

# bw ottimale tramite cv ----

library(np)

bw <- npscoefbw(Y ~ X1 + X2 + X3 | Z, data = data, bwmethod = "cv.ls")  # cross-validation
bw.opt <- bw$bw

print(bw.opt)


# starting coefficient function ----

b.tilde <- NW(data$Z,data,h=bw.opt)

par(mfrow=c(1,1))

plot(data$Z,f1(data$Z),col='indianred2', type='l',lwd=3, main=expression('True F & ' * beta[1](z)))
lines(data$Z,b.tilde[,1], col='skyblue', lwd=3)
legend("topright", 
       legend = c("Funzione Vera (f1)", "Funzione Stimata (b.tilde)"),
       col = c('indianred2', 'skyblue'), 
       lty = c(1, 1))

grid_lambda <- seq(0.1,100,length.out = 10) # creating a lambda penalty grid value

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

library(ggplot2)

plot_data <- data.frame()

for (i in seq_along(grid_lambda)) {
  df_tmp <- data.frame(
    Z = data$Z,
    f1 = b.tilde.temp[[i]][,1],
    f2 = b.tilde.temp[[i]][,2],
    f3 = b.tilde.temp[[i]][,3],
    f4 = b.tilde.temp[[i]][,4],
    Lambda = grid_lambda[i]
  )
  plot_data <- rbind(plot_data, df_tmp)
}

# porto in formato lungo per ggplot
library(tidyr)
plot_long <- pivot_longer(plot_data, cols = c("f1","f2","f3","f4"),
                          names_to = "Function", values_to = "Value")

# True functions
truth <- data.frame(
  Z = data$Z,
  f1 = f1(Z),
  f2 = f2(Z),
  f3 = f3(Z),
  f4 = f4(Z)
) |> pivot_longer(cols=c("f1","f2","f3","f4"),
                  names_to="Function", values_to="Value")

# Plot progressivo
ggplot() +
  geom_line(data = plot_long,
            aes(x = Z, y = Value, group = Lambda, color = Lambda),
            alpha = 0.7) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  scale_color_viridis_c(option = "inferno") +
  facet_wrap(~Function, scales = "free_y") +
  theme_minimal() +
  labs(y = "f(Z)", color = "Lambda",
       title = "Stime progressive delle coefficient functions")




