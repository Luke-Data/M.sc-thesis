# Simulation result for LLR + optimal bandwidth + two step estimate

# DGP ----

set.seed(1)
n = 400

X1 <- rnorm(n,mean = 2, sd = 1.5)
X2 <- rbeta(n, shape1 = .5, shape2 = .5)
X3 <- rgamma(n, shape = .5, rate = 1)
Z <- runif(n, min = 0, max = 1) # Z \in [0,1]
Z <- sort(Z)

f1 <- function(z) sin(z)
f2 <- function(z) 4*z*(1-z)
f3 <- function(z) cos(pi*z)

Y <- X1*f1(Z) + X2*f2(Z) + X3*f3(Z) + rnorm(n)

data <- cbind(Y,X1,X2,X3,Z)
data <- as.data.frame(data)

panel.d <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.05)) 
  lines(density(x), col="#6F7598", lwd = 2) 
}

pairs(data, panel = panel.smooth, diag.panel = panel.d, col="#6da7a7")


# Import LLR function  ----

kernel <- function(u,U,h=.5) (1/(sqrt(2*pi)))*exp(-(0.5)*((u-U)/h)^2)/h

LLR <- function(Z, data, h) {
  X <- as.matrix(data[, !names(data) %in% c("Y", "Z")])
  Y <- data$Y
  n <- nrow(X)
  p <- ncol(X)
  
  f.hat <- matrix(NA, nrow = n, ncol = 3 * p)
  
  for (i in 1:n) {
    z0 <- Z[i]
    w <- kernel(Z, z0, h = h)
    
    X_t1 <- sweep(X, 1, (Z - z0), "*")
    X_t2 <- sweep(X, 1, ((Z - z0)^2) / 2, "*")
    
    X_aug <- cbind(X, X_t1, X_t2)
    colnames(X_aug) <- paste0(rep(colnames(X), 3), "_t", rep(0:2, each = p))
    
    X_df <- as.data.frame(X_aug)
    X_df$Y <- Y
    fit <- lm(Y ~ . -1, data = X_df, weights = w)
    
    f.hat[i, ] <- coef(fit)
  }
  
  return(f.hat)
}



# import MSE bandwidth selector ----

p.Z <- function(z,Z,banda){
  num <- kernel(z,Z,h=banda)
  return(mean(num))
}

# matrice valori attesi E(X_ j X_ k | Z = z0)
N.e <- function(var.l,pesi){
  
  p = ncol(var.l)
  
  N <- matrix(data=NA,nrow = p,ncol = p)
  den <- sum(pesi)
  
  for (i in 1:ncol(var.l)){
    for (j in 1:ncol(var.l)){
      
      num <- sum(pesi*var.l[,i]*var.l[,j])
      N[i,j] <- num/den
    }
  }
  return(N)
}


# varianza (Y|Z) per differenti z0 --> sigma2(z)
V <- function(e2, pesi){
  
  num <- sum(pesi*e2)
  den <- sum(pesi)
  
  s.z <- num/den
  
  return (s.z)
}

asy.var <- function(n,banda,N.matrix,e2,est.density,pesi){
  
  k <- 1/(n*banda*2*sqrt(pi))
  inv_N = solve(N.matrix)
  return (k* inv_N * V(e2,pesi)/est.density)
}



# h = vettore di bande candidate
# Z = variabile di smoothing

MSE.LLR <- function(h, data, Z) {
  
  X <- as.matrix(data[, !names(data) %in% c("Y", "Z")])
  Y <- data$Y
  n <- nrow(X)
  p <- ncol(X)
  
  MSE <- numeric(length(h))
  
  for (i in seq_along(h)) {
    
    LLR.temp <- LLR(Z, data, h[i])
    bias.temp <- LLR.temp[, (ncol(LLR.temp) - p + 1):ncol(LLR.temp), drop = FALSE] * (h[i]^2 / 2)
    MSE.temp <- numeric(n)
    
    for (j in 1:n) {
      pesi <- kernel(Z[j], Z, h[i])
      est.density <- mean(pesi)
      
      Y_fitted <- sum(X[j, ] * LLR.temp[j, 1:p])
      e2 <- (Y - Y_fitted)^2
      
      N.temp <- N.e(X, pesi)
      inv_N <- tryCatch(solve(N.temp), error = function(e) NULL)
      if (is.null(inv_N)) {
        MSE.temp[j] <- NA
        next
      }
      
      asy.temp <- tryCatch(
        asy.var(n, banda = h[i], N.matrix = N.temp, e2, est.density, pesi),
        error = function(e) NULL
      )
      if (is.null(asy.temp)) {
        MSE.temp[j] <- NA
        next
      }
      
      MSE.temp[j] <- sum(diag(N.temp %*% asy.temp)) * est.density +
        sum(bias.temp[j, ] * (N.temp %*% bias.temp[j, ])) * est.density
    }
    
    MSE[i] <- sum(MSE.temp, na.rm = TRUE)
  }
  
  return(MSE)
}

h.grid=seq(0.01,2,by=0.1)

mse <- MSE.LLR(h=h,data=data, Z=data$Z)


# optimal same bandwidth for all coefficient function thorught MSE(h) ----

opt.bw <- which.min(mse)
print(h.grid[opt.bw])

s1 <- LLR(data$Z,data,h=h.grid[opt.bw])


# estimate vs true function plot ----

par(mfrow=c(3,2), oma=c(0,0,3,0))

plot(Z, f1(Z), type = 'l', main = "Funzione vera f1(z)")
plot(Z, s1[,1], type = 'l', col='indianred', main = "Stimata f1(z)")

plot(Z, f2(Z), type = 'l', main = "Funzione vera f2(z)")
plot(Z, s1[,2], type = 'l', col='indianred', main = "Stimata f2(z)")

plot(Z, f3(Z), type = 'l', main = "Funzione vera f3(z)")
plot(Z, s1[,3], type = 'l', col='indianred', main = "Stimata f3(z)")

mtext("h = 0.21", outer = TRUE, cex = 1.5, line = 1)

par(mfrow=c(1,1))


# import two-step estimate function agumented with MSE optimal bandwidth selector inside ----


LLP.s2 <- function(s1, f.tilde, model){
  
  s2.tilde <- matrix(NA, nrow = nrow(data), ncol = 1)
  
  X <- as.matrix(data[, !names(data) %in% c("Y", "Z")])
  no.var <- ncol(X)
  
  s1.f <- s1[,c(1:no.var)]
  s1.f <- s1.f[,-f.tilde]
  X.temp <- X[,-f.tilde]
  
  p.res <- data$Y - rowSums(s1.f * X.temp)
  
  data.s2 <- as.data.frame(cbind(Y=p.res, X = X[,f.tilde],Z = data$Z))
  
  h <- MSE.LLR(h=h.grid,data.s2,Z=data$Z)
  
  bw.opt <- h.grid[which.min(h)]
  print(bw.opt)
  
  if (model == 'LLP'){
    fit <- LLP(Z=data$Z,data.s2, bw.opt)[,1] 
    return (fit)
  } 
  
  if (model == 'NW'){
    fit <- NW(Z=data$Z, data.s2,bw.opt)[,1]
    return (fit)
  } 
  
  if (model == 'LLR'){
    fit <- LLR(Z=data$Z, data.s2, bw.opt)[,1]
    return(fit)
  } 
}


# reslulting function comparison with new bandwidth ----


s2 <- LLP.s2(s1,3, model = 'LLR')


par(mfrow=c(1,2))

plot(Z,s1[,3],type = 'l', col='indianred')
lines(Z,f3(Z),type = 'l')

plot(Z,s2,type = 'l', col='indianred')
lines(Z,f3(Z),type = 'l')

par(mfrow=c(1,1))

