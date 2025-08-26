# MSE park ----

#######################
# Implemented for LLR #
#######################

# empirical density estimate, parzen
p.Z <- function(z,Z,banda){
  num <- kernel(z,Z,h=banda)
  return(mean(num))
}

# matrice valori attesi E(X_ j X_ k | Z = z0)
N.e <- function(var.l,pesi){
  
  p <- ncol(var.l)
  
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
  
  #LLR.t <- LLR(Z=Z,data=data,h = banda)
  #e2 <- (data$Y - (LLR.t[,1]*data$X1 + LLR.t[,2]*data$X2 + LLR.t[,3]*data$X3))^2
  #s.z <- rep(NA,1)
  
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

h=seq(0.01,2,by=0.1)

mse <- MSE.LLR(h=h,data=data, Z=data$Z1)

which.min(mse)













