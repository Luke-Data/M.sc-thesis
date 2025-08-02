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
  
  N <- matrix(data=NA,nrow = 3,ncol = 3)
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

MSE.LLR <- function(h,data,Z){
  
  X.matrix <- as.matrix(data[,c('X1','X2','X3')])
  MSE <- numeric(length(h))
  n <- NROW(data)
  
  for (i in seq_along(h)){
    
    LLR.temp <- LLR(Z,data,h[i])
    bias.temp <- LLR.temp[,c(7:9)]*(h[i]^2 / 2)
    MSE.temp <- rep(NA,length(Z))
    
    for (j in 1:length(Z)){
      
      pesi <- kernel(Z,Z[j],h[i])
      est.density <- mean(pesi)
      
      Y_fitted <- sum(X.matrix[j, ] * LLR.temp[j, 1:3])
      e2 <- (data$Y - Y_fitted)^2
      
      #e2 <- (data$Y - (LLR.temp[,1]*data$X1 + LLR.temp[,2]*data$X2 + LLR.temp[,3]*data$X3))^2
      
      N.temp <- N.e(var.l = X.matrix, pesi=pesi)
      
      inv_N <- tryCatch(solve(N.temp), error = function(e) return(NULL))
      if (is.null(inv_N)) {
        MSE.temp[j] <- NA
        next
      }
      
      asy.temp <- asy.var(n,banda=h[i], N.matrix = N.temp ,e2,est.density,pesi=pesi)
      if (is.null(asy.temp)) {
        MSE.temp[j] <- NA
        next
      }
      MSE.temp[j] <- sum(diag(N.temp%*%asy.temp))*est.density + sum(bias.temp[j,]*(N.temp %*% bias.temp[j, ]))*est.density
    }
    MSE[i] <- sum(MSE.temp,na.rm=T)
  }
  return(MSE)
}

h=seq(0.01,2,by=0.1)

mse <- MSE.LLR(h=h,data=data, Z=data$Z1)

which.min(mse)













