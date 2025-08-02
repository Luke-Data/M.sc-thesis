# N-W & LLR

# NW vcm

kernel <- function(u,U,h=.5) (1/(sqrt(2*pi)))*exp(-(0.5)*((u-U)/h)^2)/h

NW <- function(Y, Z, data, h){
  coeff <- matrix(data = NA, nrow = nrow(data), ncol = ncol(data)) # elimino la Y e la var. di smooth
  data <- as.data.frame(data)
  for (i in 1:nrow(data)){
    z0 <- Z[i]
    w <- kernel(Z,z0,h = h)
    coeff[i,] <- lm(Y ~ . -1, data = data, weights = w)$coefficients
  }
  return(coeff)
}


# LLR vcm

##############
# TO IMPROVE #
##############

Taylor1 <- function(Z,z) Z-z
Taylor2 <- function(Z,z) ((Z-z)^2)/2

LLR <- function(Z, data, h){
  f.hat <- matrix(data = NA, nrow = NROW(data), ncol = 9)
  for (i in 1:length(Z)){
    z0 <- Z[i]
    #pesi
    w <- kernel(Z,z0,h = h)
    # matrice aumentata con serie di Taylor
    ag.data <- cbind(data[,2:4], data[,2:4] * Taylor1(Z,z0), data[,2:4] * Taylor2(Z,z0), Y = data$Y)
    colnames(ag.data)[4:9] <- c("X1t", "X2t", "X3t", "X1t2", "X2t2", "X3t2")
    
    f.hat[i,] <- lm(Y ~ X1 + X2 + X3 + X1t + X2t + X3t + X1t2 + X2t2 + X3t2 -1, data = ag.data, weights = w)$coefficients
  }
  return(f.hat)
}
