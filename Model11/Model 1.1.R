# N-W & LLR

# NW vcm ----

kernel <- function(u,U,h=.5) (1/(sqrt(2*pi)))*exp(-(0.5)*((u-U)/h)^2)/h

NW <- function(Z, data, h){
  Z <- data$Z
  data <- as.data.frame(data[, !names(data) %in%  c("Z")])
  coeff <- matrix(NA, ncol = ncol(data)-1, nrow = nrow(X))
  for (i in 1:nrow(data)){
    z0 <- Z[i]
    w <- kernel(Z,z0,h = h)
    coeff[i,] <- lm(Y ~ . -1, data = data, weights = w)$coefficients
  }
  return(coeff)
} 


# LLR vcm ----

LLP <- function(Z, data, h) {
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


# local cubic poly ----

LLP <- function(Z, data, h) {
  X <- as.matrix(data[, !names(data) %in% c("Y", "Z")])
  Y <- data$Y
  n <- nrow(X)
  p <- ncol(X)
  
  f.hat <- matrix(NA, nrow = n, ncol = 4 * p)
  
  for (i in 1:n) {
    z0 <- Z[i]
    w <- kernel(Z, z0, h = h)
    
    X_t1 <- sweep(X, 1, (Z - z0), "*")
    X_t2 <- sweep(X, 1, ((Z - z0)^2) / 2, "*")
    X_t3 <- sweep(X, 1, ((Z - z0)^3) / 6, "*")
    
    X_aug <- cbind(X, X_t1, X_t2, X_t3)
    colnames(X_aug) <- paste0(rep(colnames(X), 4), "_t", rep(0:3, each = p))
    
    X_df <- as.data.frame(X_aug)
    X_df$Y <- Y
    fit <- lm(Y ~ . -1, data = X_df, weights = w)
    
    f.hat[i, ] <- coef(fit)
  }
  
  return(f.hat)
}