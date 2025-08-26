# KLASSO

# b.tilde --> initialized mtrix of NW coeff, you must run at least one NW 
# penalty --> choosen for shrinkage; lenght(penalty) == number active var
# Y --> response var
# Z --> smooth var (single)
# X.matrix --> active variable

# starting coefficient function ----

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


