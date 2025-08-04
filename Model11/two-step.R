# two step ----

# s1 --> first step LLP | NW | LLR
# s2h --> new bandwidth

LLP.s2 <- function(s1, f.tilde, s2h, model){
  
  s2.tilde <- matrix(NA, nrow = nrow(data), ncol = 1)
  
  X <- as.matrix(data[, !names(data) %in% c("Y", "Z")])
  no.var <- ncol(X)
  
  s1.f <- s1[,c(1:no.var)]
  s1.f <- s1.f[,-f.tilde]
  X.temp <- X[,-f.tilde]
  
  p.res <- data$Y - rowSums(s1.f * X.temp)
  
  data.s2 <- as.data.frame(cbind(Y=p.res, X = X[,f.tilde],Z = data$Z))
  
  if (model == 'LLP'){
    fit <- LLP(Z=data$Z,data.s2,s2h)[,1] 
    return (fit)
  } 
  
  if (model == 'NW'){
    fit <- NW(Z=data$Z, data.s2,s2h)[,1]
    return (fit)
  } 
  
  if (model == 'LLR'){
    fit <- LLR(Z=data$Z, data.s2,s2h)[,1]
    return(fit)
  } 
}

