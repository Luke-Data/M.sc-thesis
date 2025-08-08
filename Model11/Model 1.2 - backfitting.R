# Model 1.2 

kernel <- function(u, U, h) {
  
  arg <- (u - U) / h
  (1 / (h * sqrt(2 * pi))) * exp(-0.5 * arg^2)
}




kernel_qjk_all <- function(Z, X, h) {
  
  d <- ncol(Z)   # numero variabili di smoothing
  n <- nrow(Z)
  
  Q_list <- vector("list", d)
  for (j in 1:d) {
    Q_list[[j]] <- vector("list", d)
  }
  
  kernel_vals <- vector("list", d)
  for (j in 1:d) {
    kernel_vals[[j]] <- lapply(1:n, function(u) {
      kernel(Z[u, j], Z[, j], h[j])
    })
  }
  
  for (j in 1:d) {
    for (k in 1:d) {
      
      if (j == k) {
        Q_list[[j]][[k]] <- NULL
      } else {
        Q <- matrix(0, n, n)
        for (u in 1:n) {
          kj <- kernel_vals[[j]][[u]]
          for (v in 1:n) {
            kk <- kernel_vals[[k]][[v]]
            Q[u, v] <- mean(kj * kk * X[, j] * X[, k])
          }
        }
        Q_list[[j]][[k]] <- Q
      }
    }
  }
  
  return(Q_list)
}



smooth_backfitting <- function(Y, X, Z, h, max_iter = 30, tol = 1e-6) {
  
  bande <- h
  
  X <- as.matrix(X)
  Z <- as.matrix(Z) # più variabili di smoothing
  
  n <- nrow(X)
  d <- ncol(X)
  
  m_hat_attuale <- vector("list", d)
  m_tilde <- vector("list", d)
  q_j <- vector('list', d)
  
  # un passo dell'algoritmo
  
  for (j in 1:d){
    
    kernel_vals_j <- lapply(1:n, function(i) { # calcolo q_ j ottimizzato
      kernel(Z[i, j], Z[, j], h[j])
    })
    
    q_j[[j]] <- sapply(kernel_vals_j, function(w) mean(w * X[, j]^2))
    m_tilde[[j]] <- sapply(kernel_vals_j, function(w) mean(w * X[, j] * Y)) / q_j[[j]]
    
  }
  
  q_jk <- kernel_qjk_all(Z,X,h=bande) # calcolo una volta soltanto
  
  m_hat_attuale <- m_tilde
  
  for (r in 1: max_iter){
    
    m_hat_precedente <- m_hat_attuale
    
    for(j in 1:d){
      
      res <- rep(0,n)
      
      for(k in 1:d){
        
        if(k == j) next
        
        if (j > k){
          
          res <- (q_jk[[j]][[k]] %*% m_hat_attuale[[k]]) / q_j[[j]] + res  # passo r
          
        } else {
          
          res <- (q_jk[[j]][[k]] %*% m_hat_precedente[[k]]) / q_j[[j]] + res  # passo r - 1
          
        }
        
      }
      
      m_hat_attuale[[j]] <- m_tilde[[j]] - res 
      
    }
    
    if (max(sapply(1:d, function(j) max(abs(m_hat_attuale[[j]] - m_hat_precedente[[j]])))) < tol) {
      break
    }
  }
  return (m_hat_attuale)
}

