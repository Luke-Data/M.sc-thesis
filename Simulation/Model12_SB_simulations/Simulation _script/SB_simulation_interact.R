# SB simulation graphical model

library(foreach)
library(doParallel)
library(huge)
library(wsbackfit)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(grid)

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Dependece_structure')

n <- 500
p <- 10

set.seed(6666)

# Gaussian Graphichal model; Graph = "random") ----

Gr <- huge.generator(n, p, graph="random", prob=0.4)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

# coefficient function

f1 <- function(x) x                          
f2 <- function(x) 2 * tanh(x)                
f3 <- function(x) 2 * exp(-0.5 * x^2)    
f4 <- function(x) 0.1 * x^3
f5 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)

n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

dpg_sigma <- Gr$sigma

# MC integrated bias, integrated mse, nsim = 100 ----

n <- 500
nsim <- 100

x_grid <- seq(-3,3,length=n)

f.hat <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f38 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim),
  f5 = matrix(nrow = n, ncol = nsim),
  f36 = matrix(nrow = n, ncol = nsim)
)

f.se <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f38 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim),
  f5 = matrix(nrow = n, ncol = nsim),
  f36 = matrix(nrow = n, ncol = nsim)
)

newdata <- data.frame(
  Y  = rep(0, n),
  X1 = rep(1, n), 
  X6 = rep(1, n), 
  X10 = rep(1, n), 
  X5 = rep(1,n),
  X2 = x_grid,      
  X3 = x_grid,
  X8 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))
f5.true <- f5(x_grid) - mean(f5(x_grid))


sim.results <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
                       .export = c("f1", "f2", "f3", "f4","f5", "n", "x_grid", "newdata", "dpg_sigma","p",
                                   "f1.true", "f2.true", "f3.true", "f4.true", "f5.true")) %dopar% {
                                     
  set.seed(i)
                                     
  tryCatch({
                                                                        
    data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
    
    Y <- data$X1 * f1(data$X2) + data$X6 * f2(data$X3) + data$X10 * f3(data$X8) +
      data$X5 * f4(data$X3) + f5(data$X2) + f3(data$X6) + rnorm(n)
    
    data <- data.frame(cbind(Y,data))
    
    m0 <- sback(Y ~ sb(X2,by=X1,h=-1) + sb(X3,by=X6,h=-1) + sb(X8,by=X10,h=-1) + sb(X3,by=X5,h=-1) +
                  sb(X2,h=-1) + sb(X6,h=-1), data=data)
    
    pred <- predict(m0, newdata = newdata)
    
    f.hat1 <- (pred$coeff[[grep("X1:X2|X2:X1", names(pred$coeff), value=T)]] * x_grid + 
                 pred$peffects[,grep("sb\\(X2.*by = X1", colnames(pred$peffects))])
    
    val.hat1 <- f.hat1 - mean(f.hat1)
    val.se1  <- (val.hat1 - f1.true)^2
    
    f.hat2 <- (pred$coeff[[grep("X6:X3|X3:X6", names(pred$coeff), value=T)]] * x_grid + 
                 pred$peffects[,grep("sb\\(X3.*by = X6", colnames(pred$peffects))])
    val.hat2 <- f.hat2 - mean(f.hat2)
    val.se2  <- (val.hat2 - f2.true)^2
    
    f.hat3 <- (pred$coeff[[grep("X10:X8|X8:X10", names(pred$coeff), value=T)]] * x_grid + 
                 pred$peffects[,grep("sb\\(X8.*by = X10", colnames(pred$peffects))])
    val.hat3 <- f.hat3 - mean(f.hat3)
    val.se3  <- (val.hat3 - f3.true)^2
    
    # 4. f4 (X5 * f(X3))
    f.hat4 <- (pred$coeff[[grep("X5:X3|X3:X5", names(pred$coeff), value=T)]] * x_grid + 
                 pred$peffects[,grep("sb\\(X3.*by = X5", colnames(pred$peffects))])
    val.hat4 <- f.hat4 - mean(f.hat4)
    val.se4  <- (val.hat4 - f4.true)^2
    
    lin_X2 <- if("X2" %in% names(pred$coeff)) pred$coeff[["X2"]] else 0
    cols_X2 <- grep("sb\\(X2", colnames(pred$peffects))
    sm_X2 <- pred$peffects[, cols_X2[!grepl("by", colnames(pred$peffects)[cols_X2])]]
    
    f.hat5 <- (lin_X2 * x_grid + sm_X2)
    val.hat5 <- f.hat5 - mean(f.hat5)
    val.se5  <- (val.hat5 - f5.true)^2
    
    lin_X6 <- if("X6" %in% names(pred$coeff)) pred$coeff[["X6"]] else 0
    cols_X6 <- grep("sb\\(X6", colnames(pred$peffects))
    sm_X6 <- pred$peffects[, cols_X6[!grepl("by", colnames(pred$peffects)[cols_X6])]]
    
    f.hat6 <- (lin_X6 * x_grid + sm_X6)
    val.hat6 <- f.hat6 - mean(f.hat6)
    val.se6  <- (val.hat6 - f3.true)^2
    
    list(
      hat = list(val.hat1, val.hat2, val.hat3, val.hat4, val.hat5, val.hat6),
      se  = list(val.se1, val.se2, val.se3, val.se4, val.se5, val.se6)
    )
  }, error = function(e) NULL)
}

valid_results <- sim.results[sapply(sim.results, function(x) {
  !is.null(x) && length(x$hat[[1]]) == n
})]
nsim_valid <- length(valid_results)
cat("Simulazioni valide:", nsim_valid, "su", nsim, "\n")

# Ricrea le liste con dimensione corretta
f.hat <- vector("list", 6)
f.se <- vector("list", 6)
for(j in 1:6) {
  f.hat[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
  f.se[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
}

# Popola con i risultati validi
for(i in 1:nsim_valid) {
  for(j in 1:6) {
    f.hat[[j]][, i] <- valid_results[[i]]$hat[[j]]
    f.se[[j]][, i]  <- valid_results[[i]]$se[[j]]
  }
}
  
stopCluster(cl)

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat[[4]],1,mean) - f4.true)^2)
ibs.5 <- mean((apply(f.hat[[5]],1,mean) - f5.true)^2)
ibs.6 <- mean((apply(f.hat[[6]],1,mean) - f3.true)^2)

ibs <- c(ibs.1,ibs.2,ibs.3,ibs.4,ibs.5,ibs.6)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se[[4]],1,mean)))
imse.5 <- (6 * mean(apply(f.se[[5]],1,mean)))
imse.6 <- (6 * mean(apply(f.se[[6]],1,mean)))

imse <- c(imse.1,imse.2,imse.3,imse.4,imse.5,imse.6)

# Integrated MSE plot ----

labels_list <- c("X[1]*f[1](X[2])", "X[6]*f[2](X[3])", "X[10]*f[3](X[8])", 
                 "X[5]*f[4](X[3])", "f[5](X[2])", "f[3](X[6])")
y_labels <- c("MSE~(f[list(1,2)])", "MSE~(f[list(6,3)])", "MSE~(f[list(10,8)])", 
              "MSE~(f[list(5,3)])", "MSE~(f[list(5,2)])", "MSE~(f[list(3,6)])")
plot_data <- data.frame()
for (j in 1:6) {
  temp_df <- data.frame(
    x = x_grid,
    mse = rowMeans(f.se[[j]], na.rm = TRUE),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()
for (j in 1:6) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = "#5C6BC0", alpha = 0.2) +
    geom_line(color = "#3949AB", linewidth = 0.5) +
    labs(title = parse(text = labels_list[j]), x = NULL, y = parse(text = y_labels[j])) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  plot_list[[j]] <- p
}

grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  ncol = 3,
  top = textGrob("Pointwise Mean Squared Error (MSE)", gp = gpar(fontface = "bold", fontsize = 14))
)

# Bias ----

labels_list <- c("X[1]*f[1](X[2])", "X[6]*f[2](X[3])", "X[10]*f[3](X[8])", 
                 "X[5]*f[4](X[3])", "f[5](X[2])", "f[3](X[6])")
y_labels <- c("Bias~(f[list(1,2)])", "Bias~(f[list(6,3)])", "Bias~(f[list(10,8)])", 
              "Bias~(f[list(5,3)])", "Bias~(f[list(5,2)])", "Bias~(f[list(3,6)])")

f_true_list <- list(f1.true, f2.true, f3.true, f4.true, f5.true, f3.true)

plot_data <- data.frame()
for (j in 1:6) {
  temp_df <- data.frame(
    x = x_grid,
    bias = (apply(f.hat[[j]], 1, mean) - f_true_list[[j]]),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()
for (j in 1:6) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = "#43A047", alpha = 0.2) +
    geom_line(color = "#66BB6A", linewidth = 0.5) +
    labs(title = parse(text = labels_list[j]), x = NULL, y = parse(text = y_labels[j])) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  plot_list[[j]] <- p
}

grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  ncol = 3,
  top = textGrob("Pointwise Bias", gp = gpar(fontface = "bold", fontsize = 14))
)


# Function plot ----

par(mfrow=c(2,3), mar = c(4.5, 4.5, 2, 1), oma = c(0, 0, 3, 0))

plot(x_grid, f1(x_grid), type='l', lwd=2, xlab=expression(X[2]), ylab=expression(X[1]*f[1](X[2])))
lines(x_grid, apply(f.hat[[1]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid, f2(x_grid), type='l', lwd=2, xlab=expression(X[3]),
     ylim = c(-3,3), ylab=expression(X[6]*f[2](X[3])))
lines(x_grid, apply(f.hat[[2]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid, f3(x_grid), type='l', lwd=2, xlab=expression(X[8]),
     ylim=c(-1,2),ylab=expression(X[10]*f[3](X[8])))
lines(x_grid, apply(f.hat[[3]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid, f4(x_grid), type='l', lwd=2, xlab=expression(X[3]), ylab=expression(X[5]*f[4](X[3])))
lines(x_grid, apply(f.hat[[4]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid, f5(x_grid), type='l', lwd=2, xlab=expression(X[2]), ylab=expression(f[5](X[2])))
lines(x_grid, apply(f.hat[[5]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid, f3(x_grid), type='l', lwd=2, xlab=expression(X[6]), ylab=expression(f[3](X[6])))
lines(x_grid, apply(f.hat[[6]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", 
       legend = c("True function", "Estimate"), 
       col = c("black", "indianred2"), 
       lty = c("solid", "dotted"), 
       lwd = c(2, 1.5),
       horiz = TRUE, bty = "n", inset = c(0, 0.02), cex = 1)
par(mfrow=c(1,1))

# Save results to JSON
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim,
  f1_X1X2 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X6X3 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X10X8 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X5X3 = list(ibs = unname(ibs[4]), imse = unname(imse[4])),
  f5_X2 = list(ibs = unname(ibs[5]), imse = unname(imse[5])),
  f3_X6 = list(ibs = unname(ibs[6]), imse = unname(imse[6]))
)
write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)


##############
## n = 1000 ## ----
##############


# SB simulation graphical model

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Dependece_structure\\n1000')

n <- 1000
p <- 10

set.seed(6666)

# Gaussian Graphichal model; Graph = "random") ----

Gr <- huge.generator(n, p, graph="random", prob=0.4)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

# coefficient function

f1 <- function(x) x                          
f2 <- function(x) 2 * tanh(x)                
f3 <- function(x) 2 * exp(-0.5 * x^2)    
f4 <- function(x) 0.1 * x^3
f5 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)

n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

dpg_sigma <- Gr$sigma

# MC integrated bias, integrated mse, nsim = 100 ----

nsim <- 100

x_grid1 <- seq(-3,3,length=n)

f.hat1 <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f38 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim),
  f5 = matrix(nrow = n, ncol = nsim),
  f37 = matrix(nrow = n, ncol = nsim)
)

f.se1 <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f38 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim),
  f5 = matrix(nrow = n, ncol = nsim),
  f37 = matrix(nrow = n, ncol = nsim)
)

newdata1 <- data.frame(
  Y  = rep(0, n),
  X1 = rep(1, n), 
  X6 = rep(1, n), 
  X10 = rep(1, n), 
  X5 = rep(1,n),
  X2 = x_grid1,      
  X3 = x_grid1,
  X8 = x_grid1,
  X7 = x_grid1
)

f1.true <- f1(x_grid1) - mean(f1(x_grid1))
f2.true <- f2(x_grid1) - mean(f2(x_grid1))
f3.true <- f3(x_grid1) - mean(f3(x_grid1))
f4.true <- f4(x_grid1) - mean(f4(x_grid1))
f5.true <- f5(x_grid1) - mean(f5(x_grid1))


sim.results1 <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
  .export = c("f1", "f2", "f3", "f4","f5", "n", "x_grid1", "newdata1", "dpg_sigma","p",
  "f1.true", "f2.true", "f3.true", "f4.true", "f5.true")) %dopar% {
 
 set.seed(i)
 
 tryCatch({
   
   data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
   
   Y <- data$X1 * f1(data$X2) + data$X6 * f2(data$X3) + data$X10 * f3(data$X8) +
     data$X5 * f4(data$X3) + f5(data$X2) + f3(data$X7) + rnorm(n)
   
   data <- data.frame(cbind(Y,data))
   
   m0 <- sback(Y ~ sb(X2,by=X1,h=-1) + sb(X3,by=X6,h=-1) + sb(X8,by=X10,h=-1) + sb(X3,by=X5,h=-1) +
                 sb(X2,h=-1) + sb(X7,h=-1), data=data)
   
   pred <- predict(m0, newdata = newdata1)
   
   f.hat1 <- (pred$coeff[[grep("X1:X2|X2:X1", names(pred$coeff), value=T)]] * x_grid1 + 
                pred$peffects[,grep("sb\\(X2.*by = X1", colnames(pred$peffects))])
   
   val.hat1 <- f.hat1 - mean(f.hat1)
   val.se1  <- (val.hat1 - f1.true)^2
   
   f.hat2 <- (pred$coeff[[grep("X6:X3|X3:X6", names(pred$coeff), value=T)]] * x_grid1 + 
                pred$peffects[,grep("sb\\(X3.*by = X6", colnames(pred$peffects))])
   val.hat2 <- f.hat2 - mean(f.hat2)
   val.se2  <- (val.hat2 - f2.true)^2
   
   f.hat3 <- (pred$coeff[[grep("X10:X8|X8:X10", names(pred$coeff), value=T)]] * x_grid1 + 
                pred$peffects[,grep("sb\\(X8.*by = X10", colnames(pred$peffects))])
   val.hat3 <- f.hat3 - mean(f.hat3)
   val.se3  <- (val.hat3 - f3.true)^2
   
   # 4. f4 (X5 * f(X3))
   f.hat4 <- (pred$coeff[[grep("X5:X3|X3:X5", names(pred$coeff), value=T)]] * x_grid1 + 
                pred$peffects[,grep("sb\\(X3.*by = X5", colnames(pred$peffects))])
   val.hat4 <- f.hat4 - mean(f.hat4)
   val.se4  <- (val.hat4 - f4.true)^2
   
   lin_X2 <- if("X2" %in% names(pred$coeff)) pred$coeff[["X2"]] else 0
   cols_X2 <- grep("sb\\(X2", colnames(pred$peffects))
   sm_X2 <- pred$peffects[, cols_X2[!grepl("by", colnames(pred$peffects)[cols_X2])]]
   
   f.hat5 <- (lin_X2 * x_grid1 + sm_X2)
   val.hat5 <- f.hat5 - mean(f.hat5)
   val.se5  <- (val.hat5 - f5.true)^2
   
   lin_X7 <- if("X7" %in% names(pred$coeff)) pred$coeff[["X7"]] else 0
   cols_X7 <- grep("sb\\(X7", colnames(pred$peffects))
   sm_X7 <- pred$peffects[, cols_X7[!grepl("by", colnames(pred$peffects)[cols_X7])]]
   
   f.hat6 <- (lin_X7 * x_grid1 + sm_X7)
   val.hat6 <- f.hat6 - mean(f.hat6)
   val.se6  <- (val.hat6 - f3.true)^2
   
   list(
     hat = list(val.hat1, val.hat2, val.hat3, val.hat4, val.hat5, val.hat6),
     se  = list(val.se1, val.se2, val.se3, val.se4, val.se5, val.se6)
   )
 }, error = function(e) NULL)
}

valid_results <- sim.results1[sapply(sim.results1, function(x) {
  !is.null(x) && length(x$hat[[1]]) == n
})]
nsim_valid <- length(valid_results)
cat("Simulazioni valide:", nsim_valid, "su", nsim, "\n")

# Ricrea le liste con dimensione corretta
f.hat1 <- vector("list", 6)
f.se1 <- vector("list", 6)
for(j in 1:6) {
  f.hat1[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
  f.se1[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
}

# Popola con i risultati validi
for(i in 1:nsim_valid) {
  for(j in 1:6) {
    f.hat1[[j]][, i] <- valid_results[[i]]$hat[[j]]
    f.se1[[j]][, i]  <- valid_results[[i]]$se[[j]]
  }
}

stopCluster(cl)

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat1[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat1[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat1[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat1[[4]],1,mean) - f4.true)^2)
ibs.5 <- mean((apply(f.hat1[[5]],1,mean) - f5.true)^2)
ibs.6 <- mean((apply(f.hat1[[6]],1,mean) - f3.true)^2)

ibs <- c(ibs.1,ibs.2,ibs.3,ibs.4,ibs.5,ibs.6)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se1[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se1[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se1[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se1[[4]],1,mean)))
imse.5 <- (6 * mean(apply(f.se1[[5]],1,mean)))
imse.6 <- (6 * mean(apply(f.se1[[6]],1,mean)))

imse <- c(imse.1,imse.2,imse.3,imse.4,imse.5,imse.6)

# Integrated MSE plot ----

labels_list <- c("X[1]*f[1](X[2])", "X[6]*f[2](X[3])", "X[10]*f[3](X[8])", 
                 "X[5]*f[4](X[3])", "f[5](X[2])", "f[3](X[7])")

y_labels <- c("MSE~(f[list(1,2)])", "MSE~(f[list(6,3)])", "MSE~(f[list(10,8)])", 
              "MSE~(f[list(5,3)])", "MSE~(f[list(5,2)])", "MSE~(f[list(3,7)])")
plot_data <- data.frame()
for (j in 1:6) {
  temp_df <- data.frame(
    x = x_grid1,
    mse = rowMeans(f.se1[[j]], na.rm = TRUE),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()
for (j in 1:6) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = "#5C6BC0", alpha = 0.2) +
    geom_line(color = "#3949AB", linewidth = 0.5) +
    labs(title = parse(text = labels_list[j]), x = NULL, y = parse(text = y_labels[j])) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  plot_list[[j]] <- p
}

grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  ncol = 3,
  top = textGrob("Pointwise Mean Squared Error (MSE)", gp = gpar(fontface = "bold", fontsize = 14))
)

# Bias ----

labels_list <- c("X[1]*f[1](X[2])", "X[6]*f[2](X[3])", "X[10]*f[3](X[8])", 
                 "X[5]*f[4](X[3])", "f[5](X[2])", "f[3](X[7])")

y_labels <- c("Bias~(f[list(1,2)])", "Bias~(f[list(6,3)])", "Bias~(f[list(10,8)])", 
              "Bias~(f[list(5,3)])", "Bias~(f[list(5,2)])", "Bias~(f[list(3,7)])")

f_true_list <- list(f1.true, f2.true, f3.true, f4.true, f5.true, f3.true)

plot_data <- data.frame()
for (j in 1:6) {
  temp_df <- data.frame(
    x = x_grid1,
    bias = (apply(f.hat1[[j]], 1, mean) - f_true_list[[j]]),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()
for (j in 1:6) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = "#43A047", alpha = 0.2) +
    geom_line(color = "#66BB6A", linewidth = 0.5) +
    labs(title = parse(text = labels_list[j]), x = NULL, y = parse(text = y_labels[j])) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  plot_list[[j]] <- p
}

grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  ncol = 3,
  top = textGrob("Pointwise Bias", gp = gpar(fontface = "bold", fontsize = 14))
)

# Function plot ----

par(mfrow=c(2,3), mar = c(4.5, 4.5, 2, 1), oma = c(0, 0, 3, 0))

plot(x_grid1, f1(x_grid1), type='l', lwd=2, xlab=expression(X[2]), ylab=expression(X[1]*f[1](X[2])))
lines(x_grid1, apply(f.hat1[[1]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid1, f2(x_grid1), type='l', lwd=2, xlab=expression(X[3]),
     ylim = c(-3,3), ylab=expression(X[6]*f[2](X[3])))
lines(x_grid1, apply(f.hat1[[2]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid1, f3(x_grid1), type='l', lwd=2, xlab=expression(X[8]),
     ylim=c(-1,2),ylab=expression(X[10]*f[3](X[8])))
lines(x_grid1, apply(f.hat1[[3]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid1, f4(x_grid1), type='l', lwd=2, xlab=expression(X[3]), ylab=expression(X[5]*f[4](X[3])))
lines(x_grid1, apply(f.hat1[[4]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid1, f5(x_grid1), type='l', lwd=2, xlab=expression(X[2]), ylab=expression(f[5](X[2])))
lines(x_grid1, apply(f.hat1[[5]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

plot(x_grid1, f3(x_grid1), type='l', lwd=2, xlab=expression(X[7]), ylab=expression(f[3](X[7])))
lines(x_grid1, apply(f.hat1[[6]], 1, mean), lty='dotted', col='indianred2', lwd=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", 
       legend = c("True function", "Estimate"), 
       col = c("black", "indianred2"), 
       lty = c("solid", "dotted"), 
       lwd = c(2, 1.5),
       horiz = TRUE, bty = "n", inset = c(0, 0.02), cex = 1)
par(mfrow=c(1,1))

# Save results to JSON
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim,
  f1_X1X2 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X6X3 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X10X8 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X5X3 = list(ibs = unname(ibs[4]), imse = unname(imse[4])),
  f5_X2 = list(ibs = unname(ibs[5]), imse = unname(imse[5])),
  f3_X7 = list(ibs = unname(ibs[6]), imse = unname(imse[6]))
)
write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)




