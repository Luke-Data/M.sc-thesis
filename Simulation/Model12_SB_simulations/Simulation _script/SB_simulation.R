library(wsbackfit)
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(grid)

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\uniform_Z\\n500\\outputs')

f1 <- function(x) 2.5 * sin(1.5*x) * exp(-0.2*x^2)
f2 <- function(x) x * cos(2*x)
f3 <- function(x) sin(1.5*x) 
f4 <- function(x) 3/(1+x**2) * cos(2*x)

# computazione in parallelo con più core

n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# MC integrated bias, integrated mse, nsim = 100 ----

n <- 500
nsim <- 100

x_grid <- seq(-3,3,length=500)

f.hat <- list(
  f11 = matrix(nrow = n, ncol = nsim),
  f12 = matrix(nrow = n, ncol = nsim),
  f32 = matrix(nrow = n, ncol = nsim),
  f31 = matrix(nrow = n, ncol = nsim)
)

f.se <- list(
  f11 = matrix(nrow = n, ncol = nsim),
  f12 = matrix(nrow = n, ncol = nsim),
  f32 = matrix(nrow = n, ncol = nsim),
  f31 = matrix(nrow = n, ncol = nsim)
)

newdata <- data.frame(
  Y  = rep(0, n),
  X1 = rep(1, n), 
  X2 = rep(1, n), 
  X3 = rep(1, n), 
  Z1 = x_grid,      
  Z2 = x_grid 
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))

sim.results <- foreach(i = 1:nsim, .packages = "wsbackfit", 
                       .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata",
                                   "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
                                     
   set.seed(i)
   
   X1 <- rnorm(n,2,3)
   X2 <- rexp(n,3)
   X3 <- rgamma(n,1,1)
   
   Z1 <- runif(n,-3,3) 
   Z2 <- runif(n,-3,3)
   
   Y <- X1 * f1(Z1) + X2 * f2(Z1) + X3 * (f3(Z2) + f4(Z1)) + rnorm(n,0,1)
   
   dati <- data.frame(Y=Y,X1 = X1, X2 = X2, X3 = X3, Z1 = Z1, Z2 = Z2)
   
   m0 <- sback(Y ~ sb(Z1,by=X1,h=-1) + sb(Z1,by=X2,h=-1) + sb(Z2,by=X3,h=-1) + sb(Z1,by=X3,h=-1),
               data=dati)
   
   pred <- predict(m0, newdata = newdata)
   
   f.hat1 <- (pred$coeff[[grep("X1:Z1|Z1:X1", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z1.*by = X1", colnames(pred$peffects))])
   
   val.hat1 <- f.hat1 - mean(f.hat1)
   
   val.se1 <- (val.hat1 - f1.true)^2
   
   f.hat2 <- (pred$coeff[[grep("X2:Z1|Z1:X2", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z1.*by = X2", colnames(pred$peffects))])
   
   val.hat2 <- f.hat2 - mean(f.hat2)
   
   val.se2 <- (val.hat2 - f2.true)^2
   
   f.hat3 <- (pred$coeff[[grep("X3:Z2|Z2:X3", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z2.*by = X3", colnames(pred$peffects))])
   
   val.hat3 <- f.hat3 - mean(f.hat3)
   
   val.se3 <- (val.hat3 - f3.true)^2
   
   f.hat4 <- ( pred$coeff[[grep("X3:Z1|Z1:X3", names(pred$coeff), value=T)]] * x_grid + 
                 pred$peffects[,grep("sb\\(Z1.*by = X3", colnames(pred$peffects))])
   
   val.hat4 <- f.hat4 - mean(f.hat4)
   
   val.se4 <- (val.hat4 - f4.true)^2
   
   list(
     hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
     se  = list(val.se1, val.se2, val.se3, val.se4)
   )
   
  }

for(i in 1:nsim) {
  res <- sim.results[[i]]
  
  for(j in 1:4) {
    f.hat[[j]][, i] <- res$hat[[j]]
    f.se[[j]][, i]  <- res$se[[j]]
  }
}

stopCluster(cl)

# integrated square bias (ibs)

ibs.11 <- mean((apply(f.hat[[1]],1,mean) - f1.true)^2)
ibs.12 <- mean((apply(f.hat[[2]],1,mean) - f2.true)^2)
ibs.32 <- mean((apply(f.hat[[3]],1,mean) - f3.true)^2)
ibs.31 <- mean((apply(f.hat[[4]],1,mean) - f4.true)^2)

ibs <- c(ibs.11,ibs.12,ibs.32,ibs.31)

# integrated MSE (imse)

imse.11 <- (max(x_grid) - min(x_grid)) * mean(apply(f.se[[1]],1,mean))
imse.12 <- (max(x_grid) - min(x_grid)) * mean(apply(f.se[[2]],1,mean))
imse.32 <- (max(x_grid) - min(x_grid)) * mean(apply(f.se[[3]],1,mean))
imse.31 <- (max(x_grid) - min(x_grid)) * mean(apply(f.se[[4]],1,mean))

imse <- c(imse.11,imse.12,imse.32,imse.31)


# MSE plot ----

plot_data <- data.frame()

labels_list <- c("f[1](Z[1])", "f[2](Z[1])", "f[3](Z[2])", "f[4](Z[1])")
y_labels <- c("MSE~(f[11])", "MSE~(f[21])", "MSE~(f[32])", "MSE~(f[41])")

for (j in 1:4) {
  mse_curve <- rowMeans(f.se[[j]], na.rm = TRUE)
  
  temp_df <- data.frame(
    x = x_grid,
    mse = mse_curve,
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}

plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()

for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = "#5C6BC0", alpha = 0.2) +
    geom_line(color = "#3949AB", linewidth = 0.5) +
    labs(
      title = parse(text = labels_list[j]),
      x = NULL,
      y = parse(text = y_labels[j])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[j]] <- p
}

combined_plot <- grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  ncol = 2,
  top = textGrob("Pointwise Mean Squared Error (MSE)", 
                 gp = gpar(fontface = "bold", fontsize = 14))
)

# squared bias plot ----

plot_data <- data.frame()

labels_list <- c("f[1](Z[1])", "f[2](Z[1])", "f[3](Z[2])", "f[4](Z[1])")
y_labels <- c("Bias~(f[11])", "Bias~(f[21])", "Bias~(f[32])", "Bias~(f[41])")

f_true_list <- list(f1.true, f2.true, f3.true, f4.true)

for (j in 1:4) {
  bias_curve <- (apply(f.hat[[j]], 1, mean) - f_true_list[[j]])^2
  
  temp_df <- data.frame(
    x = x_grid,
    bias = bias_curve,
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}

plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

plot_list <- list()

for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = "#43A047", alpha = 0.2) +
    geom_line(color = "#66BB6A", linewidth = 0.5) +
    labs(
      title = parse(text = labels_list[j]),
      x = NULL,
      y = parse(text = y_labels[j])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[j]] <- p
}

combined_plot <- grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  ncol = 2,
  top = textGrob("Pointwise Bias",
                 gp = gpar(fontface = "bold", fontsize = 14))
)


# function plot ----

par(mfrow=c(2,2), mar = c(4.5, 4.5, 2, 1), oma = c(0, 0, 3, 0))

plot(x_grid,f1(x_grid),type='l',lwd=2, xlab=expression(Z[1]),       
     ylab=expression(f[1](Z[1])))
lines(x_grid,apply(f.hat$f11,1,mean),lty='dotted', col='indianred2',lwd=1.5)
plot(x_grid,f2(x_grid),type='l',lwd=2, xlab=expression(Z[2]),
     ylab=expression(f[2](Z[2])))
lines(x_grid,apply(f.hat$f12,1,mean),lty='dotted', col='indianred2',lwd=1.5)
plot(x_grid,f3(x_grid),type='l',lwd=2,
     xlab=expression(Z[2]), ylab=expression(f[3](Z[2])))
lines(x_grid,apply(f.hat$f32,1,mean),lty='dotted', col='indianred2',lwd=1.5)
plot(x_grid,f4(x_grid),type='l',lwd=2,
     xlab=expression(Z[1]), ylab=expression(f[4](Z[1])), ylim = c(-1.19,3))
lines(x_grid,apply(f.hat$f31,1,mean),lty='dotted', col='indianred2',lwd=1.5)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("top", 
       legend = c("True function", "Estimate"), 
       col = c("black", "indianred2"), 
       lty = c("solid", "dotted"), 
       lwd = c(2, 1.5),
       horiz = TRUE,       
       bty = "n",          
       inset = c(0, 0.02), 
       cex = 1)

par(mfrow=c(1,1))


# Save results to JSON ----

library(jsonlite)

results_json <- list(
  sample_size = n,
  replications = nsim,
  f1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)

write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)