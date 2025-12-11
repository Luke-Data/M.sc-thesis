# Smooth variable indipendent from active variable

library(foreach)
library(doParallel)
library(huge)
library(wsbackfit)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(tictoc)

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Mixed\\500_independence')

n <- 500
p <- 9

set.seed(12)

<<<<<<< HEAD
# Gaussian Graphical model ----
=======
# Gaussian Graphichal model; Graph = "random") ----

>>>>>>> ad1967eda3b9bca63243ba5971ee7043b9240169
Gr <- huge.generator(n, p, graph="random", prob=0.4)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

dpg_sigma <- Gr$sigma

# coefficient function
f1 <- function(x) x                          
f2 <- function(x) 2 * tanh(x)                
f3 <- function(x) 2 * exp(-0.5 * x^2)
f4 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)


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
  f3 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim)
)

f.se <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f3 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim)
)

newdata <- data.frame(
  Y  = rep(0, n),
  X1 = rep(1, n), 
  X6 = rep(1, n), 
  X8 = rep(1, n), 
  X5 = rep(1,n),
  Z1 = x_grid,
  Z2 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))

# option.RNG = 123, il seme in parallelo si perde 
# una iterazione per vedere se la parallelizzazione funziona

tic()
sim.results <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
           .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata", "dpg_sigma","p",
                       "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
                                     
 set.seed(i)
 
 tryCatch({
   
  data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
   
  data$Z1 <- rnorm(n)
  data$Z2 <- rnorm(n)

  Y <- data$X1 * f1(data$Z1) + data$X6 * f2(data$Z2) + data$X8 * f3(data$Z1) +
    data$X5 * f4(data$Z2) + rnorm(n)
   
  data <- data.frame(cbind(Y,data))
   
  m0 <- sback(Y ~ sb(Z1,by=X1,h=-1) + sb(Z2,by=X6,h=-1) + sb(Z1,by=X8,h=-1) + sb(Z2,by=X5,h=-1), data=data)
   
  pred <- predict(m0, newdata = newdata)
   
  f.hat1 <- (pred$coeff[[grep("X1:Z1|Z1:X1", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z1.*by = X1", colnames(pred$peffects))])
   
  val.hat1 <- f.hat1 - mean(f.hat1)
  val.se1  <- (val.hat1 - f1.true)^2
  
  f.hat2 <- (pred$coeff[[grep("X6:Z2|Z2:X6", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z2.*by = X6", colnames(pred$peffects))])
  val.hat2 <- f.hat2 - mean(f.hat2)
  val.se2  <- (val.hat2 - f2.true)^2
   
  f.hat3 <- (pred$coeff[[grep("X8:Z1|Z1:X8", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z1.*by = X8", colnames(pred$peffects))])
  val.hat3 <- f.hat3 - mean(f.hat3)
  val.se3  <- (val.hat3 - f3.true)^2
   
   # 4. f4 (X5 * f(X3))
  f.hat4 <- (pred$coeff[[grep("X5:Z2|Z2:X5", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(Z2.*by = X5", colnames(pred$peffects))])
  val.hat4 <- f.hat4 - mean(f.hat4)
  val.se4  <- (val.hat4 - f4.true)^2
   
  list(
    hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
    se  = list(val.se1, val.se2, val.se3, val.se4),
    x_data = list(data$Z1, data$Z2, data$Z1, data$Z2)
   )
 }, error = function(e) NULL)
}

stopCluster(cl)
toc()

valid_results <- sim.results[sapply(sim.results, function(x) {
  !is.null(x) && length(x$hat[[1]]) == n
})]
nsim_valid <- length(valid_results)
cat("Simulazioni valide:", nsim_valid, "su", nsim, "\n")

# Ricrea le liste con dimensione corretta
f.hat <- vector("list", 4)
f.se <- vector("list", 4)
for(j in 1:4) {
  f.hat[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
  f.se[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
}

# Popola con i risultati validi
for(i in 1:nsim_valid) {
  for(j in 1:4) {
    f.hat[[j]][, i] <- valid_results[[i]]$hat[[j]]
    f.se[[j]][, i]  <- valid_results[[i]]$se[[j]]
  }
}

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat[[4]],1,mean) - f4.true)^2)

ibs <- round(c(ibs.1, ibs.2, ibs.3, ibs.4),3)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se[[4]],1,mean)))

imse <- round(c(imse.1, imse.2, imse.3, imse.4),3)


labels_list <- c("X[1]*f[1](Z[1])", "X[6]*f[2](Z[2])", "X[8]*f[3](Z[1])", 
                 "X[5]*f[4](Z[2])")
y_labels <- c("MSE~(f[list(1,1)])", "MSE~(f[list(6,2)])", "MSE~(f[list(8,1)])", 
              "MSE~(f[list(5,2)])")

plot_data <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    mse = rowMeans(f.se[[j]], na.rm = TRUE),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

line_col <- "#2a38d3ff"
fill_col <- "#7986CB"

theme_mse <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

plot_list <- list()
for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = fill_col, alpha = 0.3) +
    geom_line(color = line_col, linewidth = 0.7) +
    coord_cartesian(ylim = c(0, 30)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels[j])) +
    theme_mse
  
  plot_list[[j]] <- p
}

combined_mse <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise MSE (GAM P-splines)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)))
  )
print(combined_mse)



# BIAS ----

y_labels_bias <- c("Bias~(f[list(1,1)])", "Bias~(f[list(6,3)])", "Bias~(f[list(8,1)])", 
                   "Bias~(f[list(5,3)])")
f_true_list <- list(f1.true, f2.true, f3.true, f4.true)

plot_data_bias <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    bias = apply(f.hat[[j]], 1, mean) - f_true_list[[j]],
    Func_Name = labels_list[j]
  )
  plot_data_bias <- rbind(plot_data_bias, temp_df)
}
plot_data_bias$Function <- factor(plot_data_bias$Func_Name, levels = labels_list)

line_col_bias <- "#2E7D32"
fill_col_bias <- "#81C784"

theme_bias <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

plot_list_bias <- list()
for (j in 1:4) {
  subset_data <- plot_data_bias[plot_data_bias$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = fill_col_bias, alpha = 0.3) +
    geom_line(color = line_col_bias, linewidth = 0.6) +
    coord_cartesian(ylim = c(-0.5, 0.7)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels_bias[j])) +
    theme_bias
  
  plot_list_bias[[j]] <- p
}

combined_bias <- (plot_list_bias[[1]] | plot_list_bias[[2]]) / 
  (plot_list_bias[[3]] | plot_list_bias[[4]]) +
  plot_annotation(
    title = "Pointwise Bias Smooth BF",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)))
  )
print(combined_bias)



f_true_list <- list(f1.true, f2.true, f3.true, f4.true)
x_labels <- c("Z[1]", "Z[2]", "Z[1]", "Z[2]")
y_labels_func <- c("X[1]*f[1](Z[1])", "X[6]*f[2](Z[2])", "X[8]*f[3](Z[1])", "X[5]*f[4](Z[2])")

true_col <- "#a6bbff"
est_col <- "#c12068ff"
ribbon_col <- "#eae0edff"

theme_func <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(size = 10, margin = margin(r = 8), face = "italic"),
    axis.title.x = element_text(size = 10, margin = margin(t = 8), face = "italic"),
    axis.text = element_text(size = 8, color = "gray40"),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#FEFEFE", color = NA),
    plot.margin = margin(12, 18, 12, 12),
    legend.position = "none"
  )

x_rug <- valid_results[[nsim_valid]]$x_data

plot_list_func <- list()
for (j in 1:4) {
  
  plot_df <- data.frame(
    x = x_grid,
    true_f = f_true_list[[j]],
    estimate = apply(f.hat[[j]], 1, mean),
    lower = apply(f.hat[[j]], 1, quantile, 0.025),
    upper = apply(f.hat[[j]], 1, quantile, 0.975)
  )
  
  rug_df <- data.frame(x_rug = x_rug[[j]])
  
  p <- ggplot(plot_df, aes(x = x)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ribbon_col, alpha = 0.7) +
    geom_line(aes(y = true_f, color = "True function"), linewidth = 0.6) +
    geom_line(aes(y = estimate, color = "Estimate"), linewidth = 0.6, linetype = "dashed") +
    geom_rug(data = rug_df, aes(x = x_rug), sides = "b", color = "black", alpha = 0.5, linewidth = 0.3) +
    scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
    coord_cartesian(ylim = c(-10, 10)) +
    labs(x = parse(text = x_labels[j]), y = parse(text = y_labels_func[j])) +
    theme_func
  
  plot_list_func[[j]] <- p
}

legend_plot <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_line(aes(color = "True function"), linewidth = 0.6) +
  geom_line(aes(color = "Estimate"), linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, linetype = c("solid", "dashed"))))

legend_grob <- cowplot::get_legend(legend_plot)

combined_func <- (plot_list_func[[1]] | plot_list_func[[2]]) / 
  (plot_list_func[[3]] | plot_list_func[[4]]) +
  plot_annotation(
    title = "True vs Estimated Functions (GAM P-splines)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 5)))
  )

wrap_elements(legend_grob) / combined_func + plot_layout(heights = c(0.08, 1))

# Save results to JSON ----
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim,
  f1_X1Z1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X6Z2 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X8Z1 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X5Z2 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)
write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)









# Smooth variabile not indipendent from active variable ----

# SB simulation graph

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Smooth_BF\\500_dependence')

n <- 500
p <- 9

set.seed(12)

# Gaussian Graphical model ----
Gr <- huge.generator(n, p, graph="random", prob=0.4)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

dpg_sigma <- Gr$sigma
# coefficient function

f1 <- function(x) x                          
f2 <- function(x) 2 * tanh(x)                
f3 <- function(x) 2 * exp(-0.5 * x^2)
f4 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)


n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

dpg_sigma <- Gr$sigma
# MC integrated bias, integrated mse, nsim = 100 ----

n <- 500
nsim <- 100

x_grid <- seq(-3,3,length=n)

f.hat1 <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f3 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim)
)

f.se1 <- list(
  f1 = matrix(nrow = n, ncol = nsim),
  f2 = matrix(nrow = n, ncol = nsim),
  f3 = matrix(nrow = n, ncol = nsim),
  f4 = matrix(nrow = n, ncol = nsim)
)

newdata1 <- data.frame(
  Y  = rep(0, n),
  X2 = rep(1, n), 
  X6 = rep(1, n), 
  X8 = rep(1, n), 
  X5 = rep(1,n),
  X1 = x_grid,
  X3 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))

tic()
sim.results1 <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
                       .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata1", "dpg_sigma","p",
                                   "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
                                     
   set.seed(i)
   
   tryCatch({
     
     data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
     colnames(data) <- paste0("X", 1:p)
     
     Y <- data$X2 * f1(data$X1) + data$X6 * f2(data$X3) + data$X8 * f3(data$X1) +
       data$X5 * f4(data$X3) + rnorm(n)
     
     data$Y <- Y
     
     m0 <- sback(Y ~ sb(X1,by=X2,h=-1) + sb(X3,by=X6,h=-1) + sb(X1,by=X8,h=-1) + sb(X3,by=X5,h=-1), data=data)
     
     pred <- predict(m0, newdata = newdata1)
     
     f.hat1 <- (pred$coeff[[grep("X2:X1|X1:X2", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X1.*by = X2", colnames(pred$peffects))])
     
     val.hat1 <- f.hat1 - mean(f.hat1)
     val.se1  <- (val.hat1 - f1.true)^2
     
     f.hat2 <- (pred$coeff[[grep("X6:X3|X3:X6", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X3.*by = X6", colnames(pred$peffects))])
     val.hat2 <- f.hat2 - mean(f.hat2)
     val.se2  <- (val.hat2 - f2.true)^2
     
     f.hat3 <- (pred$coeff[[grep("X8:X1|X1:X8", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X1.*by = X8", colnames(pred$peffects))])
     val.hat3 <- f.hat3 - mean(f.hat3)
     val.se3  <- (val.hat3 - f3.true)^2
     
     # 4. f4 (X5 * f(X3))
     f.hat4 <- (pred$coeff[[grep("X5:X3|X3:X5", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X3.*by = X5", colnames(pred$peffects))])
     val.hat4 <- f.hat4 - mean(f.hat4)
     val.se4  <- (val.hat4 - f4.true)^2
     
     list(
       hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
       se  = list(val.se1, val.se2, val.se3, val.se4),
       x_data = list(data$X1, data$X3, data$X1, data$X3)
     )
   }, error = function(e) NULL)
 }

stopCluster(cl)
toc()


valid_results <- sim.results1[sapply(sim.results1, function(x) {
  !is.null(x) && length(x$hat[[1]]) == n
})]
nsim_valid <- length(valid_results)
cat("Simulazioni valide:", nsim_valid, "su", nsim, "\n")

# Ricrea le liste con dimensione corretta
f.hat1 <- vector("list", 4)
f.se1 <- vector("list", 4)
for(j in 1:4) {
  f.hat1[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
  f.se1[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
}

# Popola con i risultati validi
for(i in 1:nsim_valid) {
  for(j in 1:4) {
    f.hat1[[j]][, i] <- valid_results[[i]]$hat[[j]]
    f.se1[[j]][, i]  <- valid_results[[i]]$se[[j]]
  }
}

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat1[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat1[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat1[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat1[[4]],1,mean) - f4.true)^2)

ibs <- round(c(ibs.1,ibs.2,ibs.3,ibs.4),3)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se1[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se1[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se1[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se1[[4]],1,mean)))

imse <- round(c(imse.1,imse.2,imse.3,imse.4),3)

# Integrated MSE plot ----

# MSE Plot ----
labels_list <- c("X[2]*f[1](X[1])", "X[6]*f[2](X[3])", "X[8]*f[3](X[1])", 
                 "X[5]*f[4](X[3])")
y_labels <- c("MSE~(f[list(2,1)])", "MSE~(f[list(6,3)])", "MSE~(f[list(8,1)])", 
              "MSE~(f[list(5,3)])")

plot_data <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    mse = rowMeans(f.se1[[j]], na.rm = TRUE),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

line_col <- "#2a38d3ff"
fill_col <- "#7986CB"

theme_mse <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

plot_list <- list()
for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = fill_col, alpha = 0.3) +
    geom_line(color = line_col, linewidth = 0.7) +
    coord_cartesian(ylim = c(0, 8)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels[j])) +
    theme_mse
  
  plot_list[[j]] <- p
}

combined_mse <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise MSE (Smooth Backfitting)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)))
  )
print(combined_mse)

# Bias Plot ----
y_labels_bias <- c("Bias~(f[list(2,1)])", "Bias~(f[list(6,3)])", "Bias~(f[list(8,1)])", 
                   "Bias~(f[list(5,3)])")
f_true_list <- list(f1.true, f2.true, f3.true, f4.true)

plot_data_bias <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    bias = apply(f.hat1[[j]], 1, mean) - f_true_list[[j]],
    Func_Name = labels_list[j]
  )
  plot_data_bias <- rbind(plot_data_bias, temp_df)
}
plot_data_bias$Function <- factor(plot_data_bias$Func_Name, levels = labels_list)

line_col_bias <- "#2E7D32"
fill_col_bias <- "#81C784"

theme_bias <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

plot_list_bias <- list()
for (j in 1:4) {
  subset_data <- plot_data_bias[plot_data_bias$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = fill_col_bias, alpha = 0.3) +
    geom_line(color = line_col_bias, linewidth = 0.6) +
    coord_cartesian(ylim = c(-0.6, 0.6)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels_bias[j])) +
    theme_bias
  
  plot_list_bias[[j]] <- p
}

combined_bias <- (plot_list_bias[[1]] | plot_list_bias[[2]]) / 
  (plot_list_bias[[3]] | plot_list_bias[[4]]) +
  plot_annotation(
    title = "Pointwise Bias (Smooth Backfitting)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)))
  )
print(combined_bias)


# Function Plot ----
f_true_list <- list(f1.true, f2.true, f3.true, f4.true)
x_labels <- c("X[1]", "X[3]", "X[1]", "X[3]")
y_labels_func <- c("X[2]*f[1](X[1])", "X[6]*f[2](X[3])", "X[8]*f[3](X[1])", "X[5]*f[4](X[3])")

true_col <- "#a6bbff"
est_col <- "#c12068ff"
ribbon_col <- "#eae0edff"

theme_func <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(size = 10, margin = margin(r = 8), face = "italic"),
    axis.title.x = element_text(size = 10, margin = margin(t = 8), face = "italic"),
    axis.text = element_text(size = 8, color = "gray40"),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#FEFEFE", color = NA),
    plot.margin = margin(12, 18, 12, 12),
    legend.position = "none"
  )

x_rug <- valid_results[[nsim_valid]]$x_data

plot_list_func <- list()
for (j in 1:4) {
  
  plot_df <- data.frame(
    x = x_grid,
    true_f = f_true_list[[j]],
    estimate = apply(f.hat1[[j]], 1, mean),
    lower = apply(f.hat1[[j]], 1, quantile, 0.025),
    upper = apply(f.hat1[[j]], 1, quantile, 0.975)
  )
  
  rug_df <- data.frame(x_rug = x_rug[[j]])
  
  p <- ggplot(plot_df, aes(x = x)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ribbon_col, alpha = 0.7) +
    geom_line(aes(y = true_f, color = "True function"), linewidth = 0.6) +
    geom_line(aes(y = estimate, color = "Estimate"), linewidth = 0.6, linetype = "dashed") +
    geom_rug(data = rug_df, aes(x = x_rug), sides = "b", color = "black", alpha = 0.5, linewidth = 0.3) +
    scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
    coord_cartesian(ylim = c(-6, 6)) +
    labs(x = parse(text = x_labels[j]), y = parse(text = y_labels_func[j])) +
    theme_func
  
  plot_list_func[[j]] <- p
}

legend_plot <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_line(aes(color = "True function"), linewidth = 0.6) +
  geom_line(aes(color = "Estimate"), linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, linetype = c("solid", "dashed"))))

legend_grob <- cowplot::get_legend(legend_plot)

combined_func <- (plot_list_func[[1]] | plot_list_func[[2]]) / 
  (plot_list_func[[3]] | plot_list_func[[4]]) +
  plot_annotation(
    title = "True vs Estimated Functions (Smooth Backfitting)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 5)))
  )

wrap_elements(legend_grob) / combined_func + plot_layout(heights = c(0.08, 1))

# Save results to JSON ----
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim_valid,
  method = "Smooth Backfitting (kernel)",
  f1_X2X1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X6X3 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X8X1 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X5X3 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)
write_json(results_json, "results_sb.json", pretty = TRUE, auto_unbox = TRUE)


# aggiungere variabili che non ci sono. Aggiungo una variabile che sta nel grafo ma non nel modello vero
# Permutation test 


