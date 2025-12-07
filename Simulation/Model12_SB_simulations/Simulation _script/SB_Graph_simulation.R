# SB simulation indipendent smoothing variable

library(foreach)
library(doParallel)
library(huge)
library(wsbackfit)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(patchwork)

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Mixed\\n500')

n <- 500
p <- 15

set.seed(12)

# Gaussian Graphichal model; Graph = "random") ----

Gr <- huge.generator(n, p, graph="random", prob=0.23)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

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
  X3 = x_grid,
  X9 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))

# option.RNG = 123, il seme in parallelo si perde 
# una iterazione per vedere se la parallelizzazione funziona

sim.results <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
           .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata", "dpg_sigma","p",
                       "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
                                     
 set.seed(i)
 
 tryCatch({
   
   data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
   
   Y <- data$X1 * f1(data$X9) + data$X6 * f2(data$X3) + data$X8 * f3(data$X9) +
     data$X5 * f4(data$X3) + rnorm(n)
   
   data <- data.frame(cbind(Y,data))
   
   
   m0 <- sback(Y ~ sb(X9,by=X1,h=-1) + sb(X3,by=X6,h=-1) + sb(X9,by=X8,h=-1) + sb(X3,by=X5,h=-1), data=data)
   
   pred <- predict(m0, newdata = newdata)
   
   f.hat1 <- (pred$coeff[[grep("X1:X9|X9:X1", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(X9.*by = X1", colnames(pred$peffects))])
   
   val.hat1 <- f.hat1 - mean(f.hat1)
   val.se1  <- (val.hat1 - f1.true)^2
   
   f.hat2 <- (pred$coeff[[grep("X6:X3|X3:X6", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(X3.*by = X6", colnames(pred$peffects))])
   val.hat2 <- f.hat2 - mean(f.hat2)
   val.se2  <- (val.hat2 - f2.true)^2
   
   f.hat3 <- (pred$coeff[[grep("X8:X9|X9:X8", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(X9.*by = X8", colnames(pred$peffects))])
   val.hat3 <- f.hat3 - mean(f.hat3)
   val.se3  <- (val.hat3 - f3.true)^2
   
   # 4. f4 (X5 * f(X3))
   f.hat4 <- (pred$coeff[[grep("X5:X3|X3:X5", names(pred$coeff), value=T)]] * x_grid + 
                pred$peffects[,grep("sb\\(X3.*by = X5", colnames(pred$peffects))])
   val.hat4 <- f.hat4 - mean(f.hat4)
   val.se4  <- (val.hat4 - f4.true)^2
   
   list(
     hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
     se  = list(val.se1, val.se2, val.se3, val.se4)
   )
 }, error = function(e) NULL)
}

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

stopCluster(cl)

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat[[4]],1,mean) - f4.true)^2)

ibs <- c(ibs.1,ibs.2,ibs.3,ibs.4)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se[[4]],1,mean)))

imse <- c(imse.1,imse.2,imse.3,imse.4)

# Integrated MSE plot ----

labels_list <- c("X[1]*f[1](X[9])", "X[6]*f[2](X[3])", "X[8]*f[3](X[9])", 
                 "X[5]*f[4](X[3])")
y_labels <- c("MSE~(f[list(1,9)])", "MSE~(f[list(6,3)])", "MSE~(f[list(8,9)])", 
              "MSE~(f[list(5,3)])")

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

# Colori
line_col <- "#1A237E"
fill_col <- "#7986CB"

# Theme personalizzato
theme_mse <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, 
                              margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

# Genera i plot
plot_list <- list()
for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = mse)) +
    geom_area(fill = fill_col, alpha = 0.3) +
    geom_line(color = line_col, linewidth = 0.7) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels[j])) +
    theme_mse
  
  plot_list[[j]] <- p
}

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise Mean Squared Error (MSE)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 15))
    )
  )

print(combined)

################################################################################

# Preparazione dati
labels_list <- c("X[1]*f[1](X[9])", "X[6]*f[2](X[3])", "X[8]*f[3](X[9])", 
                 "X[5]*f[4](X[3])")
y_labels <- c("Bias~(f[list(1,9)])", "Bias~(f[list(6,3)])", "Bias~(f[list(8,9)])", 
              "Bias~(f[list(5,3)])")

f_true_list <- list(f1.true, f2.true, f3.true, f4.true)

plot_data <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    bias = (apply(f.hat[[j]], 1, mean) - f_true_list[[j]]),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

# Colori
line_col <- "#2E7D32"
fill_col <- "#81C784"

# Theme personalizzato
theme_bias <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, 
                              margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

# Genera i plot
plot_list <- list()
for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = fill_col, alpha = 0.3) +
    geom_line(color = line_col, linewidth = 0.6) +
    coord_cartesian(ylim = c(-0.6, 0.6)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels[j])) +
    theme_bias
  
  plot_list[[j]] <- p
}

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise Bias",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 15))
    )
  )

print(combined)

# Function plot ----
library(ggplot2)
library(patchwork)
library(cowplot)

# Preparazione dati
f_true_funcs <- list(f1, f2, f3, f4)
x_labels <- c("X[9]", "X[3]", "X[9]", "X[3]")
y_labels <- c("X[1]*f[1](X[9])", "X[6]*f[2](X[3])", "X[8]*f[3](X[9])", "X[5]*f[4](X[3])")

true_col <- "#00695C"
est_col <- "#AD1457"
ribbon_col <- "#E0E0E0"
# Theme personalizzato
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

# Genera i plot
plot_list <- list()
for (j in 1:4) {
  
  plot_df <- data.frame(
    x = x_grid,
    true_f = f_true_funcs[[j]](x_grid),
    estimate = apply(f.hat[[j]], 1, mean),
    lower = apply(f.hat[[j]], 1, quantile, 0.025),
    upper = apply(f.hat[[j]], 1, quantile, 0.975)
  )
  
  p <- ggplot(plot_df, aes(x = x)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ribbon_col, alpha = 0.7) +
    geom_line(aes(y = true_f, color = "True function"), linewidth = 1) +
    geom_line(aes(y = estimate, color = "Estimate"), linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
    coord_cartesian(ylim = c(-9, 5)) +
    labs(x = parse(text = x_labels[j]), y = parse(text = y_labels[j])) +
    theme_func
  
  plot_list[[j]] <- p
}
# Legenda separata
legend_plot <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_line(aes(color = "True function"), linewidth = 1) +
  geom_line(aes(color = "Estimate"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = c(0.8, 0.8), linetype = c("solid", "dashed"))))

# Estrai solo la legenda
legend_grob <- cowplot::get_legend(legend_plot)

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "True vs Estimated Functions",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 5))
    )
  )

# Layout finale con legenda
wrap_elements(legend_grob) / combined + plot_layout(heights = c(0.08, 1))



# Save results to JSON ----
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim,
  f1_X9X1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X3X6 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X9X8 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X3X5 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)
write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)


##################
## dep variable ## ----
##################


# SB simulation graphical model

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\Mixed\\500d')

# Gaussian Graphichal model; Graph = "random") ----

n <- 500
p <- 15

set.seed(12)

Gr <- huge.generator(n, p, graph="random", prob=0.23)
Adj.true <- as.matrix(Gr$theta)
dati <- as.data.frame(Gr$data)
qgraph::qgraph(Adj.true, type="UG", edge.color="#7A6F8E", color="#EDE8F2")
colnames(dati) <- paste0("X", 1:p)

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
  X1 = rep(1, n), 
  X6 = rep(1, n), 
  X8 = rep(1, n), 
  X5 = rep(1,n),
  X11 = x_grid,
  X13 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))


sim.results1 <- foreach(i = 1:nsim, .packages = c("wsbackfit", "mvtnorm"), 
                       .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata1", "dpg_sigma","p",
                                   "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
                                     
   set.seed(i)
   
   tryCatch({
     
     data <- data.frame(rmvnorm(n,mean = rep(0,p), sigma = dpg_sigma))
     
     Y <- data$X1 * f1(data$X11) + data$X6 * f2(data$X13) + data$X8 * f3(data$X11) +
       data$X5 * f4(data$X13) + rnorm(n)
     
     data <- data.frame(cbind(Y,data))
     
     m0 <- sback(Y ~ sb(X11,by=X1,h=-1) + sb(X13,by=X6,h=-1) + sb(X11,by=X8,h=-1) + sb(X13,by=X5,h=-1), data=data)
     
     pred <- predict(m0, newdata = newdata1)
     
     f.hat1 <- (pred$coeff[[grep("X1:X11|X11:X1", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X11.*by = X1", colnames(pred$peffects))])
     
     val.hat1 <- f.hat1 - mean(f.hat1)
     val.se1  <- (val.hat1 - f1.true)^2
     
     f.hat2 <- (pred$coeff[[grep("X6:X13|X13:X6", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X13.*by = X6", colnames(pred$peffects))])
     val.hat2 <- f.hat2 - mean(f.hat2)
     val.se2  <- (val.hat2 - f2.true)^2
     
     f.hat3 <- (pred$coeff[[grep("X8:X11|X11:X8", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X11.*by = X8", colnames(pred$peffects))])
     val.hat3 <- f.hat3 - mean(f.hat3)
     val.se3  <- (val.hat3 - f3.true)^2
     
     # 4. f4 (X5 * f(X3))
     f.hat4 <- (pred$coeff[[grep("X5:X13|X13:X5", names(pred$coeff), value=T)]] * x_grid + 
                  pred$peffects[,grep("sb\\(X13.*by = X5", colnames(pred$peffects))])
     val.hat4 <- f.hat4 - mean(f.hat4)
     val.se4  <- (val.hat4 - f4.true)^2
     
     list(
       hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
       se  = list(val.se1, val.se2, val.se3, val.se4)
     )
   }, error = function(e) NULL)
 }

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

stopCluster(cl)

# integrated square bias (ibs)

ibs.1 <- mean((apply(f.hat1[[1]],1,mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat1[[2]],1,mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat1[[3]],1,mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat1[[4]],1,mean) - f4.true)^2)

ibs <- c(ibs.1,ibs.2,ibs.3,ibs.4)

# integrated MSE (imse)

imse.1 <- (6 * mean(apply(f.se1[[1]],1,mean)))
imse.2 <- (6 * mean(apply(f.se1[[2]],1,mean)))
imse.3 <- (6 * mean(apply(f.se1[[3]],1,mean)))
imse.4 <- (6 * mean(apply(f.se1[[4]],1,mean)))

imse <- c(imse.1,imse.2,imse.3,imse.4)

# Integrated MSE plot ----

labels_list <- c("X[1]*f[1](X[11])", "X[6]*f[2](X[13])", "X[8]*f[3](X[11])", 
                 "X[5]*f[4](X[13])")
y_labels <- c("MSE~(f[list(1,11)])", "MSE~(f[list(6,13)])", "MSE~(f[list(8,11)])", 
              "MSE~(f[list(5,13)])")

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

# Colori
line_col <- "#1A237E"
fill_col <- "#7986CB"

# Theme personalizzato
theme_mse <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, 
                              margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

# Genera i plot
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

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise Mean Squared Error (MSE)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 15))
    )
  )

print(combined)

################################################################################

# Preparazione dati
labels_list <- c("X[1]*f[1](X[11])", "X[6]*f[2](X[13])", "X[8]*f[3](X[11])", 
                 "X[5]*f[4](X[13])")
y_labels <- c("Bias~(f[list(1,11)])", "Bias~(f[list(6,13)])", "Bias~(f[list(8,11)])", 
              "Bias~(f[list(5,13)])")

f_true_list <- list(f1.true, f2.true, f3.true, f4.true)

plot_data <- data.frame()
for (j in 1:4) {
  temp_df <- data.frame(
    x = x_grid,
    bias = (apply(f.hat1[[j]], 1, mean) - f_true_list[[j]]),
    Func_Name = labels_list[j]
  )
  plot_data <- rbind(plot_data, temp_df)
}
plot_data$Function <- factor(plot_data$Func_Name, levels = labels_list)

# Colori
line_col <- "#2E7D32"
fill_col <- "#81C784"

# Theme personalizzato
theme_bias <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 11, 
                              margin = margin(b = 10)),
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 8, color = "gray30"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

# Genera i plot
plot_list <- list()
for (j in 1:4) {
  subset_data <- plot_data[plot_data$Func_Name == labels_list[j], ]
  
  p <- ggplot(subset_data, aes(x = x, y = bias)) +
    geom_area(fill = fill_col, alpha = 0.3) +
    geom_line(color = line_col, linewidth = 0.6) +
    coord_cartesian(ylim = c(-0.5, 0.6)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels[j])) +
    theme_bias
  
  plot_list[[j]] <- p
}

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "Pointwise Bias",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 15))
    )
  )

print(combined)

# Function plot ----
library(ggplot2)
library(patchwork)
library(cowplot)

# Preparazione dati
f_true_funcs <- list(f1, f2, f3, f4)
x_labels <- c("X[11]", "X[13]", "X[11]", "X[13]")
y_labels <- c("X[1]*f[1](X[11])", "X[6]*f[2](X[13])", "X[8]*f[3](X[11])", "X[5]*f[4](X[13])")

# Colori
true_col <- "#00695C"
est_col <- "#AD1457"
ribbon_col <- "#E0E0E0"
# Theme personalizzato
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

# Genera i plot
plot_list <- list()
for (j in 1:4) {
  
  plot_df <- data.frame(
    x = x_grid,
    true_f = f_true_funcs[[j]](x_grid),
    estimate = apply(f.hat1[[j]], 1, mean),
    lower = apply(f.hat1[[j]], 1, quantile, 0.025),
    upper = apply(f.hat1[[j]], 1, quantile, 0.975)
  )
  
  p <- ggplot(plot_df, aes(x = x)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ribbon_col, alpha = 0.7) +
    geom_line(aes(y = true_f, color = "True function"), linewidth = 1) +
    geom_line(aes(y = estimate, color = "Estimate"), linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
    coord_cartesian(ylim = c(-10, 5)) +
    labs(x = parse(text = x_labels[j]), y = parse(text = y_labels[j])) +
    theme_func
  
  plot_list[[j]] <- p
}
# Legenda separata
legend_plot <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_line(aes(color = "True function"), linewidth = 1) +
  geom_line(aes(color = "Estimate"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("True function" = true_col, "Estimate" = est_col)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = c(0.8, 0.8), linetype = c("solid", "dashed"))))

# Estrai solo la legenda
legend_grob <- cowplot::get_legend(legend_plot)

# Combina con patchwork
combined <- (plot_list[[1]] | plot_list[[2]]) / 
  (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(
    title = "True vs Estimated Functions",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, 
                                margin = margin(b = 5))
    )
  )

# Layout finale con legenda
wrap_elements(legend_grob) / combined + plot_layout(heights = c(0.08, 1))


# Save results to JSON ----
library(jsonlite)
results_json <- list(
  sample_size = n,
  replications = nsim,
  f1_X11X1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X13X6 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X11X8 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X13X5 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)
write_json(results_json, "results.json", pretty = TRUE, auto_unbox = TRUE)


# aggiungere variabili che non ci sono. Aggiungo una variabile che sta nel grafo ma non nel modello vero
# Permutation test 


