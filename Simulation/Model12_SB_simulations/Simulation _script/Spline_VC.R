library(foreach)
library(doParallel)
library(huge)
library(mgcv)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)

setwd('C:\\Users\\gargh\\Documents\\Tesi\\Code\\Simulation\\Model12_SB_simulations\\P-spline\\500_dependence')

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

# Coefficient functions ----
f1 <- function(x) x                          
f2 <- function(x) 2 * tanh(x)                
f3 <- function(x) 2 * exp(-0.5 * x^2)
f4 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)

# Setup parallel ----
n_cores <- parallel::detectCores(logical = FALSE) - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# MC simulation parameters ----
nsim <- 100
x_grid <- seq(-3, 3, length = n)

newdata <- data.frame(
  X2 = rep(1, n), 
  X6 = rep(1, n), 
  X8 = rep(1, n), 
  X5 = rep(1, n),
  X1 = x_grid,
  X3 = x_grid
)

f1.true <- f1(x_grid) - mean(f1(x_grid))
f2.true <- f2(x_grid) - mean(f2(x_grid))
f3.true <- f3(x_grid) - mean(f3(x_grid))
f4.true <- f4(x_grid) - mean(f4(x_grid))

# MC simulation ----
sim.results <- foreach(i = 1:nsim, .packages = c("mgcv", "mvtnorm"),
                       .export = c("f1", "f2", "f3", "f4", "n", "x_grid", "newdata", 
                                   "dpg_sigma", "p", "f1.true", "f2.true", "f3.true", "f4.true")) %dopar% {
  
  set.seed(i)
  
  tryCatch({
    
    data <- data.frame(rmvnorm(n, mean = rep(0, p), sigma = dpg_sigma))
    colnames(data) <- paste0("X", 1:p)
    
    Y <- data$X2 * f1(data$X1) + data$X6 * f2(data$X3) + data$X8 * f3(data$X1) +
      data$X5 * f4(data$X3) + rnorm(n)
    
    data$Y <- Y
    
    fit <- gam(Y ~ s(X1, bs='ps', by=X2) + s(X3, bs='ps', by=X6) + 
                 s(X1, bs='ps', by=X8) + s(X3, bs='ps', by=X5), 
               data = data, method = "REML")
    
    pred <- predict(fit, newdata = newdata, type = "terms")
    
    # f1: s(X9):X1
    f.hat1 <- pred[, grep("s\\(X1\\):X2", colnames(pred))]
    val.hat1 <- f.hat1 - mean(f.hat1)
    val.se1 <- (val.hat1 - f1.true)^2
    
    # f2: s(X3):X6
    f.hat2 <- pred[, grep("s\\(X3\\):X6", colnames(pred))]
    val.hat2 <- f.hat2 - mean(f.hat2)
    val.se2 <- (val.hat2 - f2.true)^2
    
    # f3: s(X9):X8
    f.hat3 <- pred[, grep("s\\(X1\\):X8", colnames(pred))]
    val.hat3 <- f.hat3 - mean(f.hat3)
    val.se3 <- (val.hat3 - f3.true)^2
    
    # f4: s(X3):X5
    f.hat4 <- pred[, grep("s\\(X3\\):X5", colnames(pred))]
    val.hat4 <- f.hat4 - mean(f.hat4)
    val.se4 <- (val.hat4 - f4.true)^2
    
    list(
      hat = list(val.hat1, val.hat2, val.hat3, val.hat4),
      se = list(val.se1, val.se2, val.se3, val.se4),
      x_data = list(data$X1, data$X3, data$X1, data$X3)
    )
  }, error = function(e) NULL)
}

stopCluster(cl)

# Filter valid results ----
valid_results <- sim.results[sapply(sim.results, function(x) {
  !is.null(x) && length(x$hat[[1]]) == n
})]
nsim_valid <- length(valid_results)
cat("Simulazioni valide:", nsim_valid, "su", nsim, "\n")

# Populate matrices ----
f.hat <- vector("list", 4)
f.se <- vector("list", 4)
for(j in 1:4) {
  f.hat[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
  f.se[[j]] <- matrix(NA, nrow = n, ncol = nsim_valid)
}

for(i in 1:nsim_valid) {
  for(j in 1:4) {
    f.hat[[j]][, i] <- valid_results[[i]]$hat[[j]]
    f.se[[j]][, i] <- valid_results[[i]]$se[[j]]
  }
}

# Integrated squared bias (IBS) ----
ibs.1 <- mean((apply(f.hat[[1]], 1, mean) - f1.true)^2)
ibs.2 <- mean((apply(f.hat[[2]], 1, mean) - f2.true)^2)
ibs.3 <- mean((apply(f.hat[[3]], 1, mean) - f3.true)^2)
ibs.4 <- mean((apply(f.hat[[4]], 1, mean) - f4.true)^2)

ibs <- round(c(ibs.1, ibs.2, ibs.3, ibs.4),3)

# Integrated MSE (IMSE) ----
imse.1 <- 6 * mean(apply(f.se[[1]], 1, mean))
imse.2 <- 6 * mean(apply(f.se[[2]], 1, mean))
imse.3 <- 6 * mean(apply(f.se[[3]], 1, mean))
imse.4 <- 6 * mean(apply(f.se[[4]], 1, mean))

imse <- round(c(imse.1, imse.2, imse.3, imse.4),3)

# MSE Plot ----
labels_list <- c("X[2]*f[1](X[1])", "X[6]*f[2](X[3])", "X[8]*f[3](X[1])", 
                 "X[5]*f[4](X[3])")
y_labels <- c("MSE~(f[list(1,1)])", "MSE~(f[list(6,3)])", "MSE~(f[list(8,1)])", 
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
    coord_cartesian(ylim = c(0, 1.7)) +
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

# Bias Plot ----
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
    coord_cartesian(ylim = c(-0.3, 0.3)) +
    labs(title = parse(text = labels_list[j]), y = parse(text = y_labels_bias[j])) +
    theme_bias
  
  plot_list_bias[[j]] <- p
}

combined_bias <- (plot_list_bias[[1]] | plot_list_bias[[2]]) / 
  (plot_list_bias[[3]] | plot_list_bias[[4]]) +
  plot_annotation(
    title = "Pointwise Bias (GAM P-splines)",
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
    coord_cartesian(ylim = c(-4, 4)) +
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
  replications = nsim_valid,
  method = "GAM P-splines (REML)",
  f1_X2X1 = list(ibs = unname(ibs[1]), imse = unname(imse[1])),
  f2_X6X3 = list(ibs = unname(ibs[2]), imse = unname(imse[2])),
  f3_X8X1 = list(ibs = unname(ibs[3]), imse = unname(imse[3])),
  f4_X5X3 = list(ibs = unname(ibs[4]), imse = unname(imse[4]))
)
write_json(results_json, "results_gam.json", pretty = TRUE, auto_unbox = TRUE)