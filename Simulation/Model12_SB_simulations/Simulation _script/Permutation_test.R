# Permutation test

library(wsbackfit)
library(mvtnorm)

perm <- function(x){
  x_new <- sample(x,replace=F)
  return (x_new)
}

# H0: f(z2) = identicamente nulla per ogni z2 nel supporto di Z2

set.seed(666)

n <- 500
nsim <- 100

f1 <- function(x) 2 * exp(-0.5 * x^2)
f2 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)

x1 <- rbeta(n,2,3)
x2 <- rpois(n,5)
z1 <- rnorm(n,0,1)
z2 <- rnorm(n,0,1)
e <- rnorm(n,0,1)

Y <- x1*f1(z1) + x2*f2(z2) + e

df <- data.frame(cbind(Y,x1,x2,z1,z2))

# Model L
model.l <- sback(Y ~ sb(z1,by=x1,h=-1) + sb(z2,by=x2,h=-1),data=df)
mse.l <- mean(model.l$residuals^2)

x_grid <- seq(-3,3,length=n)
newdata <- data.frame(x1=rep(1,n),x2 = rep(1,n), x2.p=rep(1,n),z1=x_grid,z2=x_grid)

pred <- predict.sback(model.l,newdata = newdata)
f1.hat <- pred$coeff['x1'] + pred$coeff['x1:z1'] * x_grid + pred$peffects[,1]
f2.hat <- pred$coeff['x2'] + pred$coeff['x2:z2'] * x_grid + pred$peffects[,2]

par(mfrow=c(1,2))
plot(x_grid,f1.hat,type='l',ylim=c(-2.5,2.5))
lines(x_grid,f1(x_grid),col='red')
plot(x_grid,f2.hat,type='l',ylim=c(-2.5,2.5))
lines(x_grid,f2(x_grid),col='red')
par(mfrow=c(1,1))

f2.sim = matrix(nrow=n, ncol=nsim)

mse.b <- rep(NA,nsim)

for (i in 1:nsim){
  df.b <- df
  df.b$x2.p <- perm(df$x2)
  model.p <- sback(Y ~ sb(z1,by=x1,h=-1) + sb(z2,by=x2.p,h=-1),data=df.b)
  pred <- predict.sback(model.p,newdata = newdata)
  term_2 <- pred$coeff['x2.p']
  term_lin_2 <- pred$coeff[[grep("x2.p:z2|z2:x2.p", names(pred$coeff), value=T)]]
  term_var_2 <- pred$peffects[,2]
  f2.sim[,i] <- term_2 + term_lin_2 * x_grid + term_var_2
  mse.b[i] <- mean(model.p$residuals^2)
}

#p-value
print(sum(mse.b < mse.l)/length(mse.b)) # reject H0

plot(x_grid, x, type='l', ylim=c(-1, 1), 
     main="", ylab="f2(z2)")
abline(h=0, lwd=2, lty=2)

# smooth backfitting permutation test 
# H0: f(z2) = costante, non varia

set.seed(666)

n <- 500
nsim <- 100

f1 <- function(x) 2 * exp(-0.5 * x^2)
f2 <- function(x) 2.5 * sin(2*x) * exp(-0.2*x^2)

x1 <- rbeta(n,2,3)
x2 <- rpois(n,5)
z1 <- rnorm(n,0,1)
z2 <- rnorm(n,0,1)
e <- rnorm(n,0,1)

Y <- x1*f1(z1) + x2*f2(z2) + e

df <- data.frame(cbind(Y,x1,x2,z1,z2))

# Model L
model.l <- sback(Y ~ sb(z1,by=x1,h=-1) + sb(z2,by=x2,h=-1),data=df)
mse.l <- mean(model.l$residuals^2)

mse.b <- rep(NA,nsim)
df.b <- df

for (i in 1:nsim){
  df.b$z2.p <- perm(df$z2)
  model.p <- sback(Y ~ sb(z1,by=x1,h=-1) + sb(z2.p,by=x2,h=-1),data=df.b)
  mse.b[i] <- mean(model.p$residuals^2)
}

#p-value
print(sum(mse.b < mse.l)/length(mse.b)) # reject H0
