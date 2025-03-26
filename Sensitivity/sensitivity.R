rm(list=ls())

data.df <- read.csv("../Data/NIS-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y

library(deSolve)
library(sensitivity)
library(ggplot2)

Viro.model <- function(t, pop, param) {
  y <- pop[1]
  x <- pop[2]
  v <- pop[3]
  
  K     <- param[1]
  r     <- param[2]
  delta <- param[3]
  rho   <- param[4]
  k     <- param[5]
  a     <- param[6]
  w     <- param[7]
  eps   <- param[8]
  
  dy <- r * y * (1 - ((y+x)^eps / K^eps)) - k * y * v - rho * x * y
  dx <- k * y * v - delta * x
  dv <- a * x - w * v - k * y * v
  
  list(c(dy, dx, dv))
}

MyViro.model <- function(t, pop, param) {
  y <- pop[1]
  x <- pop[2]
  s <- pop[3]
  v <- pop[4]
  
  K     <- param[1]
  r     <- param[2]
  delta <- param[3]
  rho   <- param[4]
  k     <- param[5]
  a     <- param[6]
  w     <- param[7]
  eps   <- param[8]
  gamma <- param[9]
  l     <- param[10]
  
  dy <- r * y * (1 - ((y + x + s)^eps / K^eps)) - rho * x * y - k * y * v
  dx <- k * y * v - delta * x - rho * x * y
  ds <- 2 * rho * x * y - gamma * s
  dv <- a * x - w * v - k * y * v + l * gamma * s
  
  list(c(dy, dx, ds, dv))
}

eval_func_Viro <- function(X) {
  n <- nrow(X)
  outcomes <- numeric(n)
  dt <- 1
  Tmax <- 100
  time <- seq(0, Tmax, by = dt)
  init <- c(y = TumorData[1], x = 0, v = 0)
  for(i in 1:n){
    params <- as.numeric(X[i,])
    event.dat <- data.frame(var   = c("v"),
                            time  = 14,
                            value = 2,
                            method = "replace")
    sim <- lsoda(y = init, times = time, func = Viro.model, parms = params, events = list(data = event.dat))
    sim_df <- as.data.frame(sim)
    final <- tail(sim_df, 1)
    outcomes[i] <- final$y + final$x
  }
  return(outcomes)
}

eval_func_MyViro <- function(X) {
  n <- nrow(X)
  outcomes <- numeric(n)
  dt <- 1
  Tmax <- 100
  time <- seq(0, Tmax, by = dt)
  init <- c(y = TumorData[1], x = 0, s = 0, v = 0)
  for(i in 1:n){
    params <- as.numeric(X[i,])
    event.dat <- data.frame(var   = c("v"),
                            time  = 14,
                            value = 2,
                            method = "replace")
    sim <- lsoda(y = init, times = time, func = MyViro.model, parms = params, events = list(data = event.dat))
    sim_df <- as.data.frame(sim)
    final <- tail(sim_df, 1)
    outcomes[i] <- final$y + final$x + final$s
  }
  return(outcomes)
}

set.seed(123)
n <- 1000

X1_Viro <- data.frame(
  K     = runif(n, 1000, 3000),
  r     = runif(n, 0.1, 1),
  delta = runif(n, 0.05, 0.5),
  rho   = runif(n, 0.0001, 0.01),
  k     = runif(n, 0.0001, 0.01),
  a     = runif(n, 1, 100),
  w     = runif(n, 0.1, 1),
  eps   = runif(n, 0.5, 2)
)

X1_MyViro <- data.frame(
  K     = runif(n, 1000, 3000),
  r     = runif(n, 0.1, 1),
  delta = runif(n, 0.05, 0.5),
  rho   = runif(n, 0.0001, 0.01),
  k     = runif(n, 0.0001, 0.01),
  a     = runif(n, 1, 100),
  w     = runif(n, 0.1, 1),
  eps   = runif(n, 0.5, 2),
  gam = runif(n, 0.01, 0.1),
  l     = runif(n, 0.1, 10)
)

y_Viro <- eval_func_Viro(X1_Viro)
hist(y_Viro)
res_prcc_Viro <- pcc(X1_Viro, y_Viro, nboot = 100, rank = TRUE)
print(res_prcc_Viro)
plot(res_prcc_Viro)
abline(h = 0, col = "grey")

y_MyViro <- eval_func_MyViro(X1_MyViro)
res_prcc_MyViro <- pcc(X1_MyViro, y_MyViro, nboot = 100, rank = TRUE)
print(res_prcc_MyViro)
plot(res_prcc_MyViro)
abline(h = 0, col = "grey")
