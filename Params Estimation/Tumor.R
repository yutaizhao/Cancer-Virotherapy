rm(list=ls())

# Read CSV data
data.df <- read.csv("../Data/tumor-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y
TimeSim <- TimeData

# Tumor growth model
# y' = (g / epsilon) * y * [1 - (y / K)^epsilon]
library(deSolve)

TumorGrowth.model <- function(t, pop, param) {
  y <- pop[1]
  
  g     <- param[1]
  eps   <- param[2]
  K     <- param[3]
  
  dy <- (g / eps) * y * (1 - (y / K)^eps)
  
  list( dy )
}

sse_func <- function(par) {
  
  g_   <- par[1]
  eps_ <- par[2]
  K_   <- par[3]
  
  param_vec <- c(g_, eps_, K_)
  
  # Initial condition: y(0) = the first observed tumor size
  Init.cond <- c(y = TumorData[1])
  
  out <- lsoda(Init.cond, TimeSim, TumorGrowth.model, param_vec)
  
  model_y <- out[, 'y']
  
  diff_ <- TumorData - model_y
  
  return(sum(diff_^2))   
}


# Initial guesses for the parameters to be estimated
init_par <- c(g = 0.5, eps = 1, K = 2000)

fit <- optim(
  par   = init_par,
  fn    = sse_func,
  method= "L-BFGS-B",
  lower = c(1e-5, 1e-5, 1),
  upper =  c(10, 10, 10000)
)

# Results
best_par <- fit$par
best_par

param_best <- c(best_par["g"], best_par["eps"], best_par["K"])

Init.cond <- c(y = TumorData[1])

out_best <- lsoda(Init.cond, TimeSim, TumorGrowth.model, param_best)

out_best_df <- as.data.frame(out_best)

plot(TimeData, TumorData, pch=16, col="blue",
     xlab="Time (days)", ylab="Tumor Size", bty="n")
lines(out_best_df$time, out_best_df$y, col="red", lwd=2)

legend("topleft", legend=c("Observed data", "Fit BR model"),
       col=c("blue","red"), pch=c(16,NA), lty=c(0,1), bty="n")

