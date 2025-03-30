rm(list=ls())
library(deSolve)

##############################################################################
# Read CSV data
##############################################################################
data.df <- read.csv("../Data/NIS-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y
TimeSim <- TimeData

K_fixed   <- 2139.258  # carrying capacity (fixed)
eps_fixed <- 1.648773

##############################################################################
# Virotherapy Model with r and eps as parameters to estimate
##############################################################################
Viro.model <- function(t, pop, param) {
  # state variables
  y <- pop[1]  # uninfected
  x <- pop[2]  # infected
  v <- pop[3]  # virus
  s <- pop[4]  # syncytia

  # fixed param
  K <- K_fixed
  eps   <- eps_fixed  # logistic exponent
  

  # parameters to estimate
  r     <- param[1]  # growth rate
  w     <- param[2]  # virus elimination
  rho   <- param[3]  # infected+uninfected => syncytia
  l     <- param[4]  # release of virus
  k_inf <- param[5]  # infection consumption (y, v)
  delta <- param[6]  # infected death
  gamma <- param[7]  # syncytia death
  a     <- param[8]  # virus production
 

  # uninfected logistic
  dy <- r * y * (1 - ((y + x + s)^eps / K^eps)) - rho * x * y - k_inf * y * v
  dx <- k_inf * y * v - delta * x - rho * x * y
  ds <- 2*rho * x * y - gamma * s
  dv <- a * x - w * v - k_inf * y * v + l * gamma * s

  list(c(dy, dx, dv, ds))
}

##############################################################################
# Event: Virus injection at day = 14
##############################################################################
event.dat <- data.frame(
  var   = "v",
  time  = 14,
  value = 2,
  method= "replace"
)

##############################################################################
# Initial conditions
##############################################################################
Init.cond <- c(y=TumorData[1], x=0, v=0, s=0)

##############################################################################
# SSE function
##############################################################################
sse_func <- function(par) {
  param_vec <- par

  out <- lsoda(Init.cond, TimeSim, Viro.model, param_vec, rtol = 1e-8, atol = 1e-10, events=list(data=event.dat))

  model_tumor <- out[,"y"] + out[,"x"] + out[,"s"]
  model_tumor_data <- model_tumor[TimeSim %in% TimeData]

  SSE <- sum((TumorData - model_tumor_data)^2, na.rm=TRUE)
  return(SSE)
}

##############################################################################
# Initial guess and bounds
##############################################################################
init_par <- c(
  r=0.2,
  w=0.3,
  rho=0.141,
  l=0.1,
  k_inf=0.000591,
  delta=1.119,
  gamma=0.1,
  a=0.9
)

lower_bounds <- c(
  0.01,  0, 0.01, 0, 0.000001, 0.01, 0.01, 0
)

upper_bounds <- c(
  1,   1,  1, 1, 0.0001,    2,     10,    10
)

##############################################################################
# Optimization
##############################################################################
fit <- optim(
  par   = init_par,
  fn    = sse_func,
  method= "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds
)

best_par <- fit$par
print(best_par)

##############################################################################
# Simulate with best-fit parameters
##############################################################################
Time <- seq(from=TimeData[1], to=360, by=1)
out_best <- lsoda(Init.cond, Time, Viro.model, best_par, rtol = 1e-8, atol = 1e-10, events=list(data=event.dat))
out_best_df <- as.data.frame(out_best)
out_best_df$tumor <- out_best_df$y + out_best_df$x + out_best_df$s

##############################################################################
# Plot
##############################################################################
plot(TimeData, TumorData, pch=16, col="blue",
     xlab="Time (days)", ylab="Tumor size (y + x + s)", bty="n")
lines(out_best_df$time, out_best_df$tumor, col="red", lwd=2)

legend("topright", legend=c("Data","Tumor"),
       col=c("blue","red"),
       pch=c(16, NA), lty=c(NA, 1), lwd=c(NA, 2), bty="n")

