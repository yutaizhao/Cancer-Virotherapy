
rm(list=ls())
library(deSolve)
##############################################################################
# Read CSV data
##############################################################################
data.df <- read.csv("../Data/NIS-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y
TimeSim <- TimeData


K_fixed   <- 2139.258  # carrying capacity
r_fixed   <- 0.2062134 # growth rate day^-1
eps_fixed <- 1.648773  # logistic exponent

Viro.model <- function(t, pop, param) {
  # state variables
  y <- pop[1]  # uninfected
  x <- pop[2]  # infected
  v <- pop[3]  # virus
  s <- pop[4]  # syncytia
  
  # fixed param
   K     <- K_fixed
   r     <- r_fixed
   eps   <- eps_fixed

   # free param
   rho  <- param[1]  # infected+uninfected => syncytia
   l    <- param[2]  # release of virus
   k_inf <- param[3]  # infection consumption (y, v)
   delta <- param[4]  # infected death
   gamma <- param[5]  # syncytia death
   a     <- param[6]  # virus production
   w     <- param[7]  # virus elimination
  
  # uninfected logistic
  dy <- r * y * (1 - ((y + x + s)^eps / K^eps)) - rho* x * y -k_inf*y*v
  
  dx <- k_inf*y*v  - delta*x - rho* x * y
  
  ds <- rho*x*y - gamma*s
  
  dv <- a*x - w*v - k_inf*y*v + l*gamma*s
  
  list(c(dy, dx, dv, ds))
}

##############################################################################
# Example event: Virus injection at day = 10
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
  
  out <- lsoda(Init.cond, TimeSim, Viro.model, param_vec, events=list(data=event.dat))
  
  
  # total tumor = y + x + s
  model_tumor <- out[,"y"] + out[,"x"] + out[,"s"]
  
  # match TimeSim to TimeData
  model_tumor_data <- model_tumor[TimeSim %in% TimeData]
  
  diff_ <- TumorData - model_tumor_data
  SSE <- sum(diff_^2, na.rm=TRUE)
  
  return(SSE)
}

##############################################################################
# Parameter bounds and initial guess
##############################################################################
init_par <- c(
  rho1=0.141,
  l=0.1,
  k_inf=0.000591,
  delta=1.119,
  gamma=0.1,
  a=0.9,
  w=0.3
)

lower_bounds <- c(
  1e-10,0,1e-10,1e-10,1e-10,1e-10,1e-10
)
upper_bounds <- c(
  1, 1e6,   1,   10,   10, 1e6,   10
)

##############################################################################
# Use optim (L-BFGS-B) for local minimization
##############################################################################
fit <- optim(
  par   = init_par,
  fn    = sse_func,
  method= "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds
)

best_par <- fit$par
best_par



##############################################################################
# Simulate with best par
##############################################################################
Time=seq(from=TimeData[1],to=360,by=1)
out_best <- lsoda(Init.cond, Time, Viro.model, best_par, events=list(data=event.dat))
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
       pch=c(16,NA,NA), lty=c(NA,1,2), lwd=c(NA,2,2), bty="n")
  
