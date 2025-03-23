rm(list=ls())

##############################################################################
# Read CSV data
##############################################################################
data.df <- read.csv("plot-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y
TimeSim <- TimeData

##############################################################################
# Model
##############################################################################
library(deSolve)

Viro.model <- function(t, pop, param) {
  y <- pop[1]  # cellules tumorales non infectées
  x <- pop[2]  # cellules tumorales infectées
  v <- pop[3]  # particules virales
  
  K     <- param[1]
  r     <- param[2]
  delta <- param[3]
  rho   <- param[4]
  k     <- param[5]
  a     <- param[6]
  w     <- param[7]
  eps   <- param[8]
  
  dy <- r*y*(1 - ((y + x)^eps / K^eps)) - k*y*v - rho*x*y
  dx <- k*y*v - delta*x + rho*x*y
  dv <- a*x - w*v - k*y*v
  
  list( c(dy, dx, dv) )
}

event.dat <- data.frame(
  var   = c("y", "v"),   # Quelles variables sont affectées
  time  = 14,            # Moment de l'événement
  value = c(126.237, 2), # Nouvelles valeurs
  method= "replace"      # Remplacer la valeur
)

##############################################################################
# LSE
##############################################################################

# Fixed parameters
K_fixed   <- 2139.258
r_fixed   <- 0.2062134
eps_fixed <- 1.648773

# Initial conditions
Init.cond <- c(y=0, x=0, v=0)

# par = (delta, rho, k, a, w)
sse_func <- function(par) {
  
  delta_ <- par[1]
  rho_   <- par[2]
  k_     <- par[3]
  a_     <- par[4]
  w_     <- par[5]
  
  param_vec <- c(K_fixed, r_fixed, delta_, rho_, k_, a_, w_, eps_fixed)
  
  out <- lsoda(Init.cond, TimeSim, Viro.model, param_vec,
               events = list(data=event.dat))
  
  model_tumor <- out[, 'y'] + out[, 'x']
  
  diff_ <- TumorData - model_tumor
  sse_val <- sum(diff_^2, na.rm=TRUE)
  
  return(sse_val)
}

##############################################################################
# Optimisation des paramètres inconnus (delta, rho, k, a, w)
# On utilise la fonction optim (méthode L-BFGS-B) pour imposer des bornes.
##############################################################################

# possible init values for parameters pour le vecteur par
init_par <- c(delta=1.119, rho=0.141, k=0.000591, a=0.9, w=0.3)

# parameters estimation using 'optim'
fit <- optim(
  par   = init_par,
  fn    = sse_func,
  method= "L-BFGS-B",
  lower = c(1e-20, 1e-20, 1e-20, 1e-20, 1e-20),  # lower bound
  upper = c(100,    10,    10,    100,   100)   # upper bound
)

# Get parameters
best_par <- fit$par
best_par

param_best <- c(K_fixed, r_fixed, best_par["delta"], best_par["rho"],
                best_par["k"], best_par["a"], best_par["w"], eps_fixed)

out_best <- lsoda(Init.cond, TimeSim, Viro.model, param_best,
                  events = list(data = event.dat))
out_best_df <- as.data.frame(out_best)
out_best_df$tumor <- out_best_df$y + out_best_df$x

plot(TimeData, TumorData, pch=16, col="blue",
     xlab="Temps (jour)", ylab="Taille Tumeur (x+y)", bty="n")
lines(out_best_df$time, out_best_df$tumor, col="red", lwd=2)
legend("topright", legend=c("Données", "Modèle"),
       col=c("blue","red"), pch=c(16,NA), lty=c(0,1), bty="n")

