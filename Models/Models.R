rm(list=ls())

###################################################################
###################################################################
## Projet
###################################################################
###################################################################

data.df <- read.csv("../Data/NIS-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y



Viro.model <- function(t, pop, param) {
  
  y <- pop[1]  # Uninfected Tumors
  x <- pop[2]  # Infected Tumors
  v <- pop[3]  # Virus
  
  K=param[1] #Carrying capacity of tumors (in 10 ^6 cells)
  r=param[2] #Effective growth rate of uninfected cells (day-1)
  delta=param[3] #Effective death rate constant of infected cells (day-1)
  
  rho=param[4] #Rate constant of cell fusion (per day per 10^6 cells)
  k=param[5] #Infection rate constant (per day per 106 cells or virions) of uninfected cells by free virus v(t)
  a=param[6] #Virus production rate constant by infected cells (virions per day per cell)
  w=param[7] #Rate constant of virus elimination (day-1)
  
  eps=param[8]
  
  #x,y,v the results 
  dy <- r*y*(1-((y+x)^eps/K^eps))-k*y*v-rho*x*y
  dx <- k*y*v - delta*x
  dv <- a*x - w*v - k*y*v
  

  
  res<-c(dy, dx, dv)
  
  list(res)
}


MyViro.model <- function(t, pop, param) {
  y <- pop[1]  # cellules tumorales non infectées
  x <- pop[2]  # cellules tumorales infectées
  s <- pop[3]
  v <- pop[4]  # particules virales
  
  
  K     <- param[1]
  r     <- param[2]
  delta <- param[3]
  rho   <- param[4]
  k     <- param[5]
  a     <- param[6]
  w     <- param[7]
  eps   <- param[8]
  gamma <- param[9]  # death of syncytia
  l    <- param[10]  # release of virus
  
  dy <- r * y * (1 - ((y + x + s)^eps / K^eps)) - rho* x * y -k*y*v
  
  dx <- k*y*v  - delta*x - rho* x * y
  
  ds <- 2*rho*x*y - gamma*s
  
  dv <- a*x - w*v - k*y*v + l*gamma*s
  
  res<-c(dy, dx, ds, dv)
  
  list(res)
}

##############################
### main
##############################

require(deSolve)

Tmax = 45
dt = 1

y0=TumorData[1]
x0=0
s0=0
v0=0

#Estimated param
event.dat <- data.frame(var   = c("v"),
                        time  = 14,
                        value = 2,
                        method = "replace")

Time=seq(from=TimeData[1],to=Tmax,by=dt)

Init.cond <- c(y = y0, x = x0, v = v0)
param=c(2139.258,0.2062134,1.119, 0.141, 0.000591, 0.9, 0.3 ,1.648773)
Init_mymodel.cond <- c(y = y0, x = x0,  s = s0, v = v0)
param_mymodel=c(2139.258, 0.1820776, 1.0143859, 0.2220311, 0.0001,  0.8996148,  0.2717003, 1.648773, 1.0338861, 0.1848987)

# Execute
result <- lsoda(Init.cond, Time, Viro.model, param, events = list(data = event.dat))
result_mymodel <- lsoda(Init_mymodel.cond, Time, MyViro.model, param_mymodel, events = list(data = event.dat))

colnames(result) <- c("Time", "y", "x", "v")
colnames(result_mymodel) <- c("Time", "y", "x", "s", "v")

##############################################################################
# Plot
##############################################################################
plot(TimeData, TumorData, pch=16, col="blue",
     xlab="Time (days)", ylab="Tumor size (mm^3)", bty="n",
     xlim = c(0, Tmax),ylim=c(0, 2139.258))

lines(result[,"Time"], result[,"y"] + result[,"x"], col="green", lwd=2)
lines(result_mymodel[,"Time"], result_mymodel[,"y"] + result_mymodel[,"x"] + result_mymodel[,"s"], col="red", lwd=2, lty=2)

legend("topleft", legend=c("Experimental data", "Bajzer Model", "MyModel"),
       col=c("blue", "green", "red"),
       pch=c(16, NA, NA), lty=c(NA, 1, 2), lwd=c(NA, 2, 2), bty="n")

