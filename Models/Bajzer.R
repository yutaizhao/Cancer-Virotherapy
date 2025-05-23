rm(list=ls())

data.df <- read.csv("../Data/NIS-data.csv", header=TRUE)

TimeData <- data.df$x
TumorData <- data.df$y
###################################################################
###################################################################
## Projet
###################################################################
###################################################################

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

##############################
### main = execution du modele
##############################


require(deSolve)

Tmax = 45
dt = 1

y0=0
x0=0
v0=0
  
#Fitted param
K=  2139.258#Carrying capacity of tumors (in 10 ^6 cells)
r=  0.2062134 #Effective growth rate of uninfected cells (day-1)
eps= 1.648773
#Estimated param
delta= 1.119#Effective death rate constant of infected cells (day-1)
rho= 0.141#Rate constant of cell fusion (per day per 10^6 cells)
k= 0.000591 #Infection rate constant (per day per 106 cells or virions) of uninfected cells by free virus v(t)
a= 0.9 #Virus production rate constant by infected cells (virions per day per cell)
w= 0.3 #Rate constant of virus elimination (day-1)

event.dat <- data.frame(var   = c("y", "v"),
                        time  = 15,
                        value = c(126.237, 2),
                        method = "replace")

Time=seq(from=0,to=Tmax,by=dt)

Init.cond <- c(y = y0, x = x0, v = v0)
param=c(K,r,delta,rho,k,a,w,eps)

# Execute
result <- lsoda(Init.cond, Time, Viro.model, param, events = list(data = event.dat))
result

colnames(result) <- c("Time", "y", "x", "v")

head(result)
plot(TimeData, TumorData, pch=16, col="blue",
     xlab="Time (days)", ylab="Tumor size (x + y)", bty="n",
     xlim = c(0, Tmax),ylim=c(0, 2139.258))
lines(Time,result[,"y"] + result[,"x"],type="l",col="green",xlab="Time",ylab="",ylim=c(0,K),bty="n")
legend("topright",c("data","tumor"),col=c("blue","green"),lty=1,bty="n")

