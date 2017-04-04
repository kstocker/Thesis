#install.packages("deSolve")
library(deSolve)
SIR <- function(t, y, parms){
  ds <- -1.66*y[1]*y[2] 
  di <- 1.66*y[1]*y[2] - 0.455*y[2]
  dr <- 0.455*y[2]
  list(c(ds, di, dr))
}
yini <- c(y1=0.996, y2=0.004, y3=0)
times <- seq(from=0, to=15, by=.1)
out <- ode(times=times, y=yini, func=SIR, parms=NULL, method = "euler")

SIRdata <- head(out, n=16)
SIRdata <- data.frame(SIRdata)


##Backwards Problem###
SIRtest <- function(t, b, y, parms){
  ds <- -b*y[1]*y[2] 
  di <- b*y[1]*y[2] - 0.455*y[2]
  dr <- 0.455*y[2]
  list(c(ds,di,dr))
}
b <- 1.66
yini <- c(y1=0.996, y2=0.004, y3=0)
times <- seq(from=0, to=15, by=1)
#parms <- (b = 1.66)
out <- ode(times=times, y=yini, func=SIRtest, parms=NULL)
out




# #Without Desolve:
# func <- (function(t,x,y) {1.66*x*y - 0.455*y})
# e1 <- euler(func, h = 0.1, t0=0, x0 = 0.996, y0=0.004, tfinal=15)
# e1
# #e1 <- euler(x0=0.004, y0-0.996, times=c(0:15),)
# help("euler")
# 
# func <- function(x,y){x-y}
# ex <- euler(func, h=0.2, x0=0, y0=1, xfinal=1)
# ex

# try1 <- function(t, x, parms) {
#   with(as.list(parms), {
#     ds <- - 1.66*x[1]*x[2]
#     di <- 1.66*x[1]*x[2] - 0.455*x[2]
#     dr <- 0.455*x[2]
#     list(ds, di, dr)
#   })
# }
# 
# time  <- 0:100
# N0    <- 0.1; x<- c(0.996, 0.004, 0)
# x <- c(N = N0)
# 
# time<-seq(0, 100, 2)
# out <- as.data.frame(euler(x, time, try1))