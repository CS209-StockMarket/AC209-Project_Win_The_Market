################
# Neil's generic GARCH code
################


##################
# simulate GARCH model, mem = alpha+beta, delta = beta/(alpha+beta)
##################

sim_GARCH <- function(n,alpha,mem,delta) {
  beta = delta*mem
  gamma = (1.0-delta)*mem
  s2 = alpha/(1.0-beta-gamma)
  y = array(0,dim=c(n,2))
  y[1,1] = sqrt(s2)*rnorm(1,0,1)
  y[1,2] = sqrt(s2)
  
  for (i in 2:n) {
    s2 = alpha + (beta*(y[i-1]^2)) + (gamma*s2)
    y[i,1] = sqrt(s2)*rnorm(1,0,1)
    y[i,2] = sqrt(s2)
  }
  y
}

###################
# sd from GARCH recursion
###################

sdGARCH <- function(alpha,beta,gamma) {
  
  n = length(y)
  vol = array(0,dim=c(n,1))
  s2 = alpha/(1.0-beta-gamma); 
  
  for (i in 1:n) {
    vol[i] = sqrt(s2)
    s2 = alpha + (beta*(y[i]^2)) + (gamma*s2)    
  }  
  vol
  
}

#####################
# logL from QGARCH 
######################

LL_GARCH <- function(x) { # minus logL
  alpha=x[1];   mem=x[2];   delta=x[3]
  
  if (length(x)==4) {
    zeta = x[4]
  }
  else { zeta = 0.0 }
  
  beta = delta*mem;   gamma = (1.0-delta)*mem
  s2 = alpha/(1.0-beta-gamma); logL=0.0
  n = length(y)
  
  for (i in 1:n) {
    #    logL = logL + log(dnorm(y[i],0,sqrt(s2)))
    logL = logL -(0.5*log(s2)) - (0.5*((y[i]^2)/s2))
    s2 = alpha + (beta*((y[i]-zeta)^2)) + (gamma*s2)
    
  }
  print(cbind(mem,delta,alpha,beta,gamma,zeta,logL))
  -logL  
}


Z <- read.csv("SP500.csv")[,2] # S&P downloaded from St Louis Fed, no dividends 
Z= rev(Z)
T = length(Z)-1
x<- 1956 + (1.0/252.0)*seq(1,T) 
x1<- 1956 + (1.0/252.0)*seq(1,(T+1)) 

pdf("SP500summary.pdf")
par(mfcol=c(2,2), mar=c(4,2,3,0), oma=c(1.5,2,1,1)) # make R plot 2 plots side by side

plot(x1,Z,type="l",main="S&P level: no dividends")
plot(x1,log(Z),type="l",main="S&P log level")
r = 100*diff(log(Z))
plot(x,r,type="l",main="S&P geometric returns")
r[is.na(r)] <- 0.0 # set any missing returns to zero

xHist = hist(r,breaks=50,plot=FALSE)
plot(xHist$mid,log(xHist$density),main="Log density of returns",type="l")
dev.off()

pdf("SP500dynamics.pdf")
par(mfcol=c(2,2), mar=c(4,2,3,0), oma=c(1.5,2,1,1)) # make R plot 2 plots side by side

xAcf=acf(r,plot=FALSE)
plot(xAcf[2:(length(r)-1)],xlab="lags",main="ACF of returns")

xAcf=acf(abs(r),plot=FALSE)
plot(xAcf[2:(length(r)-1)],xlab="lags",main="ACF of absolute returns")

xAcf=acf((r^2),plot=FALSE)
plot(xAcf[2:(length(r)-1)],xlab="lags",main="ACF of squared returns")
dev.off()



y = r
print(y[1:25])

aRes = optim(c(0.05,0.96,0.2), LL_GARCH,method = "L-BFGS-B",lower=c(0.00001,0.001,0.0001),upper=c(1.0, 0.9999,0.9999)) # minimisation -LogL
aRes$par
aRes$value
aAlpha= (aRes$par)[1]
aBeta = (aRes$par)[2]*(aRes$par)[3]
aGamma= (aRes$par)[2]-((aRes$par)[2]*(aRes$par)[3])

aVol = sdGARCH(aAlpha,aBeta,aGamma)

pdf("SP500garchFit.pdf")
par(mfcol=c(2,3), mar=c(4,2,3,0), oma=c(1.5,2,1,1)) # make R plot 2 plots side by side

plot(x,y/sd(y),type="p",pch=20,main="Geometric returns",xlab="Time")

plot(x,aVol,type="l",main="Cond sd",xlab="Time")
plot(x,(y/aVol),type="p",pch=20,main="Innovations",xlab="Time")

xHist = hist((y/aVol),breaks=50,plot=FALSE)

plot(xHist$mid,log(xHist$density),main="Log density",type="l")

xAcf=acf((y/aVol),plot=FALSE)
plot(xAcf[2:(length(r)-1)],xlab="lags",main="ACF of innovations")

xAcf=acf(abs(y/aVol),plot=FALSE)
plot(xAcf[2:(length(r)-1)],xlab="lags",main="ACF of abs innovations")

dev.off()

install.packages("MASS")
library(MASS)
fitdistr(y, "t")
fitdistr(y/aVol, "t")

aRes = optim(c(0.05,0.96,0.2), LL_GARCH,method = "L-BFGS-B",lower=c(0.00001,0.001,0.0001),upper=c(1.0, 0.9999,0.9999)) # minimisation -LogL

aRes = optim(c(0.05,0.96,0.2,0.03), LL_GARCH,method = "L-BFGS-B",lower=c(0.00001,0.0,0.0,-0.2),upper=c(1.0, 0.9999,0.9999,0.2)) # minimisation -LogL

