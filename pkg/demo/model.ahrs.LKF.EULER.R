library(RAHRS)

# Data
data(anglesGyroLib)
#data(accl.coefs)
#data(magn.coefs)

accl.coefs <- matrix(c(-0.0771, -0.0817, -0.0769,  0.0066, -0.0032,  0.0001,  0.0090, -0.0023, -0.0035, -0.0281,  0.0430,  0.0306),ncol=1)
magn.coefs <- matrix(c(-0.3801, -0.3108, -0.3351,  0.0342,  0.0109, -0.0066,  0.0070, -0.0696, -0.0492,  0.0552,  0.0565,  0.0190),ncol=1)
gyro.coefs <- vector(mode='numeric',12)

# Number of simulation steps
Nsim <- dim(anglesGyroLib)[1]
m <- matrix(unlist(anglesGyroLib[,c('Mx','My','Mz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000
a <- matrix(unlist(anglesGyroLib[,c('Ax','Ay','Az')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000
w <- matrix(unlist(anglesGyroLib[,c('Wx','Wy','Wz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000

#Data length
Nsim <- dim(a)[1]

#Preallocate Data Logs
Psi_ <- matrix(0,ncol=1,nrow=Nsim)
Theta_ <- matrix(0,ncol=1,nrow=Nsim)
Gamma_ <- matrix(0,ncol=1,nrow=Nsim)
PsiMgn_ <- matrix(0,ncol=1,nrow=Nsim)
ThetaAcc_ <- matrix(0,ncol=1,nrow=Nsim)
GammaAcc_ <- matrix(0,ncol=1,nrow=Nsim)
w_ <- matrix(0,ncol=3,nrow=Nsim)
a_ <- matrix(0,ncol=3,nrow=Nsim)
m_ <- matrix(0,ncol=3,nrow=Nsim)
dw_ <- matrix(0,ncol=3,nrow=Nsim)
dB_ <- matrix(0,ncol=3,nrow=Nsim)
dG_ <- matrix(0,ncol=6,nrow=Nsim)


## Calibrate magnetometers
B = matrix(c(magn.coefs[1], magn.coefs[4], magn.coefs[5],
    magn.coefs[6], magn.coefs[2], magn.coefs[7],
    magn.coefs[8], magn.coefs[9], magn.coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(magn.coefs[10],magn.coefs[11],magn.coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) m_[n,] <- t((diag(3)-B) %*% (matrix(m[n,],3) - B0))

## Calibrate accelerometers
B = matrix(c(accl.coefs[1], accl.coefs[4], accl.coefs[5],
    accl.coefs[6], accl.coefs[2], accl.coefs[7],
    accl.coefs[8], accl.coefs[9], accl.coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(accl.coefs[10],accl.coefs[11],accl.coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) a_[n,] <- t((diag(3)-B) %*% (matrix(a[n,],3) - B0))

## Calibrate gyroscopes
B = matrix(c(gyro.coefs[1], gyro.coefs[4], gyro.coefs[5],
    gyro.coefs[6], gyro.coefs[2], gyro.coefs[7],
    gyro.coefs[8], gyro.coefs[9], gyro.coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(gyro.coefs[10],gyro.coefs[11],gyro.coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) w_[n,] <- t((diag(3)-B) %*% (matrix(w[n,],3) - B0))

Parameters<-list(mn=0,an=0,dT=0,declination=0)
Parameters$declination <- 10.2*pi/180
## Initial quaternion
ax_init <- mean(a_[1:1000,1])
ay_init <- mean(a_[1:1000,2])
az_init <- mean(a_[1:1000,3])
mx_init <- mean(m_[1:1000,1])
my_init <- mean(m_[1:1000,2])
mz_init <- mean(m_[1:1000,3])
Mb_init <- matrix(c(mx_init, my_init, mz_init),ncol=1)
Theta_init <-  atan2(ax_init, sqrt(ay_init^2+az_init^2))
Gamma_init <- -atan2(ay_init, sqrt(ax_init^2+az_init^2))
#%Horizontal projection of magnetic field
State <- list(q=0,P=0,dw=0,dB=0,dG=0)
CbnH_init <- EA2DCM(c(0, Theta_init, Gamma_init))
Mh_init <- t(CbnH_init) %*% Mb_init
Psi_init <- -atan2(Mh_init[2],Mh_init[1])+Parameters$declination
State$q <- EA2Q(Psi_init, Theta_init, Gamma_init,'zyx')
## Kalman Filter Internals
Parameters$R <- diag(.1,3)
Parameters$R[3,3] <- 10
Parameters$Q <- diag(c(1e-5,1e-5,1e-5,1e-8,1e-8,1e-8,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10))
Parameters$dT <- 1/100
Parameters$declination <- 10.2*pi/180
State$P <- diag(rep(1e-5,15))
State$dw <- matrix(0,ncol=1,nrow=3)
State$dB <- matrix(0,ncol=1,nrow=3)
State$dG <- matrix(0,ncol=1,nrow=6)
Sensors<-list(w=0,a=0,m=0)

## Main loop
for (n in 1:Nsim)
    {
#n<-1
#Sensors readings in Body frame
    Sensors$w <- (matrix(w_[n,],ncol=1))
    Sensors$a <- (matrix(a_[n,],ncol=1))
    Sensors$m <- (matrix(m_[n,],ncol=1))
    tmp <- ahrs_LKF_EULER(Sensors, State, Parameters)
Attitude <- tmp$Attitude
State <- tmp$State
    ## write logs
    Psi_[n,1] <- Attitude[1]
    Theta_[n,1] <- Attitude[2]
    Gamma_[n,1] <- Attitude[3]    
    PsiMgn_[n,1] <- Attitude[4]
    ThetaAcc_[n,1] <- Attitude[5]
    GammaAcc_[n,1] <- Attitude[6] 
    dw_[n,] <- State$dw 
    dB_[n,] <- State$dB
    dG_[n,] <- State$dG
    #list(w_=w_, a_=a_, m_=m_, R_=R_, dw_=dw_, dB_=dB_, dG_=dG_, TRIAD_=TRIAD_)
}
R_ <- cbind(Psi_, Theta_, Gamma_)

#R_$Psi_ <- Psi_ R_$Theta_ <- Theta_ R_$Gamma_ <- Gamma_
TRIAD_ <- cbind(PsiMgn_, ThetaAcc_, GammaAcc_)
str(TRIAD_)
#TRIAD_$PsiMgn_ <- PsiMgn_ TRIAD_$ThetaAcc_ <- ThetaAcc_ TRIAD_$GammaAcc_ <- GammaAcc_

#postscript('LKF_EULER.pdf')
## Plot Results
#psi - Angle around Z axis
#theta - Angle around Y axis
#gamma - Angle around X axis

yl<-c(-200,200)
plot(TRIAD_[,1]*180/pi,col='blue',ylim=yl, type='l',main='LKF EULER TRIAD / AHRS',xlab='Time 10 msec (100 Hz)', ylab='Euler angles')
par(new=T)
plot(TRIAD_[,2]*180/pi,col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(TRIAD_[,3]*180/pi,col='red',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,1]*180/pi,col='cyan',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,2]*180/pi,col='magenta',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,3]*180/pi,col='orange',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=(seq(-200,200,100)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(paste(psi,plain(TRIAD))) ,expression(paste(theta,plain(TRIAD))),
expression(paste(gamma,plain(TRIAD))),expression(paste(psi,plain(AHRS))) ,expression(paste(theta,plain(AHRS))),
expression(paste(gamma,plain(AHRS)))),col=c('blue','green','red','cyan','magenta','orange'), lty = c(1, 1, 1, 1, 1, 1),bg='white')

par(mfrow = c(3,1))
yl<- c(-0.01,0.02)
plot(dw_[,1],col='blue',ylim=yl, type='l',main='LKF EULER estimated gyroscope bias drift',xlab='Time 10 ms (100 Hz)',ylab='Estimated bias radians/sec')
par(new=T)
plot(dw_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dw_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=(seq(-0.01,0.02,0.01)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(omega[X0]) ,expression(omega[Y0]),
expression(omega[Z0])),col=c('blue','green','red'), lty = c(1, 1, 1),bg='white')

yl<- c(-2*1e-4,4*1e-4)
plot(dB_[,1],col='blue',ylim=yl, type='l',main='LKF EULER estimated gyroscope bias drift',xlab='Time 10 ms (100 Hz)',ylab='Estimated bias radians/sec')
par(new=T)
plot(dB_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dB_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=(seq(-2*1e-4,4*1e-4,1*1e-4)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(B[xx]) ,expression(B[yy]),
expression(B[zz])),col=c('blue','green','red'), lty = c(1, 1, 1),bg='white')

yl<- c(-2*1e-3,1e-3)
plot(dG_[,1],col='blue',ylim=yl, type='l',main='LKF EULER estimated gyroscope bias drift',xlab='Time 10 ms (100 Hz)',ylab='Estimated bias radians/sec')
par(new=T)
plot(dG_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dG_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dG_[,4],col='cyan',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dG_[,5],col='magenta',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dG_[,6],col='orange',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=(seq(-2*1e-3,1e-3,1e-3)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(G[xy]),expression(G[xz]),expression(G[yx]),expression(G[yz]),expression(G[zx]),expression(G[zy]) ),
col=c('blue','green','red','cyan','magenta','orange'), lty = c(1, 1, 1, 1, 1, 1),bg='white')


#dev.off()




