library(RAHRS)

# Data
data(magnetometer3k)
data(accelerometer3k)
data(gyroscope3k)
data(accl.coefs)
data(magn.coefs)
gyro.coefs <- vector(mode='numeric',12)

# Number of simulation steps
Nsim <- dim(magnetometer3k)[1]
m <- matrix(unlist(magnetometer3k),ncol=3,nrow=Nsim,byrow=TRUE)
a <- matrix(unlist(accelerometer3k),ncol=3,nrow=Nsim,byrow=TRUE)
w <- matrix(unlist(gyroscope3k),ncol=3,nrow=Nsim,byrow=TRUE)

#Number of simulation steps
Nsim <- dim(a)[1]

#Define Output data Arrays
R_ <- matrix(0,ncol=3,nrow=Nsim)
dw_ <- matrix(0,ncol=3,nrow=Nsim)
w_ <- matrix(0,ncol=3,nrow=Nsim)
a_ <- matrix(0,ncol=3,nrow=Nsim)
m_ <- matrix(0,ncol=3,nrow=Nsim)
TRIAD_ <- matrix(0,ncol=3,nrow=Nsim)

## Calibrate magnetometers
B = matrix(c(magn_coefs[1], magn_coefs[4], magn_coefs[5],
    magn_coefs[6], magn_coefs[2], magn_coefs[7],
    magn_coefs[8], magn_coefs[9], magn_coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(magn_coefs[10],magn_coefs[11],magn_coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) m_[n,] <- t((diag(3)-B) %*% (t(m[n,]) - B0))

## Calibrate accelerometers
B = matrix(c(accl_coefs[1], accl_coefs[4], accl_coefs[5],
    accl_coefs[6], accl_coefs[2], accl_coefs[7],
    accl_coefs[8], accl_coefs[9], accl_coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(accl_coefs[10],accl_coefs[11],accl_coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) a_[n,] <- t((diag(3)-B) %*% (t(a[n,]) - B0))

## Calibrate gyroscopes
B = matrix(c(gyro_coefs[1], gyro_coefs[4], gyro_coefs[5],
    gyro_coefs[6], gyro_coefs[2], gyro_coefs[7],
    gyro_coefs[8], gyro_coefs[9], gyro_coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(gyro_coefs[10],gyro_coefs[11],gyro_coefs[12]),nrow=3,ncol=1)

for (n in 1:Nsim) w_[n,] <- t((diag(3)-B) %*% (t(w[n,]) - B0))

## AHRS Parameters
#Magnetic Field Vector In Navigation Frame
Parameters<-list(mn=0,an=0,dT=0)
Parameters$mn = matrix(c(0.315777529635464, 0.057133095826051, -0.947111588535720),ncol=1)
#Acceleration vector In Navigation Frame
Parameters$an = matrix(c(0, 0, -1),ncol=1)
#Sampling Rate 1/Hz
Parameters$dT  = c(1/100)

#Initial attitude quaternion value
q = EA2Q( 0.0,0.0,0.0,'zyx' ) 
#q = single(q)
#Initial State Vector
Filter<-list(x=0,P=0,R=0,Q=0)
Filter$x = matrix(c(t(q), 0, 0, 0),ncol=1)
#Initial Covariance Matrix
Filter$P = diag(1e-5,7)
#Measuremens Noise
Filter$R = diag(1e0,6)
#Syatem Noise
Filter$Q = diag(c(1e-4, 1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-10))
Sensors<-list(w=0,a=0,m=0)
## Main loop
for (n in 1:Nsim)
    {
    Sensors$w = (matrix(w_[n,],ncol=1))
    Sensors$a = (matrix(a_[n,],ncol=1))
    Sensors$m = (matrix(m_[n,],ncol=1))

    ## AHRS Part
    Filter = ahrs_UKF_QUATERNION(Filter, Sensors, Parameters)
    #[R_(i,1) R_(i,2) R_(i,3)] = ...
  R_[n,] <- Q2EA(Filter$x[1:4])
    dw_[n,] =  c(Filter$x[5], Filter$x[6], Filter$x[7])

    ## TRIAD Part
    W1 = matrix(a_[n,],nrow=1) / norm(matrix(a_[n,],nrow=1),'f') 
    W2 = matrix(m_[n,],nrow=1)/norm(matrix(m_[n,],nrow=1),'f')
    
    V1 = (matrix(Parameters$an,nrow=1))/norm(matrix(Parameters$an,nrow=1),'f')#t
    V2 = (matrix(Parameters$mn,nrow=1))/norm(matrix(Parameters$mn,nrow=1),'f')#t
    
    Ou1 = W1
    Ou2 = cross(W1,W2)/norm(cross(W1,W2),'f')
    Ou3 = cross((W1),cross(W1,W2))/norm(cross(W1,W2),'f')

    R1 = V1
    R2 = (cross(V1,V2)/norm(cross(V1,V2),'f'))#t
    R3 = (cross((V1),cross(V1,V2))/norm(cross(V1,V2),'f'))#t
    
    Mou = cbind(t(Ou1), t(Ou2), t(Ou3))
    Mr = cbind(t(R1), t(R2), t(R3))
    
    A = Mou %*% t(Mr)

    # Calculate angles
   # [TRIAD_(i,1) TRIAD_(i,2) TRIAD_(i,3)] = DCM2EA_eml(A)
    TRIAD_[n,] <- DCM2EA(A)#rotMat2euler(A)
   
}

#pdf('EKF_QUATERNION.pdf')
## Plot Results
#psi - Angle around Z axis
#theta - Angle around Y axis
#gamma - Angle around X axis

plot(TRIAD_[,1],col='red',ylim=c(-4,4), type='l',main='UKF quaternion TRIAD / AHRS',xlab='Time 10 msec (100 Hz)', ylab='Euler angles')
par(new=T)
plot(TRIAD_[,2],col='blue',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(TRIAD_[,3],col='green',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,1],col='orange',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,2],col='cyan',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R_[,3],col='magenta',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
abline(h=(seq(-4,4,1)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(paste(psi,plain(TRIAD))) ,expression(paste(theta,plain(TRIAD))),
expression(paste(gamma,plain(TRIAD))),expression(paste(psi,plain(AHRS))) ,expression(paste(theta,plain(AHRS))),
expression(paste(gamma,plain(AHRS)))),col=c('red','blue','green','orange','cyan','magenta'), lty = c(1, 1, 1, 1, 1, 1),bg='white')



yl<-c(-6*1e-4,4*1e-4)
plot(dw_[,1],col='blue',ylim=yl, type='l',main='UKF quaternion estimated gyroscope bias drift',xlab='Time 10 ms (100 Hz)',ylab='Estimated bias radians/sec')
par(new=T)
plot(dw_[,2],col='green',ylim=yl, type='l',main='',xlab='', ylab='')
par(new=T)
plot(dw_[,3],col='red',ylim=yl, type='l',main='',xlab='', ylab='')
abline(h=(seq(-4*1e-4,4*1e-4,1*1e-4)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(paste(omega,plain(X0))) ,expression(paste(omega,plain(Y0))),
expression(paste(omega,plain(Z0)))),col=c('red','blue','green','orange','cyan','magenta'), lty = c(1, 1, 1))

#dev.off()

