###############################
#This is the script to reproduce Table 3 from the paper
#Asymptotic distribution of one-component partial least
#squares regression estimators in high dimension
#by Jeronimo Basa, Dennis Cook, Liliana Forzani and Miguel Marcos
###############################






source("howmany.R")
library(pls)
library(MASS)
library(psych)
X <- read.table("Mcal1", quote="\"", comment.char="")
y <- read.table("ycal1", quote="\"", comment.char="")
X <- as.matrix(X)
y <- y[,1]
mod <- plsr(y~X, ncomp = 1)
beta <- mod$coefficients[,,1]
tau <- sd(mod$residuals)
p <- dim(X)[2]
####################################
SIhat <- cov(X)
shat <- cov(X,y)
sts <- as.numeric((t(shat)%*%shat))
Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p]
Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat
SI <- (Omegahat/sts)*(shat%*%t(shat)) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
####################################


####################################
# First Column, first row
set.seed(14)
reps <- 1000
n <- 10
w <- matrix(NA, nrow=reps,ncol=n)
Xg <- mvrnorm(n, apply(X,2,mean), Sigma = SI)
Xgc <- scale(Xg, center = TRUE, scale=FALSE)
Ytrue <- mean(y) + Xgc%*%beta
alpha <- 0.05
z <- qnorm(1-alpha/2)
IC.vec <- NULL


for (l in 1:reps){
  cat("l =",l,"\n")
  Ysimu <- rnorm(n, Ytrue, tau)
  IC.vec <- matrix(data = NA, nrow = n, ncol = 2)
  for (i in 1:n){
    Xl <- Xg[-i,]
    XN <- Xg[i,]
    yl <- Ysimu[-i]
    yN <- Ysimu[i]
    mod <- plsr(yl~Xl, 1)
    betahat <- mod$coefficients[,1,1]
    shat <- cov(Xl,yl)
    sts <- as.numeric(crossprod(shat,shat))
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
    SIhat <- cov(Xl)
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
    SIhat2 <- (Omegahat/sts)*tcrossprod(shat,shat) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
    sigmaY2hat <- var(yl)
    tau2hat <- var(mod$residuals)
    Cp2 <- tr(Omega0hat%*%Omega0hat)
    aux <- (1/n^2)*(tau2hat + var(yl)*Cp2/(Omegahat^2))
    Vhat <- as.numeric((Omegahat^(-2)/n)*t(XN-colMeans(Xl))%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%(XN-colMeans(Xl)))
    IC <- mean(yl) + as.numeric(t(betahat)%*%(XN-colMeans(Xl))) + c(-1,1)*z*(tau2hat/n + Vhat + aux)^(1/2)
    IC.vec[i,] <- IC
  }
  w[l,] <- howmany(IC.vec,Ytrue)$vec
}

output <- round(mean(w),3)
print(output)
####################################


####################################
# First Column, second row
reps <- 1000
n <- 100
w <- matrix(NA, nrow=reps,ncol=n)
Xg <- mvrnorm(n, apply(X,2,mean), Sigma = SI)
Xgc <- scale(Xg, center = TRUE, scale=FALSE)
Ytrue <- mean(y) + Xgc%*%beta
alpha <- 0.05
z <- qnorm(1-alpha/2)
IC.vec <- NULL


for (l in 1:reps){
  cat("l =",l,"\n")
  Ysimu <- rnorm(n, Ytrue, tau)
  IC.vec <- matrix(data = NA, nrow = n, ncol = 2)
  for (i in 1:n){
    Xl <- Xg[-i,]
    XN <- Xg[i,]
    yl <- Ysimu[-i]
    yN <- Ysimu[i]
    mod <- plsr(yl~Xl, 1)
    betahat <- mod$coefficients[,1,1]
    shat <- cov(Xl,yl)
    sts <- as.numeric(crossprod(shat,shat))
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
    SIhat <- cov(Xl)
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
    SIhat2 <- (Omegahat/sts)*tcrossprod(shat,shat) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
    sigmaY2hat <- var(yl)
    tau2hat <- var(mod$residuals)
    Cp2 <- tr(Omega0hat%*%Omega0hat)
    aux <- (1/n^2)*(tau2hat + var(yl)*Cp2/(Omegahat^2))
    Vhat <- as.numeric((Omegahat^(-2)/n)*t(XN-colMeans(Xl))%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%(XN-colMeans(Xl)))
    IC <- mean(yl) + as.numeric(t(betahat)%*%(XN-colMeans(Xl))) + c(-1,1)*z*(tau2hat/n + Vhat + aux)^(1/2)
    IC.vec[i,] <- IC
  }
  w[l,] <- howmany(IC.vec,Ytrue)$vec
}

output <- round(mean(w),3)
print(output)
####################################


####################################
# Second Column, first row
reps <- 1000
n <- 10
w <- matrix(NA, nrow=reps,ncol=n)
Xg <- mvrnorm(n, apply(X,2,mean), Sigma = SI)
Xgc <- scale(Xg, center = TRUE, scale=FALSE)
Ytrue <- mean(y) + Xgc%*%beta
alpha <- 0.05
z <- qnorm(1-alpha/2)
IC.vec <- NULL


for (l in 1:reps){
  cat("l =",l,"\n")
  Ysimu <- rnorm(n, Ytrue, tau)
  IC.vec <- matrix(data = NA, nrow = n, ncol = 2)
  for (i in 1:n){
    Xl <- Xg[-i,]
    XN <- Xg[i,]
    yl <- Ysimu[-i]
    yN <- Ysimu[i]
    mod <- plsr(yl~Xl, 1)
    betahat <- mod$coefficients[,1,1]
    shat <- cov(Xl,yl)
    sts <- as.numeric(crossprod(shat,shat))
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p]
    SIhat <- cov(Xl)
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat
    SIhat2 <- (Omegahat/sts)*tcrossprod(shat,shat) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
    sigmaY2hat <- var(yl)
    tau2hat <- var(mod$residuals)
    Cp2 <- tr(Omega0hat%*%Omega0hat)
    aux <- (1/n^2)*(tau2hat + var(yl)*Cp2/(Omegahat^2))
    Vhat <- as.numeric((Omegahat^(-2)/n)*t(XN-colMeans(Xl))%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%(XN-colMeans(Xl)))
    IC <- mean(yl) + as.numeric(t(betahat)%*%(XN-colMeans(Xl))) + c(-1,1)*z*(tau2hat/n + tau2hat + Vhat + aux)^(1/2)
    IC.vec[i,] <- IC
  }
  w[l,] <- howmany(IC.vec,Ysimu)$vec
}

output <- round(mean(w),3)
print(output)
####################################

####################################
# Second Column, second row
reps <- 1000
n <- 100
w <- matrix(NA, nrow=reps,ncol=n)
Xg <- mvrnorm(n, apply(X,2,mean), Sigma = SI)
Xgc <- scale(Xg, center = TRUE, scale=FALSE)
Ytrue <- mean(y) + Xgc%*%beta
alpha <- 0.05
z <- qnorm(1-alpha/2)
IC.vec <- NULL


for (l in 1:reps){
  cat("l =",l,"\n")
  Ysimu <- rnorm(n, Ytrue, tau)
  IC.vec <- matrix(data = NA, nrow = n, ncol = 2)
  for (i in 1:n){
    Xl <- Xg[-i,]
    XN <- Xg[i,]
    yl <- Ysimu[-i]
    yN <- Ysimu[i]
    mod <- plsr(yl~Xl, 1)
    betahat <- mod$coefficients[,1,1]
    shat <- cov(Xl,yl)
    sts <- as.numeric(crossprod(shat,shat))
    sigma.0hat <- qr.Q(qr(shat),complete=TRUE)[,2:p] # px(p-1)
    SIhat <- cov(Xl)
    Omegahat <- as.numeric(t(shat)%*%SIhat%*%shat/sts)
    Omega0hat <- t(sigma.0hat)%*%SIhat%*%sigma.0hat # (p-1)x(p-1)
    SIhat2 <- (Omegahat/sts)*tcrossprod(shat,shat) + sigma.0hat%*%Omega0hat%*%t(sigma.0hat)
    sigmaY2hat <- var(yl)
    tau2hat <- var(mod$residuals)
    Cp2 <- tr(Omega0hat%*%Omega0hat)
    aux <- (1/n^2)*(tau2hat + var(yl)*Cp2/(Omegahat^2))
    Vhat <- as.numeric((Omegahat^(-2)/n)*t(XN-colMeans(Xl))%*%(SIhat2*sigmaY2hat - shat%*%t(shat))%*%(XN-colMeans(Xl)))
    IC <- mean(yl) + as.numeric(t(betahat)%*%(XN-colMeans(Xl))) + c(-1,1)*z*(tau2hat/n + Vhat + aux)^(1/2)
    IC.vec[i,] <- IC
  }
  w[l,] <- howmany(IC.vec,Ytrue)$vec
}

output <- round(mean(w),3)
print(output)
####################################



