#Hypothesis: mu2 > {mu4, mu1} > mu3
#mu_1: DID patient
#mu_2: mimic normal
#mu_3: symptom simulated
#mu_4: true amnesic

#Required Packages & Data & script
rm(list = ls(all = TRUE))
library(mvtnorm)
library(coda)
library(MASS)
library(MCMCpack)
DID <- read.table("/.../DID.txt", header = FALSE, sep = "")
mcmc <- read.table("/.../mcmc.txt", header = FALSE, sep = "")
source("/.../PWMDE_pgm.R")

mcmc <- as.matrix(mcmc)
mcmc[,5] <- sqrt(mcmc[,5])
DID <- as.matrix(DID)

sum_y1 <- sum(DID[which(DID[,1]==1),2])
sum_y2 <- sum(DID[which(DID[,1]==2),2])
sum_y3 <- sum(DID[which(DID[,1]==3),2])
sum_y4 <- sum(DID[which(DID[,1]==4),2])
n1 <- length(which(DID[,1]==1))
n2 <- length(which(DID[,1]==2))
n3 <- length(which(DID[,1]==3))
n4 <- length(which(DID[,1]==4))

#parameter setting
#a is a_0 in Section 5
#nop: no. of interested points
#ncutslice: no. of partition
a <- 0.01
n <- dim(DID)[1]
nop <- 1000
ncutslice <- 20

mean_ref <- apply(mcmc, 2, mean)
sd_ref <- apply(mcmc, 2, sd)

#interested points for marginal posterior density estimation
targetpoint_mu2 <- seq(from=0.1, to=15, length.out=nop)
targetpoint_mu1 <- seq(from=0.1, to=15, length.out=nop)
targetpoint_mu4 <- seq(from=0.1, to=15, length.out=nop)
targetpoint_mu3 <- seq(from=0.1, to=15, length.out=nop)

##############################################################################
################################## mu1 #######################################
##############################################################################
#TT: the size of MCMC sample
TT <- 1000

######CMDE
den_est_cmde <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n1+0.19)*(mcmc[j, 1] - sum_y1/(n1+0.19))^2/(2*mcmc[j,5]^2)
  term1 <- pnorm(mcmc[j,2], mean=(sum_y1/(n1+0.19)), sd=(mcmc[j,5]/sqrt(n1+0.19)) )
  term2 <- pnorm(mcmc[j,3], mean=(sum_y1/(n1+0.19)), sd=(mcmc[j,5]/sqrt(n1+0.19)) )
  val[j] <- val[j] - log(term1 - term2)
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n1+0.19)))
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(targetpoint_mu1[1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1[i], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est_cmde[i] <- floating_control(val_, TT)  
}

######KDE
den_est_kde_y <- density(mcmc[1:TT,1], from=0.1, to=15,n=nop)$y
den_est_kde_x <- density(mcmc[1:TT,1], from=0.1, to=15,n=nop)$x

######PWMDE
inside <- 0
den_est <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu1 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu1 <- (slice_mu1[-ncutslice] + slice_mu1[-1])/2
  vol_mu1 <- diff(slice_mu1)
  position <- intersect(length(which(mcmc[j,1]>=slice_mu1)), c(1:(ncutslice-1)) )
  rep_updv_mu1 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu1[g] <- updv(rep_mu1[g], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu1) + rep_updv_mu1 )
    vo_den <- log(sum(exp(log(vol_mu1) + rep_updv_mu1 - lcons)) ) + lcons
    val[j] <- rep_updv_mu1[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1[1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1[i], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est[i] <- floating_control(val_, TT)  
}

######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu1 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu1 <- (slice_mu1[-ncutslice] + slice_mu1[-1])/2
  vol_mu1 <- diff(slice_mu1)
  position <- intersect(length(which(mcmc[j,1]>=slice_mu1)), c(1:(ncutslice-1)) )
  rep_updv_mu1 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu1[g] <- updv(rep_mu1[g], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu1) + rep_updv_mu1 )
    vo_den <- log(sum(exp(log(vol_mu1) + rep_updv_mu1 - lcons)) ) + lcons
    val[j] <- rep_updv_mu1[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1[1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est1[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1[i], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est1[i] <- floating_control(val_, TT)  
}

yticks <- seq(0, 1.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu1_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y, type="l", xlab=expression(mu[1]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu1_K20_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est), type="l", xlab=expression(mu[1]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu1_K50_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1), type="l", xlab=expression(mu[1]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

yticks <- seq(-0.5, 0.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu1_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y-exp(den_est_cmde), type="l", xlab=expression(mu[1]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu1_K20_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est)-exp(den_est_cmde), type="l", xlab=expression(mu[1]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu1_K50_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1)-exp(den_est_cmde), type="l", xlab=expression(mu[1]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

###L-distance
VA <- diff(den_est_kde_x)[1]
LD_kde <- sum(abs(den_est_kde_y-exp(den_est_cmde) )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#[1] 0.065 0.014 0.001


##############################################################################
################################## mu2 #######################################
##############################################################################
######CMDE
den_est_cmde <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n2+0.25)*(mcmc[j, 2] - sum_y2/(n2+0.25))^2/(2*mcmc[j,5]^2)
  val[j] <- val[j] - log(1-pnorm(max(c(mcmc[j,4],mcmc[j,1])),mean=(sum_y2/(n2+0.25)), sd=(mcmc[j,5]/sqrt(n2+0.25)) ))
  print(1-pnorm(max(c(mcmc[j,4],mcmc[j,1])),mean=(sum_y2/(n2+0.25)), sd=(mcmc[j,5]/sqrt(n2+0.25)) ))
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n2+0.25)))
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[1], mcmc[j,3], mcmc[j,4], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[i], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est_cmde[i] <- floating_control(val_, TT)  
}

######KDE
den_est_kde_y <- density(mcmc[1:TT,2], from=0.1, to=15,n=nop)$y
den_est_kde_x <- density(mcmc[1:TT,2], from=0.1, to=15,n=nop)$x

######PWMDE
ncutslice <- 20
inside <- 0
den_est <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  lb <- max(mcmc[j,1],mcmc[j,4])
  ub <- 15
  slice_mu2 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu2 <- (slice_mu2[-ncutslice] + slice_mu2[-1])/2
  vol_mu2 <- diff(slice_mu2)
  position <- intersect(length(which(mcmc[j,2]>=slice_mu2)), c(1:(ncutslice-1)) )
  rep_updv_mu2 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu2[g] <- updv(mcmc[j,1], rep_mu2[g], mcmc[j,3], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu2) + rep_updv_mu2 )
    vo_den <- log(sum(exp(log(vol_mu2) + rep_updv_mu2 - lcons)) ) + lcons
    val[j] <- rep_updv_mu2[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[1], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[i], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est[i] <- floating_control(val_, TT)  
}

######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  lb <- max(mcmc[j,1],mcmc[j,4])
  ub <- 15
  slice_mu2 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu2 <- (slice_mu2[-ncutslice] + slice_mu2[-1])/2
  vol_mu2 <- diff(slice_mu2)
  position <- intersect(length(which(mcmc[j,2]>=slice_mu2)), c(1:(ncutslice-1)) )
  rep_updv_mu2 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu2[g] <- updv(mcmc[j,1], rep_mu2[g], mcmc[j,3], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu2) + rep_updv_mu2 )
    vo_den <- log(sum(exp(log(vol_mu2) + rep_updv_mu2 - lcons)) ) + lcons
    val[j] <- rep_updv_mu2[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[1], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est1[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], targetpoint_mu2[i], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  }
  den_est1[i] <- floating_control(val_, TT)  
}
yticks <- seq(0, 1.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu2_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y, type="l", xlab=expression(mu[2]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu2_K20_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est), type="l", xlab=expression(mu[2]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu2_K50_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1), type="l", xlab=expression(mu[2]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

yticks <- seq(-0.5, 0.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu2_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y-exp(den_est_cmde), type="l", xlab=expression(mu[2]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu2_K20_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est)-exp(den_est_cmde), type="l", xlab=expression(mu[2]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu2_K50_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1)-exp(den_est_cmde), type="l", xlab=expression(mu[2]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

###L-distance
VA <- diff(den_est_kde_x)[1]
LD_kde <- sum(abs(den_est_kde_y-exp(den_est_cmde) )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#[1] 0.069 0.002 0.006


##############################################################################
################################## mu3 #######################################
##############################################################################
######CMDE
den_est_cmde <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n3+0.25)*(mcmc[j, 3] - sum_y3/(n3+0.25))^2/(2*mcmc[j,5]^2)
  val[j] <- val[j] - log(pnorm(min(c(mcmc[j,4],mcmc[j,1])),mean=(sum_y3/(n3+0.25)), sd=(mcmc[j,5]/sqrt(n3+0.25)) ))
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n3+0.25)))
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[1], mcmc[j,4], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[i], mcmc[j,4], mcmc[j,5])
  }
  den_est_cmde[i] <- floating_control(val_, TT)  
}

######KDE
den_est_kde_y <- density(mcmc[1:TT,3], from=0.1, to=15,n=nop)$y
den_est_kde_x <- density(mcmc[1:TT,3], from=0.1, to=15,n=nop)$x

######PWMDE
ncutslice <- 20
inside <- 0
den_est <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- min(mcmc[j,1],mcmc[j,4])
  lb <- 0
  slice_mu3 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu3 <- (slice_mu3[-ncutslice] + slice_mu3[-1])/2
  vol_mu3 <- diff(slice_mu3)
  position <- intersect(length(which(mcmc[j,3]>=slice_mu3)), c(1:(ncutslice-1)) )
  rep_updv_mu3 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu3[g] <- updv(mcmc[j,1], mcmc[j,2], rep_mu3[g], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu3) + rep_updv_mu3 )
    vo_den <- log(sum(exp(log(vol_mu3) + rep_updv_mu3 - lcons)) ) + lcons
    val[j] <- rep_updv_mu3[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[1], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[i], mcmc[j,4], mcmc[j,5])
  }
  den_est[i] <- floating_control(val_, TT)  
}

######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- min(mcmc[j,1],mcmc[j,4])
  lb <- 0
  slice_mu3 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu3 <- (slice_mu3[-ncutslice] + slice_mu3[-1])/2
  vol_mu3 <- diff(slice_mu3)
  position <- intersect(length(which(mcmc[j,3]>=slice_mu3)), c(1:(ncutslice-1)) )
  rep_updv_mu3 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu3[g] <- updv(mcmc[j,1], mcmc[j,2], rep_mu3[g], mcmc[j,4], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu3) + rep_updv_mu3 )
    vo_den <- log(sum(exp(log(vol_mu3) + rep_updv_mu3 - lcons)) ) + lcons
    val[j] <- rep_updv_mu3[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[1], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est1[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], targetpoint_mu3[i], mcmc[j,4], mcmc[j,5])
  }
  den_est1[i] <- floating_control(val_, TT)  
}
yticks <- seq(0, 1.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu3_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y, type="l", xlab=expression(mu[3]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu3_K20_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est), type="l", xlab=expression(mu[3]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu3_K50_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1), type="l", xlab=expression(mu[3]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

yticks <- seq(-0.5, 0.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu3_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y-exp(den_est_cmde), type="l", xlab=expression(mu[3]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu3_K20_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est)-exp(den_est_cmde), type="l", xlab=expression(mu[3]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu3_K50_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1)-exp(den_est_cmde), type="l", xlab=expression(mu[3]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

###L-distance
VA <- diff(den_est_kde_x)[1]
LD_kde <- sum(abs(den_est_kde_y-exp(den_est_cmde) )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#[1] 0.055 0.000 0.000


##############################################################################
################################## mu4 #######################################
##############################################################################
######CMDE
den_est_cmde <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n4+0.25)*(mcmc[j, 4] - sum_y4/(n4+0.25))^2/(2*mcmc[j,5]^2)
  term1 <- pnorm(mcmc[j,2], mean=(sum_y4/(n4+0.25)), sd=(mcmc[j,5]/sqrt(n4+0.25)) )
  term2 <- pnorm(mcmc[j,3], mean=(sum_y4/(n4+0.25)), sd=(mcmc[j,5]/sqrt(n4+0.25)) )
  val[j] <- val[j] - log(term1 - term2)
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n4+0.25)))
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[1], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[i], mcmc[j,5])
  }
  den_est_cmde[i] <- floating_control(val_, TT)  
}

######KDE
den_est_kde_y <- density(mcmc[1:TT,4], from=0.1, to=15,n=nop)$y
den_est_kde_x <- density(mcmc[1:TT,4], from=0.1, to=15,n=nop)$x

######PWMDE
ncutslice <- 20
inside <- 0
den_est <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu4 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu4 <- (slice_mu4[-ncutslice] + slice_mu4[-1])/2
  vol_mu4 <- diff(slice_mu4)
  position <- intersect(length(which(mcmc[j,4]>=slice_mu4)), c(1:(ncutslice-1)) )
  rep_updv_mu4 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu4[g] <- updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], rep_mu4[g], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu4) + rep_updv_mu4 )
    vo_den <- log(sum(exp(log(vol_mu4) + rep_updv_mu4 - lcons)) ) + lcons
    val[j] <- rep_updv_mu4[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[1], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[i], mcmc[j,5])
  }
  den_est[i] <- floating_control(val_, TT)  
}

######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(NA, nop)
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu4 <- seq(from=lb, to=ub, length.out=ncutslice)
  rep_mu4 <- (slice_mu4[-ncutslice] + slice_mu4[-1])/2
  vol_mu4 <- diff(slice_mu4)
  position <- intersect(length(which(mcmc[j,4]>=slice_mu4)), c(1:(ncutslice-1)) )
  rep_updv_mu4 <- rep(NA, ncutslice-1)
  if (sum(position)>0){
    for (g in 1:(ncutslice-1)){
      rep_updv_mu4[g] <- updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], rep_mu4[g], mcmc[j,5] )
    }
    lcons <- max(log(vol_mu4) + rep_updv_mu4 )
    vo_den <- log(sum(exp(log(vol_mu4) + rep_updv_mu4 - lcons)) ) + lcons
    val[j] <- rep_updv_mu4[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[1], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est1[1] <- floating_control(val_, TT)
for (i in 2:nop) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], targetpoint_mu4[i], mcmc[j,5])
  }
  den_est1[i] <- floating_control(val_, TT)  
}
yticks <- seq(0, 1.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu4_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y, type="l", xlab=expression(mu[4]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu4_K20_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est), type="l", xlab=expression(mu[4]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu4_K50_a.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1), type="l", xlab=expression(mu[4]), ylab="", ylim=c(0,1.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()

yticks <- seq(-0.5, 0.5, by=0.5)
xticks <- seq(0, 15, by=5)
pdf("KDE_mu4_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, den_est_kde_y-exp(den_est_cmde), type="l", xlab=expression(mu[4]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
legend(13, 0.5, "L-distance=0.082")
dev.off()
pdf("PWMDE_mu4_K20_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est)-exp(den_est_cmde), type="l", xlab=expression(mu[4]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()
pdf("PWMDE_mu4_K50_d.pdf")
par(mar=c(6,5,1,1))
plot(den_est_kde_x, exp(den_est1)-exp(den_est_cmde), type="l", xlab=expression(mu[4]), ylab="", ylim=c(-0.5,0.5), xlim=c(0,15), axes=F, cex.lab = 3.5, lwd=5, frame.plot = T)
axis(2, at =yticks, las=2, cex.axis = 2.3)
axis(1, at =xticks, las=1, cex.axis = 2.3)
dev.off()


###L-distance
VA <- diff(den_est_kde_x)[1]
LD_kde <- sum(abs(den_est_kde_y-exp(den_est_cmde) )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#[1] 0.082 0.029 0.001


##############################################################################
################################## mu1 mu4 ###################################
##############################################################################
nop <- 100
targetpoint_mu1 <- seq(from=1, to=5, length.out=nop)
targetpoint_mu4 <- seq(from=3, to=6, length.out=nop)

targetpoint_mu1mu4 <- NULL
for (i in 1:nop){
  temp <- cbind(rep(targetpoint_mu1[i], nop), targetpoint_mu4)
  targetpoint_mu1mu4 <- rbind(targetpoint_mu1mu4, temp)
}
#######CMDE
den_est_cmde <- rep(NA, dim(targetpoint_mu1mu4)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n4+0.25)*(mcmc[j, 4] - sum_y4/(n4+0.25))^2/(2*mcmc[j,5]^2)
  val[j] <- val[j] -(n1+0.19)*(mcmc[j, 1] - sum_y1/(n1+0.19))^2/(2*mcmc[j,5]^2)
  
  term1 <- pnorm(mcmc[j,2], mean=(sum_y4/(n4+0.25)), sd=(mcmc[j,5]/sqrt(n4+0.25)) )
  term2 <- pnorm(mcmc[j,3], mean=(sum_y4/(n4+0.25)), sd=(mcmc[j,5]/sqrt(n4+0.25)) )
  val[j] <- val[j] - log(term1 - term2)
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n4+0.25)))  
  print(log(term1 - term2))
  term1 <- pnorm(mcmc[j,2], mean=(sum_y1/(n1+0.19)), sd=(mcmc[j,5]/sqrt(n1+0.19)) )
  term2 <- pnorm(mcmc[j,3], mean=(sum_y1/(n1+0.19)), sd=(mcmc[j,5]/sqrt(n1+0.19)) )
  val[j] <- val[j] - log(term1 - term2)
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n1+0.19)))  
  
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(targetpoint_mu1mu4[1, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[1, 2], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:dim(targetpoint_mu1mu4)[1]) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1mu4[i, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[i, 2], mcmc[j,5])
  }
  den_est_cmde[i] <- floating_control(val_, TT)  
}

######KDE
kkpure <- cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], c(kde2d(mcmc[1:TT,4], mcmc[1:TT,1], n = 100, lims = c(3,6, 1,5))$z) )
kkdif <- cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], c(kde2d(mcmc[1:TT,4], mcmc[1:TT,1], n = 100, lims = c(3,6, 1,5))$z)-exp(den_est_cmde) )
write.table(kkpure, file = ".../pgm/temp.csv", sep=',')
write.table(kkdif, file = ".../pgm/temp_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], kkdif[,3], type = "p")

#######PWMDE
ncutslice <- 20
inside <- 0
den_est <- rep(NA, dim(targetpoint_mu1mu4)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu4 <- seq(from=lb, to=ub, length.out=ncutslice)
  slice_mu1 <- seq(from=lb, to=ub, length.out=ncutslice)
  
  rep_mu4 <- (slice_mu4[-ncutslice] + slice_mu4[-1])/2
  rep_mu1 <- (slice_mu1[-ncutslice] + slice_mu1[-1])/2
  vol_mu4 <- diff(slice_mu4)[1]
  vol_mu1 <- diff(slice_mu1)[1]
  position1 <- intersect(length(which(mcmc[j,4]>=slice_mu4)), c(1:(ncutslice-1)) )
  position2 <- intersect(length(which(mcmc[j,1]>=slice_mu1)), c(1:(ncutslice-1)) )
  position <- (position2-1)*(ncutslice-1) + position1

  rep_updv <- NULL
  if (sum(position1)>0 & sum(position2)>0){
    for (g1 in 1:(ncutslice-1)){
      for (g2 in 1:(ncutslice-1)){
        rep_updv <- c(rep_updv, updv(rep_mu1[g1], mcmc[j,2], mcmc[j,3], rep_mu4[g2], mcmc[j,5] ) )
      }
    }  
    lcons <- max(log(vol_mu4) + log(vol_mu1) + rep_updv )
    vo_den <- log(sum(exp(log(vol_mu4) + log(vol_mu1) + rep_updv - lcons)) ) + lcons
    val[j] <- rep_updv[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1mu4[1, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[1, 2], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est[1] <- floating_control(val_, TT) 

for (i in 2:dim(targetpoint_mu1mu4)[1]) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1mu4[i, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[i, 2], mcmc[j,5])
  }
  den_est[i] <- floating_control(val_, TT)  
}
write.table(cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], exp(den_est) ), file = ".../pgm/temp1.csv", sep=',')
write.table(cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], exp(den_est)-exp(den_est_cmde) ), file = ".../pgm/temp1_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], (exp(den_est)-exp(den_est_cmde)), type = "p")

#######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(NA, dim(targetpoint_mu1mu4)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub <- mcmc[j,2]
  lb <- mcmc[j,3]
  slice_mu4 <- seq(from=lb, to=ub, length.out=ncutslice)
  slice_mu1 <- seq(from=lb, to=ub, length.out=ncutslice)
  
  rep_mu4 <- (slice_mu4[-ncutslice] + slice_mu4[-1])/2
  rep_mu1 <- (slice_mu1[-ncutslice] + slice_mu1[-1])/2
  vol_mu4 <- diff(slice_mu4)[1]
  vol_mu1 <- diff(slice_mu1)[1]
  position1 <- intersect(length(which(mcmc[j,4]>=slice_mu4)), c(1:(ncutslice-1)) )
  position2 <- intersect(length(which(mcmc[j,1]>=slice_mu1)), c(1:(ncutslice-1)) )
  position <- (position2-1)*(ncutslice-1) + position1
  
  rep_updv <- NULL
  if (sum(position1)>0 & sum(position2)>0){
    for (g1 in 1:(ncutslice-1)){
      for (g2 in 1:(ncutslice-1)){
        rep_updv <- c(rep_updv, updv(rep_mu1[g1], mcmc[j,2], mcmc[j,3], rep_mu4[g2], mcmc[j,5] ) )
      }
    }  
    lcons <- max(log(vol_mu4) + log(vol_mu1) + rep_updv )
    vo_den <- log(sum(exp(log(vol_mu4) + log(vol_mu1) + rep_updv - lcons)) ) + lcons
    val[j] <- rep_updv[position] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1mu4[1, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[1, 2], mcmc[j,5])
    inside <- inside + 1
  }
}
den_est1[1] <- floating_control(val_, TT) 

for (i in 2:dim(targetpoint_mu1mu4)[1]) {
  for (j in 1:TT){
    val_[j] <- val[j] + updv(targetpoint_mu1mu4[i, 1], mcmc[j,2], mcmc[j,3], targetpoint_mu1mu4[i, 2], mcmc[j,5])
  }
  den_est1[i] <- floating_control(val_, TT)  
}
write.table(cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], exp(den_est1) ), file = ".../pgm/temp2.csv", sep=',')
write.table(cbind(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], exp(den_est1)-exp(den_est_cmde) ), file = ".../pgm/temp2_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu4[, 1], targetpoint_mu1mu4[, 2], (exp(den_est1)-exp(den_est_cmde)), type = "p")

###L-distance
VA <- diff(targetpoint_mu1)[1]*diff(targetpoint_mu4 )[1]
LD_kde <- sum(abs(kkdif[, 3] )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#[1] 0.174 0.016 0.000


##############################################################################
################################## mu1 mu3 ###################################
##############################################################################
nop <- 100
targetpoint_mu1 <- seq(from=0, to=5, length.out=nop)
targetpoint_mu3 <- seq(from=0, to=5, length.out=nop)

targetpoint_mu1mu3 <- NULL
for (i in 1:nop){
  temp <- cbind(rep(targetpoint_mu1[i], nop), targetpoint_mu3)
  targetpoint_mu1mu3 <- rbind(targetpoint_mu1mu3, temp)
}
######CMDE
den1 <- function(mu1){
  denin <- NA
  meanformu3 <- sum_y3/(n3+0.25)
  sdformu3 <- mcmc[j,5]*sqrt(1/(n3+0.25))
  meanformu1 <- sum_y1/(n1+0.19)
  sdformu1 <- mcmc[j,5]*sqrt(1/(n1+0.19))
  
  part1formu3 <- pnorm(mcmc[j,4], mean=meanformu3, sd=sdformu3 )
  denin <- dnorm(mu1, mean=meanformu1, sd=sdformu1 )*part1formu3
}
den2 <- function(mu1){
  denin <- NA
  meanformu3 <- sum_y3/(n3+0.25)
  sdformu3 <- mcmc[j,5]*sqrt(1/(n3+0.25))
  meanformu1 <- sum_y1/(n1+0.19)
  sdformu1 <- mcmc[j,5]*sqrt(1/(n1+0.19))
  
  part1formu3 <- pnorm(mu1, mean=meanformu3, sd=sdformu3 )
  denin <- dnorm(mu1, mean=meanformu1, sd=sdformu1 )*part1formu3
}
den_est_cmde <- rep(-Inf, dim(targetpoint_mu1mu3)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  val[j] <- -(n3+0.25)*(mcmc[j, 3] - sum_y3/(n3+0.25))^2/(2*mcmc[j,5]^2)
  val[j] <- val[j] -(n1+0.19)*(mcmc[j, 1] - sum_y1/(n1+0.19))^2/(2*mcmc[j,5]^2)
  
  val[j] <- val[j] - log(integrate(den1, mcmc[j,4], mcmc[j,2])$value + integrate(den2, -Inf, mcmc[j,4])$value )
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n3+0.25)))  
  val[j] <- val[j] - log(mcmc[j,5]*sqrt(2*pi/(n1+0.19)))  
  
  val[j] <- val[j] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
  val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
}
den_est_cmde[1] <- floating_control(val_, TT)
for (i in 2:dim(targetpoint_mu1mu3)[1]) {
  if (targetpoint_mu1mu3[i, 1]>targetpoint_mu1mu3[i, 2]){
    for (j in 1:TT){
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[i, 1], mcmc[j,2], targetpoint_mu1mu3[i, 2], mcmc[j,4], mcmc[j,5])
    }
    den_est_cmde[i] <- floating_control(val_, TT)  
  }
}

######KDE
kkpure <- cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], c(kde2d(mcmc[1:TT,3], mcmc[1:TT,1], n = 100, lims = c(0,5, 0,5))$z) )
kkdif <- cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], c(kde2d(mcmc[1:TT,3], mcmc[1:TT,1], n = 100, lims = c(0,5, 0,5))$z)-exp(den_est_cmde) )
mark0 <- which(kkpure[,1]<=kkpure[,2])
kkpure[mark0,3] <- 0
kkdif[mark0,3] <- 0
write.table(kkpure, file = ".../pgm/temp.csv", sep=',')
write.table(kkdif, file = ".../pgm/temp_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], kkdif[,3], type = "p")

######PWMDE
ncutslice <- 20
inside <- 0
den_est <- rep(-Inf, dim(targetpoint_mu1mu3)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub1_1 <- mcmc[j,2]
  lb1_1 <- mcmc[j,4]
  ub3_1 <- mcmc[j,4]
  lb3_1 <- 0
  
  #rectangular
  slice_mu3_1 <- seq(from=lb3_1, to=ub3_1, length.out=ncutslice)
  slice_mu1_1 <- seq(from=lb1_1, to=ub1_1, length.out=ncutslice)
  
  rep_mu3_1 <- (slice_mu3_1[-ncutslice] + slice_mu3_1[-1])/2
  rep_mu1_1 <- (slice_mu1_1[-ncutslice] + slice_mu1_1[-1])/2
  vol_mu3_1 <- diff(slice_mu3_1)[1]
  vol_mu1_1 <- diff(slice_mu1_1)[1]
  
  rep_updv_1 <- NULL
  for (g1 in 1:(ncutslice-1)){
    for (g2 in 1:(ncutslice-1)){
      rep_updv_1 <- c(rep_updv_1, updv(rep_mu1_1[g1], mcmc[j,2], rep_mu3_1[g2], mcmc[j,4], mcmc[j,5] ) )
    }
  }  
  #triangular  
  ub1_2 <- mcmc[j,4]
  lb1_2 <- 0
  ub3_2 <- mcmc[j,4]
  lb3_2 <- 0
  slice_mu3_2 <- seq(from=lb3_2, to=ub3_2, length.out=ncutslice)
  slice_mu1_2 <- seq(from=lb1_2, to=ub1_2, length.out=ncutslice)
  slice_mu3_2 <- slice_mu3_2[-ncutslice] 
  slice_mu1_2 <- slice_mu1_2[-1] 
  rep_mu3_2 <- (slice_mu3_2[-(ncutslice-1)] + slice_mu3_2[-1])/2
  rep_mu1_2 <- (slice_mu1_2[-(ncutslice-1)] + slice_mu1_2[-1])/2
  
  vol_mu3_2 <- diff(slice_mu3_2)[1]
  vol_mu1_2 <- diff(slice_mu1_2)[1]
  
  rep_updv_2 <- NULL
  for (g1 in 1:(ncutslice-2)){
    for (g2 in g1:(ncutslice-2)){
      rep_updv_2 <- c(rep_updv_2, updv(rep_mu1_2[g2], mcmc[j,2], rep_mu3_2[g1], mcmc[j,4], mcmc[j,5] ) )
    }
  }  
  
  v_rect <- log(vol_mu3_1) + log(vol_mu1_1) + rep_updv_1
  v_triag <- log(vol_mu3_2) + log(vol_mu1_2) + rep_updv_2
  lcons <- max( c(v_rect, v_triag) )
  vo_den <- log(sum(exp(c(v_rect, v_triag) - lcons) ) ) + lcons

  if (mcmc[j,1]>mcmc[j,4]){
    position1 <- intersect(length(which(mcmc[j,3]>=slice_mu3_1)), c(1:(ncutslice-1)) )
    position2 <- intersect(length(which(mcmc[j,1]>=slice_mu1_1)), c(1:(ncutslice-1)) )
    position_1 <- (position2-1)*(ncutslice-1) + position1
    
    val[j] <- rep_updv_1[position_1] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
  else {
    pa1 <- intersect(which(mcmc[j,1]>=slice_mu1_2), c(1:(ncutslice-2)) ) 
    lpa1 <- length(pa1)
    pa2 <- intersect(which(mcmc[j,3]>=slice_mu3_2), c(1:(ncutslice-2)) ) 
    lpa2 <- length(pa2)
    if (lpa1 >= lpa2 & lpa2>1){
      position_2 <- (ncutslice-1)*(lpa2-1) - sum(1:(lpa2-1) ) + lpa1 - (lpa2-1)
      val[j] <- rep_updv_2[position_2] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
      val[j] <- val[j] - vo_den 
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
      inside <- inside + 1
    }
    if (lpa1 >= lpa2 & lpa2==1){
      position_2 <- lpa1
      
      val[j] <- rep_updv_2[position_2] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
      val[j] <- val[j] - vo_den 
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
      inside <- inside + 1
    }
  }
}
den_est[1] <- floating_control(val_, TT) 

for (i in 2:dim(targetpoint_mu1mu3)[1]) {
  if (targetpoint_mu1mu3[i, 1]>targetpoint_mu1mu3[i, 2]){
    for (j in 1:TT){
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[i, 1], mcmc[j,2], targetpoint_mu1mu3[i, 2], mcmc[j,4], mcmc[j,5])
    }
    den_est[i] <- floating_control(val_, TT)  
  }
  print(i)
}
write.table(cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], exp(den_est) ), file = ".../pgm/temp1.csv", sep=',')
write.table(cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], exp(den_est)-exp(den_est_cmde) ), file = ".../pgm/temp1_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], (exp(den_est)-exp(den_est_cmde)), type = "p")

######PWMDE
ncutslice <- 50
inside <- 0
den_est1 <- rep(-Inf, dim(targetpoint_mu1mu3)[1])
val <- rep(NA, TT)
val_ <- rep(NA, TT)
for (j in 1:TT){
  ub1_1 <- mcmc[j,2]
  lb1_1 <- mcmc[j,4]
  ub3_1 <- mcmc[j,4]
  lb3_1 <- 0
  
  #rectangular
  slice_mu3_1 <- seq(from=lb3_1, to=ub3_1, length.out=ncutslice)
  slice_mu1_1 <- seq(from=lb1_1, to=ub1_1, length.out=ncutslice)
  
  rep_mu3_1 <- (slice_mu3_1[-ncutslice] + slice_mu3_1[-1])/2
  rep_mu1_1 <- (slice_mu1_1[-ncutslice] + slice_mu1_1[-1])/2
  vol_mu3_1 <- diff(slice_mu3_1)[1]
  vol_mu1_1 <- diff(slice_mu1_1)[1]
  
  rep_updv_1 <- NULL
  for (g1 in 1:(ncutslice-1)){
    for (g2 in 1:(ncutslice-1)){
      rep_updv_1 <- c(rep_updv_1, updv(rep_mu1_1[g1], mcmc[j,2], rep_mu3_1[g2], mcmc[j,4], mcmc[j,5] ) )
    }
  }  
  #triangular  
  ub1_2 <- mcmc[j,4]
  lb1_2 <- 0
  ub3_2 <- mcmc[j,4]
  lb3_2 <- 0
  slice_mu3_2 <- seq(from=lb3_2, to=ub3_2, length.out=ncutslice)
  slice_mu1_2 <- seq(from=lb1_2, to=ub1_2, length.out=ncutslice)
  slice_mu3_2 <- slice_mu3_2[-ncutslice] 
  slice_mu1_2 <- slice_mu1_2[-1] 
  rep_mu3_2 <- (slice_mu3_2[-(ncutslice-1)] + slice_mu3_2[-1])/2
  rep_mu1_2 <- (slice_mu1_2[-(ncutslice-1)] + slice_mu1_2[-1])/2
  
  vol_mu3_2 <- diff(slice_mu3_2)[1]
  vol_mu1_2 <- diff(slice_mu1_2)[1]
  
  rep_updv_2 <- NULL
  for (g1 in 1:(ncutslice-2)){
    for (g2 in g1:(ncutslice-2)){
      rep_updv_2 <- c(rep_updv_2, updv(rep_mu1_2[g2], mcmc[j,2], rep_mu3_2[g1], mcmc[j,4], mcmc[j,5] ) )
    }
  }  
  
  v_rect <- log(vol_mu3_1) + log(vol_mu1_1) + rep_updv_1
  v_triag <- log(vol_mu3_2) + log(vol_mu1_2) + rep_updv_2
  lcons <- max( c(v_rect, v_triag) )
  vo_den <- log(sum(exp(c(v_rect, v_triag) - lcons) ) ) + lcons

  if (mcmc[j,1]>mcmc[j,4]){
    position1 <- intersect(length(which(mcmc[j,3]>=slice_mu3_1)), c(1:(ncutslice-1)) )
    position2 <- intersect(length(which(mcmc[j,1]>=slice_mu1_1)), c(1:(ncutslice-1)) )
    position_1 <- (position2-1)*(ncutslice-1) + position1
    
    val[j] <- rep_updv_1[position_1] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
    val[j] <- val[j] - vo_den 
    val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
    inside <- inside + 1
  }
  else {
    pa1 <- intersect(which(mcmc[j,1]>=slice_mu1_2), c(1:(ncutslice-2)) ) 
    lpa1 <- length(pa1)
    pa2 <- intersect(which(mcmc[j,3]>=slice_mu3_2), c(1:(ncutslice-2)) ) 
    lpa2 <- length(pa2)
    if (lpa1 >= lpa2 & lpa2>1){
      position_2 <- (ncutslice-1)*(lpa2-1) - sum(1:(lpa2-1) ) + lpa1 - (lpa2-1)
      val[j] <- rep_updv_2[position_2] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
      val[j] <- val[j] - vo_den 
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
      inside <- inside + 1
    }
    if (lpa1 >= lpa2 & lpa2==1){
      position_2 <- lpa1
      
      val[j] <- rep_updv_2[position_2] - updv(mcmc[j,1], mcmc[j,2], mcmc[j,3], mcmc[j,4], mcmc[j,5])
      val[j] <- val[j] - vo_den 
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[1, 1], mcmc[j,2], targetpoint_mu1mu3[1, 2], mcmc[j,4], mcmc[j,5])
      inside <- inside + 1
    }
  }
}
den_est1[1] <- floating_control(val_, TT) 

for (i in 2:dim(targetpoint_mu1mu3)[1]) {
  if (targetpoint_mu1mu3[i, 1]>targetpoint_mu1mu3[i, 2]){
    for (j in 1:TT){
      val_[j] <- val[j] + updv(targetpoint_mu1mu3[i, 1], mcmc[j,2], targetpoint_mu1mu3[i, 2], mcmc[j,4], mcmc[j,5])
    }
    den_est1[i] <- floating_control(val_, TT)  
  }
  print(i)
}
write.table(cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], exp(den_est1) ), file = ".../pgm/temp2.csv", sep=',')
write.table(cbind(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], exp(den_est1)-exp(den_est_cmde) ), file = ".../pgm/temp2_d.csv", sep=',')
#scatterplot3d(targetpoint_mu1mu3[, 1], targetpoint_mu1mu3[, 2], (exp(den_est)-exp(den_est_cmde)), type = "p")

###L-distance
VA <- diff(targetpoint_mu1)[1]*diff(targetpoint_mu3 )[1]
LD_kde <- sum(abs(kkdif[, 3] )*VA )
LD_PWMDE <- sum(abs(exp(den_est)-exp(den_est_cmde) )*VA )
LD_PWMDE1 <- sum(abs(exp(den_est1)-exp(den_est_cmde) )*VA )
print(round(c(LD_kde, LD_PWMDE, LD_PWMDE1), 3) )
#0.141 0.003 0.002
