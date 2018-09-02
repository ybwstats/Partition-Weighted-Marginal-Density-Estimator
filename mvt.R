##Required packages
rm(list = ls(all = TRUE))
library(mvtnorm)

######################################################################
###Functions
#purpose: numerical control when summing up exp() terms
#input: "ratios" includes all MCMC sample ratios, which have been assigned to weighted averages
#input: "tot" is the size of an MCMC sample
#output is the marginal posterior density estimate in log scale
control <- function(ratios, tot)
{
  b <- ratios
  est_d_r <- 0
  b.max <- max(b)
  est_d_r <- log(sum(exp(b-b.max))) + b.max - log(tot) 
  return(est_d_r)
}

#purpose: forming the partition subsets to get a weighted average
#input: "paramcmc" is a point of the MCMC sample + its position 
#input: "nslice" is no. of partition subsets
#output is a weighted average: q(\theta_k^*,\xi_t)/\sum{q(\theta_k^*,\xi_t)V(A_k(\xi_t)) } 
LOR_partition <- function(paramcmc, nslice){
  loc_ma <- paramcmc[1]
  paramcmc <- paramcmc[-1]

  rpoints <- matrix(paramcmc, nrow = nslice, ncol = length(paramcmc), byrow = TRUE)   #matrix(c(1,2,3), nrow=100, ncol=3, byrow=T)
  rpoints[, 1:nofp] <- matrix(M_mean[1:nofp], nrow = nslice, ncol = nofp, byrow = TRUE) + reprp%*%scale_vec
  kreprp <- apply(rpoints, 1, dmvt, delta = rep(0, length(paramcmc) ), sigma = sigma_, df = df_, log = TRUE, type = "shifted")
   
  rvoll <- rvol + kreprp
  rvoll <- log(sum(exp(rvoll - max(rvoll) ) ) ) + max(rvoll)  
  result <- kreprp[loc_ma] - rvoll 
  return(result)
}

######################################################################
###Settings
#DUMZ: a random seed, which will be assigned as a value from 1 to 100
#nofp: no. of focused points 
#nonfp: no. of non-focused points
#tpoint: the point of interest on the marginal density
#rcover: the specified radius value forming the conditional working parameter space
#delta_: location parameter of the multivariate t distribution
#sigma_: scale parameter (covariance) of the multivariate t distribution
#df_: degree of freedom of the multivariate t distribution
#TT: the size of an MCMC sample
#nos: no. of simulation repeats under the chosen random seed
set.seed(DUMZ)
setwd("/.../")
nofp <- 3
nonfp <- 2
tpoint <- rep(0, nofp)
rcover <- round(sqrt(qchisq(0.95, nofp) ), 1)
delta_ <- rep(0, (nofp+nonfp) )
sigma_ <- NULL
for (i in 1:(nofp+nonfp) ){
  sigma_ <- rbind(sigma_, abs((1-i):(nofp+nonfp-i)) )
}  
sigma_ <- 3*0.5^sigma_
df_ <- 3
TT <- 10000
nos <- 10

######################################################################
###Computation & Results
#den_est: the IWMDE estimate
#den_est_pwmde: the PWMDE estimate
#iwmdetime: IWMDE CPU running time for each simulation repeat 
#pwmdetime: PWMDE CPU running time for each simulation repeat 
den_est <- rep(NA, nos)
den_est_pwmde1 <- rep(NA, nos)
den_est_pwmde2 <- rep(NA, nos)
den_est_pwmde3 <- rep(NA, nos)
den_est_pwmde4 <- rep(NA, nos)
iwmdetime <- NULL
pwmdetime1 <- NULL
pwmdetime2 <- NULL
pwmdetime3 <- NULL
pwmdetime4 <- NULL
for (k in 1:nos){
  mcmc <- rmvt(n=TT, sigma = sigma_, df = df_, delta = delta_, type = "shifted")
  mcmc <- as.matrix(mcmc, ncol=(nofp+nonfp), byrow=TRUE)
  target <- mcmc
  target[, 1:nofp] <- matrix(tpoint, nrow=TT, ncol=nofp, byrow=TRUE) #matrix(c(1,2,3), nrow=100, ncol=3, byrow=T)
  
  #IWMDE 
  time_iwmde1 <- proc.time()[1:3]
  M_mean <- apply(mcmc, 2, mean)
  D <- mcmc - matrix(M_mean, nrow=TT, ncol=(nofp+nonfp), byrow=TRUE)
  sigma_all <- t(D)%*%D/(TT-1)
  val <- apply(mcmc, 1, dmvnorm, mean = M_mean, sigma = sigma_all, log = TRUE)
  val <- val - apply(mcmc[,-c(1:nofp)], 1, dmvnorm, mean = M_mean[-c(1:nofp)], sigma = sigma_all[-c(1:nofp),-c(1:nofp)], log = TRUE)
  val <- val + 
         apply(target, 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") -
         apply(mcmc, 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted")

  den_est[k] <- control(val, TT)
  time_iwmde2 <- proc.time()[1:3]
  iwmdetime <- rbind(iwmdetime, time_iwmde2 - time_iwmde1)
  
  
  
  #PWMDE1 
  ncutslice <- 5
  time_pwmde1 <- proc.time()[1:3]
  M_mean <- apply(mcmc, 2, mean)
  D <- mcmc - matrix(M_mean, nrow=TT, ncol=(nofp+nonfp), byrow=TRUE)
  #create the covariance matrix of those parameters of interest in the margianl dist'n
  D <- D[, 1:nofp]
  tar_sigma <- t(D)%*%D/(TT-1)
  #find eigenvalues & eigenvectors
  evaluescov <- diag(eigen(tar_sigma)$values)
  evectorscov <- eigen(tar_sigma)$vectors
  PsqrtD <- evectorscov%*%(evaluescov)^0.5
  inevectorscov_sqrtdiag <- solve(PsqrtD)
  #partjo1 is the Jocobi term
  partjo1 <- det(PsqrtD)
  lpartjo1 <- log(abs(partjo1) )
  scale_vec <- t(PsqrtD%*%evectorscov[,1] )
  #standardized sample of parameters of interested 
  #used to get a distance from the center (which marks which elliptical ring it belongs to)
  newpara_mcmc <- t(inevectorscov_sqrtdiag%*%t(D))
  r_ <- sqrt(apply(newpara_mcmc^2, 1, sum))  
  temp <- which(r_ < rcover)
  interval <- seq(0, rcover, length=(ncutslice+1) )
  reprp <- matrix(interval[-(ncutslice+1)] + rcover/ncutslice/2, ncol=1 )
  #https://en.wikipedia.org/wiki/Volume_of_an_n-ball
  rarea <- pi^(nofp/2)*interval^nofp/gamma(nofp/2+1)  
  rvol <- log(diff(rarea) ) + lpartjo1
  
  position <- ceiling(r_[temp]*ncutslice/rcover) 
  mcmc_p <- cbind(position, mcmc[temp, ])
  rings_ <- apply(mcmc_p, 1, LOR_partition, nslice=ncutslice)
  val <- rings_ + 
         apply(target[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") - 
         apply(mcmc[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") 

  den_est_pwmde1[k] <- control(val, TT)
  time_pwmde2 <- proc.time()[1:3]
  pwmdetime1 <- rbind(pwmdetime1, time_pwmde2 - time_pwmde1)
  
  
  
  #PWMDE2 
  ncutslice <- 10
  time_pwmde1 <- proc.time()[1:3]
  M_mean <- apply(mcmc, 2, mean)
  D <- mcmc - matrix(M_mean, nrow=TT, ncol=(nofp+nonfp), byrow=TRUE)
  #create the covariance matrix of those parameters of interest in the margianl dist'n
  D <- D[, 1:nofp]
  tar_sigma <- t(D)%*%D/(TT-1)
  #find eigenvalues & eigenvectors
  evaluescov <- diag(eigen(tar_sigma)$values)
  evectorscov <- eigen(tar_sigma)$vectors
  PsqrtD <- evectorscov%*%(evaluescov)^0.5
  inevectorscov_sqrtdiag <- solve(PsqrtD)
  #partjo1 is the Jocobi term
  partjo1 <- det(PsqrtD)
  lpartjo1 <- log(abs(partjo1) )
  scale_vec <- t(PsqrtD%*%evectorscov[,1] )
  #standardized sample of parameters of interested 
  #used to get a distance from the center (which marks which elliptical ring it belongs to)
  newpara_mcmc <- t(inevectorscov_sqrtdiag%*%t(D))
  r_ <- sqrt(apply(newpara_mcmc^2, 1, sum))  
  temp <- which(r_ < rcover)
  interval <- seq(0, rcover, length=(ncutslice+1) )
  reprp <- matrix(interval[-(ncutslice+1)] + rcover/ncutslice/2, ncol=1 )
  #https://en.wikipedia.org/wiki/Volume_of_an_n-ball
  rarea <- pi^(nofp/2)*interval^nofp/gamma(nofp/2+1)  
  rvol <- log(diff(rarea) ) + lpartjo1
  
  position <- ceiling(r_[temp]*ncutslice/rcover) 
  mcmc_p <- cbind(position, mcmc[temp, ])
  rings_ <- apply(mcmc_p, 1, LOR_partition, nslice=ncutslice)
  val <- rings_ + 
         apply(target[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") - 
         apply(mcmc[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") 
  
  den_est_pwmde2[k] <- control(val, TT)
  time_pwmde2 <- proc.time()[1:3]
  pwmdetime2 <- rbind(pwmdetime2, time_pwmde2 - time_pwmde1)
  
  
  
  #PWMDE3 
  ncutslice <- 15
  time_pwmde1 <- proc.time()[1:3]
  M_mean <- apply(mcmc, 2, mean)
  D <- mcmc - matrix(M_mean, nrow=TT, ncol=(nofp+nonfp), byrow=TRUE)
  #create the covariance matrix of those parameters of interest in the margianl dist'n
  D <- D[, 1:nofp]
  tar_sigma <- t(D)%*%D/(TT-1)
  #find eigenvalues & eigenvectors
  evaluescov <- diag(eigen(tar_sigma)$values)
  evectorscov <- eigen(tar_sigma)$vectors
  PsqrtD <- evectorscov%*%(evaluescov)^0.5
  inevectorscov_sqrtdiag <- solve(PsqrtD)
  #partjo1 is the Jocobi term
  partjo1 <- det(PsqrtD)
  lpartjo1 <- log(abs(partjo1) )
  scale_vec <- t(PsqrtD%*%evectorscov[,1] )
  #standardized sample of parameters of interested 
  #used to get a distance from the center (which marks which elliptical ring it belongs to)
  newpara_mcmc <- t(inevectorscov_sqrtdiag%*%t(D))
  r_ <- sqrt(apply(newpara_mcmc^2, 1, sum))  
  temp <- which(r_ < rcover)
  interval <- seq(0, rcover, length=(ncutslice+1) )
  reprp <- matrix(interval[-(ncutslice+1)] + rcover/ncutslice/2, ncol=1 )
  #https://en.wikipedia.org/wiki/Volume_of_an_n-ball
  rarea <- pi^(nofp/2)*interval^nofp/gamma(nofp/2+1)  
  rvol <- log(diff(rarea) ) + lpartjo1
  
  position <- ceiling(r_[temp]*ncutslice/rcover) 
  mcmc_p <- cbind(position, mcmc[temp, ])
  rings_ <- apply(mcmc_p, 1, LOR_partition, nslice=ncutslice)
  val <- rings_ + 
         apply(target[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") - 
         apply(mcmc[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") 
  
  den_est_pwmde3[k] <- control(val, TT)
  time_pwmde2 <- proc.time()[1:3]
  pwmdetime3 <- rbind(pwmdetime3, time_pwmde2 - time_pwmde1)
  
  
  
  #PWMDE4 
  ncutslice <- 20
  time_pwmde1 <- proc.time()[1:3]
  M_mean <- apply(mcmc, 2, mean)
  D <- mcmc - matrix(M_mean, nrow=TT, ncol=(nofp+nonfp), byrow=TRUE)
  #create the covariance matrix of those parameters of interest in the margianl dist'n
  D <- D[, 1:nofp]
  tar_sigma <- t(D)%*%D/(TT-1)
  #find eigenvalues & eigenvectors
  evaluescov <- diag(eigen(tar_sigma)$values)
  evectorscov <- eigen(tar_sigma)$vectors
  PsqrtD <- evectorscov%*%(evaluescov)^0.5
  inevectorscov_sqrtdiag <- solve(PsqrtD)
  #partjo1 is the Jocobi term
  partjo1 <- det(PsqrtD)
  lpartjo1 <- log(abs(partjo1) )
  scale_vec <- t(PsqrtD%*%evectorscov[,1] )
  #standardized sample of parameters of interested 
  #used to get a distance from the center (which marks which elliptical ring it belongs to)
  newpara_mcmc <- t(inevectorscov_sqrtdiag%*%t(D))
  r_ <- sqrt(apply(newpara_mcmc^2, 1, sum))  
  temp <- which(r_ < rcover)
  interval <- seq(0, rcover, length=(ncutslice+1) )
  reprp <- matrix(interval[-(ncutslice+1)] + rcover/ncutslice/2, ncol=1 )
  #https://en.wikipedia.org/wiki/Volume_of_an_n-ball
  rarea <- pi^(nofp/2)*interval^nofp/gamma(nofp/2+1)  
  rvol <- log(diff(rarea) ) + lpartjo1
  
  position <- ceiling(r_[temp]*ncutslice/rcover) 
  mcmc_p <- cbind(position, mcmc[temp, ])
  rings_ <- apply(mcmc_p, 1, LOR_partition, nslice=ncutslice)
  val <- rings_ + 
         apply(target[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") - 
         apply(mcmc[temp, ], 1, dmvt, delta = delta_, sigma = sigma_, df = df_, log = TRUE, type = "shifted") 
  
  den_est_pwmde4[k] <- control(val, TT)
  time_pwmde2 <- proc.time()[1:3]
  pwmdetime4 <- rbind(pwmdetime4, time_pwmde2 - time_pwmde1)
}

output <- cbind(den_est, den_est_pwmde1, den_est_pwmde2, den_est_pwmde3, den_est_pwmde4)
apply(output, 2, mean)
write.table(output, file=(paste('output_', DUMZ, '.txt', sep='')), sep=' ', col.names = FALSE, row.names = FALSE)

time <- cbind(iwmdetime[, 3], pwmdetime1[, 3], pwmdetime2[, 3], pwmdetime3[, 3], pwmdetime4[, 3])
write.table(time, file=(paste('time_', DUMZ, '.txt', sep='')), sep=' ', col.names = FALSE, row.names = FALSE)

