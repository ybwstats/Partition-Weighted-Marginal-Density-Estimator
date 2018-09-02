###Required packages
library(mvtnorm)

######################################################################
###Settings
setwd("/.../")
nofp <- 3         #no. of focused parameters
nonfp <- 2        #no. of non-focused parameters
tpoint <- rep(0, nofp)      #a point of interest

######################################################################
###Computation & Results
#calculate the true marginal posterior density of an interested point 
delta_ <- rep(0, nofp+nonfp)
sigma_ <- NULL
for (i in 1:(nofp+nonfp) ){
  sigma_ <- rbind(sigma_, abs((1-i):(nofp+nonfp-i)) )
}  
sigma_ <- 3*0.5^sigma_
df_ <- 3
True <- dmvt(tpoint, delta = delta_[1:nofp], sigma = sigma_[1:nofp,1:nofp], df = df_, log = TRUE,
             type = "shifted")

#create lists which contain all simulation results (estimates and time)
output_ <- NULL
time_ <- NULL
for (i in 1:100){
  output_ <- rbind(output_, read.table(paste("output_", i, ".txt", sep = "")) )
  time_ <- rbind(time_, read.table(paste("time_", i, ".txt", sep = "")) )
}

#calculate means, MCSEs, and RMSEs of the IWMDE and PWMDE
result_mean <- round(apply(output_, 2, mean) ,3)
result_sd <- round(apply(output_, 2, sd), 3)
result_mse <- round(sqrt(apply((output_-True)^2, 2, mean) ), 3)
aaa <- cbind(result_mean, result_sd, result_mse, apply(time_, 2, mean) )
write.table(aaa, file='aaa.csv', sep=',', col.names = FALSE, row.names = FALSE)

