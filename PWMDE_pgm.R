#unnormalizing posterior density value here ignores some constant parts 
#due to the cancellation of those terms in the PWMDE formular (ratio)
updv <- function(xx1, xx2, xx3, xx4, xx5){
  #H1a2: mu2 > {mu4, mu1} > mu3
  DID_ <- DID[which(DID[,1]==1),2]
  Control <- DID[which(DID[,1]==2),2]
  Control_sim <- DID[which(DID[,1]==3),2]
  Control_amne <- DID[which(DID[,1]==4),2]
  part1 <- sum((DID_-xx1)^2)+sum((Control-xx2)^2)+sum((Control_sim-xx3)^2)+sum((Control_amne-xx4)^2)
  part1 <- -part1/(2*xx5^2) - 47*log(2*pi*xx5^2)
  part2 <- dnorm(xx1, mean = 0, sd=(xx5/sqrt(0.19)), log=TRUE)
  part2 <- part2 + dnorm(xx2, mean = 0, sd=(xx5/0.5), log=TRUE)  
  part2 <- part2 + dnorm(xx3, mean = 0, sd=(xx5/0.5), log=TRUE)  
  part2 <- part2 + dnorm(xx4, mean = 0, sd=(xx5/0.5), log=TRUE)  
  part2 <- part2 + dgamma(xx5^2, shape=0.0001, rate=0.0001, log=TRUE)
  updv_ <- part1 + part2
  return(updv_)
} 

#purpose: numerical control when summing up exp() terms
#input: "ratios" includes all MCMC sample ratios, which have been assigned to weighted averages
#input: "tot" is the size of an MCMC sample
#output is the marginal posterior density estimate in log scale
floating_control <- function(ratios,  tot)
{
  b <- ratios[!is.na(ratios)]
  print(length(b))
  est_d_r <- 0
  b.max <- max(b)
  est_d_r <- log(sum(exp(b-b.max)))+b.max - log(tot) 
  return(est_d_r)
}


