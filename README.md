# Partition-Weighted-Marginal-Density-Estimator

This is the readme file for the scripts "mvt.R" and "evaluation.R", which are used to compute the marginal densities in Section 4.
#########################################################################################################################################
"mvt.R" contains four parts: (I) Required Packages; (II) Functions (III) Settings; and (IV) Computation & Results. We list them in details in the following. 

(I) Required Packages:
   1. This code requires a package: mvtnorm.   

(II) Functions:
   1. control().
   2. LOR_partition().
   
(III) Settings:
   1. Set the random seed, DUMZ, equal to 1~100 ("seeds.txt").
   2. Set the parameters (ex. no. of focused parameters, no. of non-focused parameters, the conditional working parameter space (that is, radius)). 

(IV) Computation & Results:
   1. 10 IWMDE estimates (because each random seed will execute 10 simulation repeats) are stored in "den_est".
   2. 10 PWMDE estimates with K=5 are stored in "den_est_pwmde1".
   3. 10 PWMDE estimates with K=10 are stored in "den_est_pwmde2".  
   4. 10 PWMDE estimates with K=15 are stored in "den_est_pwmde3".
   5. 10 PWMDE estimates with K=20 are stored in "den_est_pwmde4".   
   6. Results of all estimates (10 repeats, each with 5 estimates) are saved as "output_DUMZ.txt".
   
 
#########################################################################################################################################
"evaluation.R" contains three parts: (I) Required Packages; (II) Settings; and (III) Computation & Results. We list them in details in the following. 

(I) Required Packages:
   1. This code requires a package: mvtnorm.   
   
(II) Settings:   
   1. Set the no. of focused parameters, the no. of non-focused parameters, and the point of interest. 

(III) Computation & Results:
   1. Compute the true marginal posterior density of the point of interest.
   2. Read in all "output_DUMZ.txt" and creat a list containing all estimates (1000 repeats, each with 5 estimates).
   3. Based on the true value and the list, compute means, MCSEs, and RMSEs of the IWMDE and PWMDE.
   
