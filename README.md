# Section 4
This is the readme file for the scripts "mvt.R" and "evaluation.R", which are used to compute the marginal densities in Section 4.

# mvt.R
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
   
# evaluation.R
"evaluation.R" contains three parts: (I) Required Packages; (II) Settings; and (III) Computation & Results. We list them in details in the following. 

(I) Required Packages:
   1. This code requires a package: mvtnorm.   
   
(II) Settings:   
   1. Set the no. of focused parameters, the no. of non-focused parameters, and the point of interest. 

(III) Computation & Results:
   1. Compute the true marginal posterior density of the point of interest.
   2. Read in all "output_DUMZ.txt" and creat a list containing all estimates (1000 repeats, each with 5 estimates).
   3. Based on the true value and the list, compute means, MCSEs, and RMSEs of the IWMDE and PWMDE.



# Section 5
This is the readme file for the scripts "Main_pgm.R" and "PWMDE_pgm.R", which are used to compute the marginal and joint marginal posterior distributions in Section 5.

# Main_pgm.R
"Main_pgm.R" is the main code for the computation. It contains (I) Required Packages & Data & script; (II) Setting; and (III) Computation. We list them in details in the following. 

(I) Required Packages & Data & script:
   1. The main code requires packages: mvtnorm, coda, MASS, and MCMCpack.
   2. "DID" denotes the DID data set. The size of DID is 94x2 (94 subjects with 2 variables), where the first column is about the DID group a subject belongs to ("DID patients" coded 
      as 1, "mimic normal" coded as 2, "symptom simulated" coded as 3, and "true amnesic" coded as 4.); and the second one records the memory scores of the subjects.   
   3. "mcmc" is an MCMC sample (100,000x5) generated in Fortran. In computation, we actually only use the first 1,000 MCMC sample points to generate Figures 1, 2, 3, and 5.
   4. "PWMDE_pgm.R" is used to provide (1) the unnormalized posterior density evaluated at the assigned point, say an MCMC sample point or a representative point for a partition subset, and (2) the floating control for the summation of the weighted-averaging ratios of the unnormalized posterior densities in the PWMDE estimator (or in the CMDE estimator).   		

(II) Setting:
   1. Set the parameters.
   2. Decide the number of partition subsets for PWMDE.

(III) Computation:
   After setting up (I) and (II), the marginal and joint marginal posterior distributions of \mu1, \mu2, \mu3, \mu4, (\mu1,\mu4), and (\mu1,\mu3) are estimated in order by using CMDE, KDE, and PWMDE, respectively. For the joint marginal posterior distributions of (\mu1,\mu4) and (\mu1,\mu3), we further use Mathematica to generate Figure 3 and 5 based on the estimated values and differences of values (KDE versus CMDE, and PWMDE versus CMDE). Piecewise L_1 distance is also calculated for each case. 

# PWMDE_pgm.R  
"PWMDE_pgm.R" is the script providing (1) updv(): for the calculation of the unnormalized posterior density in log scale, and (2) floating_control(): for the floating control.


  
