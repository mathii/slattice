########################################################################################################
##    Copyright 2012 Iain Mathieson

##    Licensed under the Apache License, Version 2.0 (the "License");
##    you may not use this file except in compliance with the License.
##    You may obtain a copy of the License at

##        http://www.apache.org/licenses/LICENSE-2.0

##    Unless required by applicable law or agreed to in writing, software
##    distributed under the License is distributed on an "AS IS" BASIS,
##    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##    See the License for the specific language governing permissions and
##    limitations under the License.
########################################################################################################

## This script runs some simple examples of the estimator described in Mathieson & Mcvean 2013.

## load all the required libraries
source("include.R")
set.seed(11235)

## 1) Single population case. 
cat(paste(paste(rep("*", 20), collapse=""), " Starting single population example ", paste(rep("*", 20), collapse=""), "\n", sep=""))

## Parameters to simulate a 1-d example
Ne <- 1000                              #N_e
g<-100                                  #Number of generations
p0 <- 0.1                               #Initial freq
s <- 0.05                               #Selection coefficient

## Simulate. data$obs is the data required - it has two columns named N and N.A giving the total
## number of samples and the number of the selected type at each generation. i.e. n_t and a_t
## in the paper. 
data<-generate.observations(Ne, g, p0, s, missing.p=0.9, size.params=list(N=100,p=0.5))

dev.new()
par(mfrow=c(2,2))
## Run EM estimator. Set method to different values for different estimators. "Soft EM" is the EM algorithm
## described in the paper. You should see a plot showing the observations and the posterior for f_t (set verbose=FALSE to avoid)
## the resulting estimate contains the estimated value of s, and the log-likelihood from the final model
estimate <- estimate.s(data$obs, Ne, method="Soft EM", verbose=TRUE)
cat(paste("Final s =", estimate$s, "\n\n"))

## 2) lattice case

cat(paste(paste(rep("*", 20), collapse=""), " Starting lattice population example ", paste(rep("*", 20), collapse=""), "\n", sep="" ))

## Change these to try out different grid shapes. 
k1 <- 5                                  #Number of rows of demes
k2 <- 3                                 #Number of cols of demes

## Parameters for simulatuon
Ne <- 1000                               #N_e in each deme
g<-100                                  #Number of generations
p0 <- 0.1                               #Initial frequency
s <- matrix(0.06*seq(1,-1,length.out=k1), k1, k2, byrow=FALSE) #S^{ij} - matrix of selection coefficients
m <- 0.04                                #Scaled migration rate

## Simulate - note we need to give the absolute migration rate, only for the simulations. 
lattice.data<-generate.lattice.observations(Ne, g, p0, s, k1, k2, Ne*m, missing.p=0.9, size.params=list(N=100, p=0.5))

## Run EM estimator. Again, method selects differnt estimators.
## if you specify M to be non-null, that will be used as the fixed value
## You do not need to specify initial.M, but the estimator will perform better if you are close.
lattice.estimate<-estimate.s.m(lattice.data$obs, Ne, M=NULL, update="Soft EM", max.iters=10, verbose=TRUE, initial.M=m)

cat(paste("Final M =", lattice.estimate$M, "\n"))
cat("Final s =\n")
print(lattice.estimate$s)

## plot the results - true frequency(green), observations (blue), viterbi path (orange), and estimated selection coefficients (background colours).
dev.new()
plot.wright.fisher.lattice.observations(lattice.data$obs, lattice.data$f, lattice.estimate$f, est.s=lattice.estimate$s, error.bars=TRUE, main="Lattice Example")
