################################################################################
##    Copyright 2014 Iain Mathieson

##    Licensed under the Apache License, Version 2.0 (the "License");
##    you may not use this file except in compliance with the License.
##    You may obtain a copy of the License at

##        http://www.apache.org/licenses/LICENSE-2.0

##    Unless required by applicable law or agreed to in writing, software
##    distributed under the License is distributed on an "AS IS" BASIS,
##    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##    See the License for the specific language governing permissions and
##    limitations under the License.
################################################################################

## Profile some simple examples of the estimator.

source("include.R")
set.seed(11235)

## Profile 1-d example

Ne <- 1000                              #N_e
g<-100                                  #Number of generations
p0 <- 0.1                               #Initial freq
s <- 0.02                               #Selection coefficient
data<-generate.observations(Ne, g, p0, s, missing.p=0.9, size.params=list(N=100,p=0.5))

Rprof("profile.out")
t <- system.time(estimate <- estimate.s(data$obs, Ne, method="Soft EM", verbose=FALSE, tol=1e-4))
Rprof(NULL)

cat("Single population\n")
cat("Expected: s = 0.0129619540782165\n\n")
cat(paste("Final s =", estimate$s, "\n"))
cat(paste("Time  t =", t[1], "s\n\n"))

summaryRprof("profile.out")

## profile 2-d example

Ne <- 1000                               #N_e in each deme
g<-100                                  #Number of generations
p0 <- 0.1                               #Initial frequency
s <- matrix(c(0.06, 0.02, -0.02, -0.06), 4, 4, byrow=FALSE) #S^{ij} - matrix of selection coefficients
k <- 4                                                      #square root of number of demes
m <- 0.04                                                   #Scaled migration rate

## Simulate - note we need to give the absolute migration rate, only for the simulations. 
lattice.data<-generate.lattice.observations(Ne, g, p0, s, k, Ne*m, missing.p=0.9, size.params=list(N=100, p=0.5))

Rprof("profile_lattice.out")
t <- system.time(lattice.estimate<-estimate.s.m(lattice.data$obs, Ne, M=NULL, update="Soft EM", max.iters=10, initial.M=m))
Rprof(NULL)

cat("Lattice population\n")
cat("Expected: M = 0.0146838924164589\n 
s =
             [,1]         [,2]        [,3]        [,4]
[1,]  0.058661691  0.039383786  0.07276261  0.04639734
[2,]  0.036322768  0.025556467  0.02950090  0.01355296
[3,] -0.008798924 -0.006494724 -0.04084026 -0.01429270
[4,] -0.133296095 -0.065719422  0.04420272 -0.09474492
\n\n")
cat(paste("Final M =", lattice.estimate$M, "\n"))
cat("Final s =\n")
print(lattice.estimate$s)
cat(paste("Time  t =", t[1], "s\n\n"))

summaryRprof("profile_lattice.out")


