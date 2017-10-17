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

## Functions for simulating selection in the lattice Wright-Fisher model
## and generating observations from the model.


## Simulate a wf population on a k1 by k2 grid with a fixed size n in each deme,
## initial allele frequency p0 and selection coefficient s - can be either k1 by k2
## matrices or constants, and migration rate m individuals per deme per neigbour
## IN this simulation, you get born, then you migrate. 

simulate.wright.fisher.lattice <- function(n, g, p0, s, k1, k2, m){
  if(!is.matrix(p0)==1){p0=matrix(p0, nrow=k1, ncol=k2)}
  if(!is.matrix(s)==1){s=matrix(s, nrow=k1, ncol=k2)}

  mrm <- mig.rate.mat(k1, k2)

  counts <- array(0, dim=c(g,k1,k2))
  f <- array(0, dim=c(k1,k2,g))

  f[,,1] <- p0

  ## I'm sorry R, but too stupid to vectorise
  for(t in 2:(g)){
    ## born counts the number of A individuals born in that square
    ## incoming counts the number of A individuals migrating into that square
    born <- matrix(0, nrow=k1, ncol=k2)
    incoming <- matrix(0, nrow=k1, ncol=k2)
    outgoing <- matrix(0, nrow=k1, ncol=k2)
    
    ## first select n inidivduals which are born in that square 
    ## And also the number who migrate out. 
    for(i in 1:k1){
      for(j in 1:k2){
        p.sel <- f[i,j,t-1]*(1+s[i,j])/(1+f[i,j,t-1]*s[i,j])
        born[i,j] <- rbinom(1, n, p.sel)
        outgoing[i,j] <- sum(sample(c(rep(1,born[i,j]),rep(0,n-born[i,j])),mrm[i,j]*m, replace=FALSE))
      }
    }

    ## Now select the migrants - we do N,S,E,W separately so it's easier to count (though it means more loops)
    ## This is a horrible horrible horrible so horrible way to do it.
    probs <- mrm
    temp.outgoing <- outgoing

    if(k1>1){
      for(i in 2:k1){
        for(j in 1:k2){
          n.mig <- sum(sample(c(rep(1,temp.outgoing[i,j]),rep(0,probs[i,j]*m-temp.outgoing[i,j])),m), replace=FALSE)
          temp.outgoing[i,j] <- temp.outgoing[i,j]-n.mig
          incoming[i-1,j] <- incoming[i-1,j]+n.mig
          probs[i,j] <- probs[i,j]-1
        }
      }

      for(i in 1:(k1-1)){
        for(j in 1:k2){
          n.mig <- sum(sample(c(rep(1,temp.outgoing[i,j]),rep(0,probs[i,j]*m-temp.outgoing[i,j])),m), replace=FALSE)
          temp.outgoing[i,j] <- temp.outgoing[i,j]-n.mig
          incoming[i+1,j] <- incoming[i+1,j]+n.mig
          probs[i,j] <- probs[i,j]-1
        }
      }
    }

    if(k2>1){
      for(i in 1:k1){
        for(j in 2:k2){
          n.mig <- sum(sample(c(rep(1,temp.outgoing[i,j]),rep(0,probs[i,j]*m-temp.outgoing[i,j])),m), replace=FALSE)
          temp.outgoing[i,j] <- temp.outgoing[i,j]-n.mig
          incoming[i,j-1] <- incoming[i,j-1]+n.mig
          probs[i,j] <- probs[i,j]-1
        }
      }
    
      for(i in 1:k1){
        for(j in 1:(k2-1)){
          n.mig <- sum(sample(c(rep(1,temp.outgoing[i,j]),rep(0,probs[i,j]*m-temp.outgoing[i,j])),m), replace=FALSE)
          temp.outgoing[i,j] <- temp.outgoing[i,j]-n.mig
          incoming[i,j+1] <- incoming[i,j+1]+n.mig
          probs[i,j] <- probs[i,j]-1
        }
      }
    }
    
    f[,,t] <- (born + incoming - outgoing)/n
  }

  return(f)
}                  

## generate observations from path, with N observations
## at each time point. f and N are k1 x k2 x g arrays where
## k1*k2 is the number of demes and g is the number of generations

generate.lattice.observations.from.path <- function(f, N){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  g <- dim(f)[3]

  N.A <- 0*N

  for(t in 1:g){
    for(i in 1:k1){
      for(j in 1:k2){
        N.A[i,j,t] <- rbinom(1,N[i,j,t],f[i,j,t])
      }
    }
  }
  return(list(N=N, N.A=N.A))
}

## Generate a path, and some observations, return the whole thing

generate.lattice.observations <- function(Ne, g, p0, s, k1, k2, m, missing.p=0.1, size.params=list(N=100, p=0.5)){
  f <- simulate.wright.fisher.lattice(Ne, g, p0, s, k1, k2, m)
  N <- array(0,dim=c(k1,k2,g))
  for(t in 1:g){
    entries <- rbinom(k1*k2, size.params$N, size.params$p)
    entries[as.logical(rbinom(k1*k2, 1, missing.p))] <- 0
    N[,,t] <- entries
  }

  obs <- generate.lattice.observations.from.path(f, N)
  return(list(f=f, obs=obs))
}

