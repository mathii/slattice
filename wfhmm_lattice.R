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

## Just like wfhmm.R, but for the lattice model. We assume that the
## frequency increments from t to t+1 are independent, conditional
## on the frequency history up to t.

source("wfhmm.R")
source("lattice.R")

## Discretise and call the lattice viterbi algorithm

wfhmm.lattice.call <- function(obs, s, M, Ne, estimated.f, h=0.5, viterbi=TRUE, paths=0, grid=100, extend=0.1, likelihood=TRUE, forward.backward=TRUE){
  obs.f <- obs$N.A/obs$N
  k1 <- dim(obs.f)[1]
  k2 <- dim(obs.f)[2]
  g <- dim(obs.f)[3]
  f.max <- min(1,extend+max(obs.f, na.rm=TRUE))
  f.min <- max(0,min(obs.f-extend, na.rm=TRUE))
  disc.f <- seq(f.min, f.max, length.out=grid)
  interval <- (f.max-f.min)/(grid-1)
  params <- list(s=s, M=M, Ne=Ne, interval=interval, obs=obs, est.f=estimated.f, h=h, states=disc.f)
  
  results <- list()
  if(viterbi){
    tb <- wfhmm.lattice.viterbi(disc.f, obs.f, wfhmm.lattice.emission, wfhmm.lattice.transmission, params)
    f.path <- array(disc.f[tb], dim=c(k1,k2,g))
    results$viterbi <- f.path
  }
  if(paths>0){
    all.paths <- array(0, dim=c(k1,k2,g,paths))
    fwd <- wfhmm.lattice.forward(disc.f, obs.f, wfhmm.lattice.emission, wfhmm.lattice.transmission, params)
    for( i in 1:paths){
      all.paths[,,,i] <- wfhmm.lattice.simulate(disc.f, fwd$forward.matrix, wfhmm.lattice.emission, wfhmm.lattice.transmission, params)
    }
    results$paths <- all.paths
  }
  if(forward.backward|likelihood){
    results$fb <- wfhmm.lattice.forward.backward(disc.f, obs.f, wfhmm.lattice.emission, wfhmm.lattice.transmission, params)
  }  
  if(likelihood){
    results$log.likelihood <- sum(log(results$fb$scaling.factor))
  }

  results$params <- params
  results$emission.func=wfhmm.lattice.emission
  results$transition.func=wfhmm.lattice.transmission
  
  return(results)
}

## Copy parameters and get the right obs
## Helper function for the loops in the hmm algorithm

wfhmm.lattice.copy.params <- function(params, i, j){
  params$i <- i
  params$j <- j
  params$obs$N <- params$obs$N[i,j,]
  params$obs$N.A <- params$obs$N.A[i,j,]
  return(params)
}

## lattice Viterbi tb<-wfhmm.lattice.viterbi(disc.f, obs.f, wfhmm.lattice.emission, wfhmm.lattice.transmission, params)
## Todo, just call the wfhmmm lattice function.

wfhmm.lattice.viterbi <- function(states, observations, emission.func, transition.func, params=list()){
  n.states <- length(states)
  k1 <- dim(observations)[1]
  k2 <- dim(observations)[2]
  n.obs <- dim(observations)[3]

  vit.mat <- array(0, dim=c(k1,k2,n.states, n.obs))
  tb <- array(0,dim=c(k1,k2, n.obs))

  for(i in 1:k1){
    for(j in 1:k2){
      these.params <- wfhmm.lattice.copy.params(params, i, j)
      tb[i,j,] <- wfhmm.viterbi(states, observations[i,j,], emission.func, transition.func, these.params)
    
    }
  }

  return(tb)
}

## Forward algorithm for the lattice hmm; This, and the backwards
## algorithm below are just wrappers for the 1-d functions
wfhmm.lattice.forward <- function(states, observations, emission.func, transition.func, params=list()){
  n.states <- length(states)
  k1 <- dim(observations)[1]
  k2 <- dim(observations)[2]

  n.obs <- dim(observations)[3]

  f.mat <- array(0, dim=c(k1,k2,n.states, n.obs))
  scale.factors.mat <-  array(0, dim=c(k1,k2, n.obs))

  for(i in 1:k1){
    for(j in 1:k2){
      these.params <- wfhmm.lattice.copy.params(params, i, j)
      this.fwd <- wfhmm.forward(states, observations[i,j,], emission.func, transition.func, these.params)
      f.mat[i,j,,] <- this.fwd$forward.matrix
      scale.factors.mat[i,j,] <- this.fwd$scaling.factor
    }
  }

  return(list(forward.matrix=f.mat, scaling.factor=scale.factors.mat))
}

## backward algorithm
wfhmm.lattice.backwards <- function(states, observations, emission.func, transition.func, params=list(), scaling.factors=NULL){
  n.states <- length(states)
  k1 <- dim(observations)[1]
  k2 <- dim(observations)[2]
  n.obs <- dim(observations)[3]

  b.mat <- array(0, dim=c(k1,k2,n.states, n.obs))
  if(all(is.null(scaling.factors))){scaling.factors <- array(1, dim=c(k1,k2, n.obs))}

  for(i in 1:k1){
    for(j in 1:k2){
      these.params <- wfhmm.lattice.copy.params(params, i, j)
      these.scaling.factors <- scaling.factors[i,j,]
      this.bwd <- wfhmm.backwards(states, observations[i,j,], emission.func, transition.func, these.params, these.scaling.factors)
      b.mat[i,j,,] <- this.bwd$backwards.matrix
    }
  }

  return(list(backwards.matrix=b.mat))
}

## well...
wfhmm.lattice.forward.backward <- function(states, observations, emission.func, transition.func, params){
  fwd <- wfhmm.lattice.forward(states, observations, emission.func, transition.func, params)
  bwd <- wfhmm.lattice.backwards( states, observations, emission.func, transition.func, params, scaling.factors=fwd$scaling.factor)
    fwd.bwd <- fwd$forward.matrix*bwd$backwards.matrix
  return(list(f.mat=fwd$forward.matrix, b.mat=bwd$backwards.matrix, fb.mat=fwd.bwd, scaling.factor=fwd$scaling.factor))
}

## Simulate a path from the posterior distribution of paths
wfhmm.lattice.simulate <- function( states, fwd, emission.func, transition.func, params){
  k1 <- dim(fwd)[1]
  k2 <- dim(fwd)[2]
  n.obs <- dim(fwd)[4]
  if(length(states)!=dim(fwd)[3]){stop("lenth of states does not match dimensions of forward matrix")}
    
  path <- array(0, dim=c(k1,k2, n.obs))

  for(i in 1:k1){
    for(j in 1:k2){
      these.params <- wfhmm.lattice.copy.params(params, i, j)
      path[i,j,] <- wfhmm.simulate(states, fwd[i,j,,], emission.func, transition.func, these.params)
    }
  }
  return(path)
}

## binomial emission func, conditional on frequency f
## of seeing N.A out of N observations 

wfhmm.lattice.emission <- function( N.A, N, f){
  return(dbinom(N.A, N, f))
}

## transition func, taking into account migration.

wfhmm.lattice.transmission <- function(f.from, f.to, t, params){
  i <- params$i
  j <- params$j
  k1 <- dim(params$est.f)[1]
  k2 <- dim(params$est.f)[2]
  f.previous <- params$est.f[,,t]
  if(is.null(dim(f.previous))){dim(f.previous) <- c(k1,k2)} #for k1/k2=1 case
  mrm <- mig.rate.mat(k1,k2)
  h <- params$h

  mu <- (1-mrm[i,j]*params$M)*f.from+2*params$s[i,j]*f.from*(1-f.from)*(f.from+h*(1-2*f.from))
  sig <- f.from*(1-f.from)/params$Ne

  if(i<k1){
    f.add <- f.previous[i+1,j]
    if(is.na(f.add)){f.add <- f.from+2*params$s[i,j]}
    mu <- mu+params$M*f.add
  }
  if(i>1){
    f.add <- f.previous[i-1,j]
    if(is.na(f.add)){f.add <- f.from+2*params$s[i,j]}
    mu <- mu+params$M*f.add
  }
  if(j<k2){
    f.add <- f.previous[i,j+1]
    if(is.na(f.add)){f.add <- f.from+2*params$s[i,j]}
    mu <- mu+params$M*f.add
  }
  if(j>1){
    f.add <- f.previous[i,j-1]
    if(is.na(f.add)){f.add <- f.from+2*params$s[i,j]}
    mu <- mu+params$M*f.add
  }
 
  sig <- pmax(1e-8,sqrt(sig))

  int <- params$interval/2
  return(pnorm(f.to+int,mu,sig)-pnorm(f.to-int,mu,sig))
}

## path likelihood
wfhmm.model.path.lattice.likelihood <- function(f, transition.func, params){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  ll <- 0
  for(i in 1:k1){
    for(j in 1:k2){
      this.params <- params
      this.params$i <- i
      this.params$j <- j
      ll <- ll+wfhmm.model.path.log.likelihood(f[i,j,], transition.func, this.params)
    }
  }
  return(ll)
}

## observation likelihood
wfhmm.observation.lattice.likelihood <- function(f, obs){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  ll <- 0
  for(i in 1:k1){
    for(j in 1:k2){
      this.obs <- list(N=obs$N[i,j,], N.A=obs$N.A[i,j,]) 
      ll <- ll+wfhmm.observation.log.likelihood(f[i,j,], this.obs)
    }
  }
  return(ll)
}

