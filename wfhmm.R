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

## Solving the HMM for the best frequency path, given observed data counts.
## We discretise the underlying frequency space, and use the binomial emission
## probabilities.

## discrectieze the data and get the Viterbi path

wfhmm.call <- function(obs, s, Ne, h=0.5, grid=100, extend=0.1, viterbi=TRUE, forward.backward=FALSE, paths=0, transitions="normal", likelihood="Model"){
  transition.func <- switch(transitions, normal=wfhmm.transition, binomial=wfhmm.transition.binomial, poisson=wfhmm.transition.poisson, NULL)
  if(is.null(transition.func)){stop("Unknown transition func")}

  obs.f <- obs$N.A/obs$N
  f.max <- min(1,extend+max(obs.f, na.rm=TRUE))
  f.min <- max(0,min(obs.f-extend, na.rm=TRUE))

  if(length(grid)==1){                  #If a single number, use as the length of the grid, otherwise, use the actual points
    disc.f <- seq(f.min, f.max, length.out=grid)
    interval <- (f.max-f.min)/(grid-1)
  }
  else{
    disc.f <- grid
    interval <- grid[2]-grid[1]         #Assuming you are evenly spaced. 
  }
  
  params <- list(Ne=Ne, s=s, obs=obs, interval=interval, states=disc.f, h=h)
  
  results=list(params=params)
  results$transition.func <- transition.func
  results$emission.func <- wfhmm.emission
  
  if(viterbi){
    path <- wfhmm.viterbi(disc.f, obs.f, wfhmm.emission, transition.func, params )
    f.viterbi.path <- disc.f[path]
    results$viterbi=f.viterbi.path
  }
  if(paths>0){
    all.paths <- matrix(0, nrow=paths, ncol=length(obs.f))
    fwd<-wfhmm.forward(disc.f, obs.f, wfhmm.emission, transition.func, params )
    for(i in 1:paths){
      all.paths[i,] <- wfhmm.simulate( disc.f, fwd$forward.matrix, wfhmm.emission, transition.func, params)
    }
    results$paths<-all.paths
  }
  if(forward.backward|(likelihood=="Model")){
    results$fb <- wfhmm.forward.backwards(disc.f, obs.f, wfhmm.emission, transition.func, params )
  }
  if(likelihood[1]=="Model"){
    results$log.likelihood <- sum(log(results$fb$scaling.factor))
  }
  
  return(results)
}

## Find the Viterbi path through discretised frequency space.

wfhmm.viterbi <- function(states, observations, emission.func, transition.func, params=list()){
  n.states <- length(states)
  n.obs <- length(observations)

  vit.mat <- matrix(0, nrow=n.states, ncol=n.obs)
  tb.mat <- matrix(0,nrow=n.states, ncol=n.obs)

  ## uniform prior on initial frequency
  vit.mat[,1] <- 1/n.states * emission.func(params$obs$N.A[1], params$obs$N[1], states)
  if(max(vit.mat,1)==0){vit.mat[,1] <- 1/n.states}
  
  for(j in 2:n.obs){
    for(i in 1:n.states){
      ## In general replace this with the generic viterbi calculation.
      this.args <- vit.mat[,j-1] * transition.func(states, states[i], j, params)
      if(any(is.na(this.args))){
        print(j)
        print(states[i])
        print("vit")
        print(vit.mat[,j-1])
        print(transition.func(states, states[i], j, params))
      }
      if(all(this.args==0)){
        best.from <- i
      }else{
        best.from <- which.max(this.args)
      }
      vit.mat[i,j] <- emission.func(params$obs$N.A[j], params$obs$N[j], states[i]) * this.args[best.from]
      tb.mat[i,j] <- best.from
    }
    ## normalise viterbi matrix so the largest entry is 1
    if(max(vit.mat[,j])>0){
      vit.mat[,j] <- vit.mat[,j]/max(vit.mat[,j])
    } else{
      vit.mat[,j] <-  1/n.states
    }
    
  }

  tb <- rep(0,n.obs)
  tb[n.obs] <- which.max(vit.mat[,n.obs])
  for(t in rev(1:(n.obs-1))){
    tb[t] <- tb.mat[tb[t+1],t+1]
  }
  return(tb)
}

## Forward algorithm for the hmm

wfhmm.forward <- function(states, observations, emission.func, transition.func, params=list()){
  n.states <- length(states)
  n.obs <- length(observations)

  N.A <- params$obs$N.A
  N <- params$obs$N
  
  f.mat <- matrix(0, nrow=n.states, ncol=n.obs)
  scaling.factor <- rep(0,n.obs)
  f.mat[,1] <- 1/n.states

  for(i in 1:n.states){
    f.mat[i,1] <- emission.func(N.A[1], N[1], states[i])*f.mat[i,1]
  }
  scaling.factor[1] <- sum(f.mat[,1])
  f.mat[,1] <- f.mat[,1]/scaling.factor[1]

  for(j in 2:n.obs){
    for(i in 1:n.states){
      this.args <- f.mat[,j-1] * transition.func(states, states[i], j-1, params)
      this.elem <- emission.func(N.A[j], N[j], states[i])*sum(this.args)
      f.mat[i,j] <- this.elem
    }
    scaling.factor[j] <- sum(f.mat[,j])
    f.mat[,j] <-   f.mat[,j]/scaling.factor[j]
  }

  return(list(scaling.factor=scaling.factor, forward.matrix=f.mat, P=sum(f.mat[,n.obs])))
}

## backwards algorithm, including scaling factors

wfhmm.backwards <- function(states, observations, emission.func, transition.func, params, scaling.factor=NULL){
  n.obs <- length(observations)
  n.states <- length(states)
  if(all(is.null(scaling.factor))){scaling.factor <- rep(1,n.obs)}

  N.A <- params$obs$N.A
  N <- params$obs$N
  
  b.mat <- matrix(0, nrow=n.states, ncol=n.obs)
  b.mat[,n.obs] <- 1
  for(j in rev(1:(n.obs-1))){
    this.em <- emission.func(N.A[j+1], N[j+1], states)
    for(i in 1:n.states){
      this.args <- b.mat[,j+1] * transition.func(states[i], states, j, params) * this.em
      b.mat[i,j] <- sum(this.args)/scaling.factor[j+1]
    }
  }

  return(list(backwards.matrix=b.mat))
}

## Well, obviously.

wfhmm.forward.backwards <- function(states, observations, emission.func, transition.func, params){
  fwd <- wfhmm.forward(states, observations, emission.func, transition.func, params)
  bwd <- wfhmm.backwards( states, observations, emission.func, transition.func, params, scaling.factor=fwd$scaling.factor)
  fwd.bwd <- fwd$forward.matrix*bwd$backwards.matrix/fwd$P
  return(list(f.mat=fwd$forward.matrix, b.mat=bwd$backwards.matrix, fb.mat=fwd.bwd, scaling.factor=fwd$scaling.factor, P=fwd$P))
}

## simulate a path, given the forward matrix, and emission and transmission functions.
wfhmm.simulate <- function( states, fwd, emission.func, transition.func, params){
  if(!all(dim(fwd)>1)){stop("Not all dimensions of forward matrix are >1")}
  if(NROW(fwd)!=length(states)){stop("Size of forward matrix does not match states")}

  n.states <- NROW(fwd)
  n.obs <- NCOL(fwd)
  
  path <- rep(0, n.obs)
  last.state.idx <- sample(n.states, 1, prob=fwd[,n.obs]) 
  path[n.obs] <- states[last.state.idx]

  for( i in rev(1:(n.obs-1))){
    top <- fwd[,i]*transition.func(states, path[i+1], i, params)
    for(l in 1:n.states){
      top[l] <-top[l]*emission.func(params$obs$N.A[i+1], params$obs$N[i+1], states[l])
    }
    probs <- top/fwd[last.state.idx,i+1]
    last.state.idx <- sample(n.states, 1, prob=probs) 
    path[i] <- states[last.state.idx]

  }
  return(path)
}
## Binomial emission func, conditional on frequency f
## of seeing N.A out of N observations 

wfhmm.emission <- function( N.A, N, f){
    return(dbinom(N.A, N, f))
}

## Normal approximation to the Binomial transition func, conditional on s and Ne
## of going from f.from to somewhere in the range (f.to.lower,f.to.upper)

wfhmm.transition <- function( f.from, f.to, t, params){
  s <- params$s
  Ne <- params$Ne
  int <- params$interval/2
  h <- params$h
  mu <- f.from+2*s*f.from*(1-f.from)*(f.from+h*(1-2*f.from))
  si <- sqrt(f.from*(1-f.from)/Ne)
  return(pnorm(f.to+int,mu,si)-pnorm(f.to-int,mu,si))
}

## Poission transition func

wfhmm.transition.poisson <- function(f.from, f.to, t, params){
  s <- params$s
  Ne <- params$Ne
  h <- params$h
  lambda <- Ne* f.from+2*s*f.from*(1-f.from)*(f.from+h*(1-2*f.from))
  int <- params$interval/2
  return(ppois(Ne*(f.to+int),lambda)-ppois(Ne*(f.to-int), lambda))
}

## Binomial (i.e. exact transition func)

wfhmm.transition.binomial <- function(f.from, f.to, t, params){
  s <- params$s
  Ne <- params$Ne
  int <- params$interval/2
  if(params$h!=0.5){stop("exact transition only supported for h=0.5")}
  p <- pmin(pmax(0,f.from*(1+s)/(1+f.from*s)),1) #can't go outside [0,1]...
  return(pbinom(Ne*(f.to+int), Ne, p) - pbinom(Ne*(f.to-int), Ne, p))
}

## given a selection coefficient and Ne, return the joint likelihood of the
## path and observations. 

wfhmm.log.likelihood <- function(f, s, Ne, obs){
  return(wfhmm.path.log.likelihood(f, Ne, s)+wfhmm.observation.log.likelihood(f, obs))
}

## given a selection coefficient and Ne, return the joint likelihood of the
## path and observations. 

wfhmm.model.log.likelihood <- function(f, obs, params){
  return(wfhmm.model.path.log.likelihood(f, params)+wfhmm.observation.log.likelihood(f, obs))
}

## Compute binomial log-likelihood of path - just the path, given the
## selection coefficient, ignoring any observations.

wfhmm.path.log.likelihood <- function(f, Ne, s){
  N.A <- round(f*Ne)
  N <- rep(Ne, length(N.A))

  g <- length(N.A)

  f.old <- N.A[1:(g-1)]/Ne
  ps <- f.from*(1+s)/(1+f.from*s)
  
  return(sum(dbinom(N.A[2:g],Ne,prob=ps, log=TRUE)))
}

## Compute log-likelihood of path using the normal approximation- just the path, given the
## selection coefficient, ignoring any observations.

wfhmm.model.path.log.likelihood <- function(f, transition.func, params){
  g <- length(f)
  ll <- 0
  for(i in 2:g){
    ll <- ll+log(transition.func(f[i-1],f[i],i-1,params))
  }
  
  return(ll)
}

## compute the log-likelihood of the observations.

wfhmm.observation.log.likelihood <- function(f, obs){
  return(sum(dbinom(obs$N.A, obs$N, f, log=TRUE)))
}

## The expected value of f_{t+1}, given f_t=i

wfhmm.conditional.expectation <- function(call, i, t){
  fwd <- call$fb$f.mat
  bwd <- call$fb$b.mat
  fb <- call$fb$fb.mat

  params <- call$params
  tf <- call$transition.func
  ef <- call$emission.func
  states <- params$states

  probs <- fwd[i,t]*tf(states[i], states, t, params)*bwd[,t+1]/fb[i,t]
  for(j in 1:length(states)){
    probs[j]=probs[j]*ef(params$obs$N.A[t+1], params$obs$N[t+1], states[j])
  }
  probs <- probs/sum(probs)             
  return(sum(probs*states))
}
