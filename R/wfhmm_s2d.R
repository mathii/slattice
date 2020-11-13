########################################################################################################
##    Copyright 2020 Iain Mathieson

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

## Wright-Fisher hmm, allowing selection coefficent to vary between n different values
## s will be estimated using the em algorithm. Basically we encode the n selection
## states by multiplying the number of frequency states by n. 
## Here s is a vector of length n

wfhmm_s2d.call <- function(obs, s, Ne, h=0.5, s.trans.p=1/NROW(obs), grid=100, extend=0.1, viterbi=TRUE, forward.backward=FALSE, paths=0, transitions="normal", likelihood="Model"){
  transition.func <- switch(transitions, normal=wfhmm_s2d.transition,  NULL)
  if(is.null(transition.func)){stop("Unknown transition func")}
  obs.f <- obs$N.A/obs$N
  if(!(length(Ne)==1|length(Ne)==length(obs.f))){
      stop("Ne must be a constant or the same length as the observations")
  }
  params <- wfhmm.setup.params(obs, grid, Ne, s, h, extend)
  params$s.trans.p <- s.trans.p
  
  results=list(params=params)
  results$transition.func <- transition.func
  results$emission.func <- wfhmm.emission

  if(viterbi){
      path <- wfhmm_s2d.viterbi(params$states, obs.f, wfhmm.emission, transition.func, params )
      f.viterbi.path <- rep(params$states, times=length(s))[path]
      s.viterbi.path <- rep(1:length(s), each=length(params$states))[path]
      results$viterbi=f.viterbi.path
      results$s.viterbi=s.viterbi.path
  }
  if(forward.backward|(likelihood=="Model")){
      results$fb <- wfhmm_s2d.forward.backwards(params$states, obs.f, wfhmm.emission, transition.func, params )
      #collapse forward-backward matrix
      N.s.states <- length(params$s) 

      i=1
      start <- grid*(i-1)+1
      end <- grid*i
      results$fb$f.fb <- results$fb$fb[start:end,]
      results$fb$s.fb <- matrix(0, nrow=N.s.states, ncol=NCOL(results$fb$fb))
      results$fb$s.fb[1,] <- colSums(results$fb$fb[start:end,])
      if(N.s.states>1){
          for(i in 1:N.s.states){
              start <- grid*(i-1)+1
              end <- grid*i
              results$fb$f.fb <- results$fb$f.fb + results$fb$fb[start:end,]
              results$fb$s.fb[i,] <- colSums(results$fb$fb[start:end,])
          }
      }
  }
  if(likelihood=="Model"){
      results$log.likelihood <- sum(log(results$fb$scaling.factor))
  }
  
  return(results)
}

wfhmm_s2d.viterbi <- function(states, observations, emission.func, transition.func, params=list()){
    #states is just the frequency states here
    N.s.states <- length(params$s) 
    n.states <- length(states)
    n.obs <- length(observations)

    vit.mat <- matrix(0, nrow=n.states*N.s.states, ncol=n.obs)
    tb.mat <- matrix(0,nrow=n.states*N.s.states, ncol=n.obs)

    ## uniform prior on initial frequency
    vit.mat[,1] <- (1/n.states/N.s.states) * rep(emission.func(params$obs$N.A[1], params$obs$N[1], states), times=N.s.states)
    if(max(vit.mat,1)==0){vit.mat[,1] <- 1/n.states}
  
    for(j in 2:n.obs){
        for(i in 1:(n.states*N.s.states)){
            ## Are we looking at the same selection coefficient?
            same.s.state <- rep(0,n.states*N.s.states)

            block <- floor((i-1)/params$grid)
            same.s.state[(params$grid*block+1):(params$grid*block+params$grid)] <- 1
            this.args <- vit.mat[,j-1] * transition.func(rep(states, times=N.s.states), states[(i-1)%%n.states+1], j-1, params, same.s.state)
            if(all(this.args==0)){
                best.from <- i
            }else{
                best.from <- which.max(this.args)
            }
            vit.mat[i,j] <- emission.func(params$obs$N.A[j], params$obs$N[j], states[(i-1)%%n.states+1]) * this.args[best.from]
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

wfhmm_s2d.forward <- function(states, observations, emission.func, transition.func, params=list()){
    N.s.states <- length(params$s) 
    n.states <- length(states)
    n.obs <- length(observations)

    N.A <- params$obs$N.A
    N <- params$obs$N
  
    f.mat <- matrix(0, nrow=n.states*N.s.states, ncol=n.obs)
    scaling.factor <- rep(0,n.obs)
    f.mat[,1] <- (1/n.states/N.s.states) * rep(emission.func(params$obs$N.A[1], params$obs$N[1], states), times=N.s.states)

    scaling.factor[1] <- sum(f.mat[,1])
    f.mat[,1] <- f.mat[,1]/scaling.factor[1]

    for(j in 2:n.obs){
        for(i in 1:(n.states*N.s.states)){

            same.s.state <- rep(0,n.states*N.s.states)
            block <- floor((i-1)/params$grid)
            same.s.state[(params$grid*block+1):(params$grid*block+params$grid)] <- 1

            this.args <- f.mat[,j-1] * transition.func(rep(states, times=N.s.states), states[(i-1)%%n.states+1], j, params, same.s.state, rep(params$s, each=length(params$states)))
            this.elem <- emission.func(N.A[j], N[j], states[(i-1)%%n.states+1])*sum(this.args)
            f.mat[i,j] <- this.elem
        }
        scaling.factor[j] <- sum(f.mat[,j])
        f.mat[,j] <-   f.mat[,j]/scaling.factor[j]
    }

    return(list(scaling.factor=scaling.factor, forward.matrix=f.mat, P=sum(f.mat[,n.obs])))
}

wfhmm_s2d.forward.backwards <- function(states, observations, emission.func, transition.func, params){
  fwd <- wfhmm_s2d.forward(states, observations, emission.func, transition.func, params)
  bwd <- wfhmm_s2d.backwards( states, observations, emission.func, transition.func, params, scaling.factor=fwd$scaling.factor)
  fwd.bwd <- fwd$forward.matrix*bwd$backwards.matrix/fwd$P
  return(list(f.mat=fwd$forward.matrix, b.mat=bwd$backwards.matrix, fb.mat=fwd.bwd, scaling.factor=fwd$scaling.factor, P=fwd$P))
}

wfhmm_s2d.backwards <- function(states, observations, emission.func, transition.func, params, scaling.factor=NULL){
    n.obs <- length(observations)
    N.s.states <- length(params$s) 
    n.states <- length(states)
    if(all(is.null(scaling.factor))){scaling.factor <- rep(1,n.obs)}

    N.A <- params$obs$N.A
    N <- params$obs$N
  
    b.mat <- matrix(0, nrow=n.states*N.s.states, ncol=n.obs)
    b.mat[,n.obs] <- 1
    for(j in rev(1:(n.obs-1))){
        this.em <- emission.func(N.A[j+1], N[j+1], states)
        for(i in 1:(n.states*N.s.states)){

            same.s.state <- rep(0,n.states*N.s.states)
            block <- floor((i-1)/params$grid)
            same.s.state[(params$grid*block+1):(params$grid*block+params$grid)] <- 1
            
            this.args <- b.mat[,j+1] * transition.func(states[(i-1)%%n.states+1], rep(states, times=N.s.states), j, params, same.s.state, params$s[block+1]) * rep(this.em, times=N.s.states)
            b.mat[i,j] <- sum(this.args)/scaling.factor[j+1]
        }
    }

    return(list(backwards.matrix=b.mat))
}


wfhmm_s2d.transition <- function( f.from, f.to, t, params, same, s){
    s <- rep(params$s, each=length(params$states))
    if(length(params$Ne)==1){
        Ne <- params$Ne
    }else{
        Ne <- params$Ne[t]
    }
    int <- params$interval/2
    h <- params$h
    mu <- f.from+2*s*f.from*(1-f.from)*(f.from+h*(1-2*f.from))
    si <- sqrt(f.from*(1-f.from)/Ne)
    probs <- pnorm(f.to+int,mu,si)-pnorm(f.to-int,mu,si)
    all.probs <- probs*ifelse(same==1, 1-params$s.trans.p, params$s.trans.p)
    return(all.probs)
}

## The expected value of f_{t+1}, given f_t=i
## i is the index in the full table, not just the frequency table
wfhmm_s2d.conditional.expectation <- function(call, i, t){
  fwd <- call$fb$f.mat
  bwd <- call$fb$b.mat
  fb <- call$fb$fb.mat
  N.s.states <- length(call$params$s) 
  n.states <- length(call$params$states)
      
  params <- call$params
  tf <- call$transition.func
  ef <- call$emission.func
  states <- params$states

  N.A <- call$params$obs$N.A
  N <- call$params$obs$N
  
  same.s.state <- rep(0,n.states*N.s.states)
  block <- floor((i-1)/params$grid)
  same.s.state[(params$grid*block+1):(params$grid*block+params$grid)] <- 1

  probs <- fwd[i,t]*tf(states[(i-1)%%n.states+1], rep(states, times=N.s.states), t, params, same.s.state, params$s[block+1])*bwd[,t+1]/fb[i,t]

  this.em <- ef(N.A[t+1], N[t+1], states)
  
  probs <- probs*rep(this.em, times=N.s.states)
  probs <- probs/sum(probs)             
  return(sum(probs*rep(states, times=N.s.states)))
}
