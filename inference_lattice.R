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

## Inference in the lattice model, both exact and using the hmm
## includes general wrapper functions.

source("wfhmm_lattice.R")
source("lattice.R")

EPSILON <- 10e-100                       #Used as a minimum value in some calculations

## Approximate lattice likelihood, as a function of s and m
## Computed assuming that frequency increments are normally
## distributed.

s.m.lattice.likelihood <- function(s, m, N, fr, h=0.5){
  g <- dim(fr)[3]
  k1 <- dim(fr)[1]
  k2 <- dim(fr)[2]

  if(!is.matrix(s)==1){s=matrix(s, nrow=k1, ncol=k2)}
  mrm <- mig.rate.mat(k1,k2)

  mu <- array(0, c(k1,k2,g))
  s.2 <- array(0,c(k1,k2,g))

  hi <- fr[,,t-1]+h*(1-2*fr[,,t-1])
  
  for(t in 2:g){
    ## term from births in that square
    mu[,,t] <- (1-m*mrm/N) * (fr[,,t-1]+2*hi*s*fr[,,t-1]*(1-fr[,,t-1]))
    s.2[,,t] <- fr[,,t-1]*(1-fr[,,t-1])/N

    ## Sorry!
    for(i in 2:k1){
      for( j in 1:k2){
        mu[i,j,t]=mu[i,j,t]+(m/N)*fr[i-1,j,t-1]
      }
    }

    for(i in 1:(k1-1)){
      for( j in 1:k2){
        mu[i,j,t]=mu[i,j,t]+(m/N)*fr[i+1,j,t-1]
      }

    }
    for(i in 1:k1){
      for( j in 2:k2){
        mu[i,j,t]=mu[i,j,t]+(m/N)*fr[i,j-1,t-1]
      }
    }

    for(i in 1:k1){
      for( j in 1:(k2-1)){
        mu[i,j,t]=mu[i,j,t]+(m/N)*fr[i,j+1,t-1]
      }
    }
  }

  nonzero.var <- s.2[,,2:g]>0
##   s.2[!nonzero.var] <- 1e-8             #To avoid warnings - we will remove these later though
  
  terms <- -log(s.2[,,2:g]) - (fr[,,2:g]-mu[,,2:g])*(fr[,,2:g]-mu[,,2:g])/(2*s.2[,,2:g])

  ## We need to exlcude zero variance terms
  log.likelihood <- sum(terms[nonzero.var])
    
  return(log.likelihood)
}

## wrapper to be used for numerical optimization
## m is known (returns negative ll)
s.lattice.likelihood.vector <- function(s, m, N, fr, h=0.5){
  k1 <- dim(fr)[1]
  k2 <- dim(fr)[2]

  s <- matrix(s, nrow=k1, ncol=k2, byrow=TRUE)
  return(-s.m.lattice.likelihood(s, m, N, fr, h=h))
}

## wrapper for numerical optimisation
## m is unknown. The first vector should of length
## 2 or k1*k2+1 with the last element being m
s.m.lattice.likelihood.vector <- function(s.m, N, fr, h=0.5){
  k1 <- dim(fr)[1]
  k2 <- dim(fr)[2]
  m <- s.m[length(s.m)]
  s <- s.m[1:(length(s.m)-1)]
  s <- matrix(s, nrow=k1, ncol=k2, byrow=TRUE)
  return(-s.m.lattice.likelihood(s, m, N, fr, h=h))
}  

## s.lattice.likelihood with the first two args flipped
m.lattice.likelihood <- function(m, s, N, fr, h=0.5){
  return(s.m.lattice.likelihood(s, m, N, fr, h=h))
}

## Numerical MLE for s in the lattice
lattice.s.hat <- function(fr, m, N, h=0.5){
  k1 <- dim(fr)[1]
  k2 <- dim(fr)[2]
  result <- optim(rep(0.05,k1*k2), s.lattice.likelihood.vector, m=m, N=N, fr=fr, h=h)$par
  s.hat <- matrix(result, nrow=k1, byrow=TRUE)
  return(s.hat)
}  

## approximate MLE for s, given M
## in a lattice model (M=m/N)

approx.lattice.s.hat <- function(f, M, h=0.5){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]

  s.hat <- matrix(0,nrow=k1,ncol=k2)

  for(i in 1:k1){
    for(j in 1:k2){
      s.hat[i,j] <- approx.lattice.s.hat.ij(i,j,f,M, h=h)
    }
  }
  return(s.hat)
}

approx.lattice.s.hat.ij <- function(i,j,f,M, h=0.5){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]

  G <- dim(f)[3]-1
  mrm <- mig.rate.mat(k1, k2)

  fi <- f[i,j,1:G]
  hi <- fi+h*(1-2*fi)

  term1 <- approx.MLE(f[i,j,], h=h)
  term2 <- mrm[i,j]*sum(hi*f[i,j,1:G])
  if(i<k1){term2=term2-sum(hi*f[i+1,j,1:G])}
  if(i>1){term2=term2-sum(hi*f[i-1,j,1:G])}
  if(j<k2){term2=term2-sum(hi*f[i,j+1,1:G])}
  if(j>1){term2=term2-sum(hi*f[i,j-1,1:G])}
  
  het <- sum(2*hi*hi*f[i,j,1:G]*(1-f[i,j,1:G]))
  
  if(het){
    return((term1+M*term2/het))
  } else{
    return(0)
  }
  
}

## MLE for s, assuming that it is constant
approx.lattice.s.hat.constant <- function(f, M, h=0.5){
##   if(h!=0.5){stop("non-genic selection not implemented for constant s")}
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  G <- dim(f)[3]-1
  ht <- f+h*(1-2*f)
  
  mrm <- mig.rate.mat(k1, k2)

  term1 <- matrix(0,k1,k2)                #1-m|k|
  term2 <- matrix(0,k1,k2)                #f_G - f_0 + sum(m|k|fij - sum(fi'j'))
  term3 <- matrix(0,k1,k2)                #f(1-f)

  for( i in 1:k1){
    for(j in 1:k2){
      term2[i,j] <- h*(f[i,j,G+1]-f[i,j,1])
      term3[i,j] <- sum(ht[i,j,1:G]*ht[i,j,1:G]*f[i,j,1:G]*(1-f[i,j,1:G]))
      if(h!=0.5){
        term2[i,j] <- term2[i,j]+(1-2*h)*sum(f[i,j,1:G]*(f[i,j,2:(G+1)]-f[i,j,1:G]))
      }
    }
  }

  return(matrix(sum(term2)/sum(term3)/2,k1,k2))
}

## EM MLE for s, assuming that it is constant
approx.lattice.s.hat.constant.em <- function(call, M, h=0.5){
  fb <- call$fb$fb.mat
  k1 <- dim(fb)[1]
  k2 <- dim(fb)[2]
  G <- dim(fb)[4]-1
  D <- dim(fb)[3]
  states <- call$params$states
  ht <- states+h*(1-2*states)
  
  term2 <- matrix(0,k1,k2)                #f_G - f_0 + sum(m|k|fij - sum(fi'j'))
  term3 <- matrix(0,k1,k2)                #f(1-f)

  for( i in 1:k1){
    for(j in 1:k2){
      term2[i,j] <- h*(sum(fb[i,j,,G+1]*states)-sum(fb[i,j,,1]*states))      
      term3[i,j] <-sum(apply(fb[i,j,,1:G], 2, '*', ht*ht*states*(1-states)))
      if(h!=0.5){
        this.sum <- matrix(0,nrow=D, ncol=G)
        for(s in 1:D){
          for(t in 1:(G)){
            if(0==fb[i,j,s,t]){
              this.sum[s,t] <- 0
            }else{
              this.call <- call
              this.call$fb <- list(fb.mat=call$fb$fb.mat[i,j,,], f.mat=call$fb$f.mat[i,j,,], b.mat=call$fb$b.mat[i,j,,])
              this.call$params <- wfhmm.lattice.copy.params(call$params,i,j)
              this.sum[s,t] <- fb[i,j,s,t]*states[s]*(sum(wfhmm.conditional.expectation(this.call, s, t)-states[s]))
            }
          }
          ## sometimes some of the conditional exp terms can underflow and become NaN. Ignore those there (they are 0)
        }
        term2[i,j] <- term2[i,j]+sum(this.sum[!is.nan(this.sum)])*(1-2*h)
      }
    }
  }
  return(matrix(sum(term2)/sum(term3)/2,k1,k2))
}

## approximate MLE for s, given M
## in a lattice model (M=m/N)

approx.lattice.s.hat.em <- function(call, M, h=0.5){
  k1 <- dim(call$fb$fb.mat)[1]
  k2 <- dim(call$fb$fb.mat)[2]
  s.hat <- matrix(0,nrow=k1,ncol=k2)

  for(i in 1:k1){
    for(j in 1:k2){
      s.hat[i,j] <- approx.lattice.s.hat.ij.em(i,j,call,M, h=h)
    }
  }
  return(s.hat)
}

approx.lattice.s.hat.ij.em <- function(i,j,call,M, h=0.5){
  k1 <- dim(call$fb$fb.mat)[1]
  k2 <- dim(call$fb$fb.mat)[2]
  mrm <- mig.rate.mat(k1, k2)
  fb <- call$fb$fb.mat
  G <- dim(fb)[4]-1
  states <- call$params$states

  ht <- states+h*(1-2*states)

  this.call <- call
  this.call$fb <- list(fb.mat=call$fb$fb.mat[i,j,,], f.mat=call$fb$f.mat[i,j,,], b.mat=call$fb$b.mat[i,j,,])
  this.call$params <- wfhmm.lattice.copy.params(call$params,i,j)
  term1 <- approx.em.MLE(this.call, h=h)

  f <- call$params$est.f
  
  term2 <- mrm[i,j]*sum(apply(fb[i,j,,(1:G)], 2, '*', ht*states))

  exp.ht <- colSums(apply(fb[i,j,,(1:G)], 2, '*', ht))
  if(i<k1){term2=term2-sum(exp.ht*f[i+1,j,(1:G)])}
  if(i>1){term2=term2-sum(exp.ht*f[i-1,j,(1:G)])}
  if(j<k2){term2=term2-sum(exp.ht*f[i,j+1,(1:G)])}
  if(j>1){term2=term2-sum(exp.ht*f[i,j-1,(1:G)])}

  het <- sum(apply(fb[i,j,,1:G], 2, '*', ht*ht*states*(1-states))) #expected heterozygosity  

  if(het){
    return((term1+M*term2/het/2))
  } else{
    return(0)
  }
  
}

## approximate MLE for m
approx.lattice.m.hat <- function(f,s, h=0.5){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  G <- dim(f)[3]
  mrm <- mig.rate.mat(k1, k2)
  if(!is.matrix(s)){s <- matrix(s,k1,k2)}

  term1 <- array(0, dim=c(k1,k2,G-1))     #f_{ij}^t-\mu_{ij}^t
  term2 <- array(0, dim=c(k1,k2,G-1))     #|k|\mu_{ij}^t-\sum \mu_{i'j'}
  term3 <- array(0, dim=c(k1,k2,G-1))     #f(1-f)

  fij <- f[,,1:(G-1)]
  hij <- fij+h*(1-2*fij)
  
  het <- f[,,1:(G-1)]*(1-f[,,1:(G-1)])
  term1 <- f[,,2:G]-f[,,1:(G-1)]-array(apply(het,3,'*',s),dim=c(k1,k2,G-1))*hij*2
  term2 <- array(apply(f[,,1:(G-1)],3,'*',mrm),dim=c(k1,k2,G-1))
  for(i in 1:k1){
    for(j in 1:k2){
      if(i<k1){term2[i,j,]=term2[i,j,]-f[i+1,j,1:(G-1)]}
      if(i>1){term2[i,j,]=term2[i,j,]-f[i-1,j,1:(G-1)]}
      if(j<k2){term2[i,j,]=term2[i,j,]-f[i,j+1,1:(G-1)]}
      if(j>1){term2[i,j,]=term2[i,j,]-f[i,j-1,1:(G-1)]}
    }
  }

  non.zero <- het>0
  m.hat <- -sum((term1*term2/het)[non.zero])/sum((term2*term2/het)[non.zero])
  if(m.hat<0){
    cat("M.hat < 0; setting to 0\n")
    m.hat <- 0
  }
  return(m.hat)
}

## The em version of the MLE for m
approx.lattice.m.hat.em <- function(call, h=0.5){
  k1 <- dim(call$fb$fb.mat)[1]
  k2 <- dim(call$fb$fb.mat)[2]
  mrm <- mig.rate.mat(k1, k2)
  fb <- call$fb$fb.mat
  G <- dim(fb)[4]-1
  states <- call$params$states
  ht <- states+h*(1-2*states)
  s <- call$params$s
  
  term1 <- array(0,dim=c(k1,k2,G))
  term2 <- array(0,dim=c(k1,k2,G))

  C <- states*(1-states)
  for(i in 1:k1){
    for(j in 1:k2){
      A <- states+2*h*s[i,j]*states*(1-states)
      for(t in 1:G){
        B <- states*mrm[i,j]
        if(i<k1){B <- B-sum(states*fb[i+1,j,,t])}
        if(i>1){B <- B-sum(states*fb[i-1,j,,t])}
        if(j<k2){B <- B-sum(states*fb[i,j+1,,t])}
        if(j>1){B <- B-sum(states*fb[i,j-1,,t])}

        this.call <- call
        this.call$fb <- list(fb.mat=call$fb$fb.mat[i,j,,], f.mat=call$fb$f.mat[i,j,,], b.mat=call$fb$b.mat[i,j,,])
        this.call$params <- wfhmm.lattice.copy.params(call$params,i,j)

        eft <- 0*states
        for(l in 1:length(states)){
          if(fb[i,j,l,t]>EPSILON){
            eft[l] <- wfhmm.conditional.expectation(this.call,l,t) 
          }
        }
        term1[i,j,t] <- sum((fb[i,j,,t]*(A-eft)*B/C)[C>0])
        term2[i,j,t] <- sum((fb[i,j,,t]*B*B/C)[C>0])
      }                
    }
  }

  m.hat <- sum(term1)/sum(term2)
##   if(m.hat<0){
##     cat("M.hat < 0; setting to 0\n")
##     m.hat <- 0
##   }
  return(m.hat)
}

## if s.free is FALSE, then we assume that s is constant over the range
## Otherwise it is free to vary. 
iterate.approx.mle <- function(f, initial.s=0, h=0.5, max.iters=100, eps=0.001, s.free=TRUE){
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]

  if(!is.matrix(initial.s)){initial.s <- matrix(initial.s,k1,k2)}
  M.hat <- approx.lattice.m.hat(f, initial.s, h=h)
  if(s.free){s.hat <- approx.lattice.s.hat(f, M.hat, h=h)}
  else{s.hat <- approx.lattice.s.hat.constant(f, M.hat, h=h)}
  old.M.hat <- M.hat+2*eps
  old.s.hat <- s.hat+2*eps

  iter=1
  while(abs(old.M.hat-M.hat)>eps | max(abs(s.hat-old.s.hat))>eps){
    old.M.hat <- M.hat
    old.s.hat <- s.hat
    M.hat <- approx.lattice.m.hat(f, s.hat, h=h)
    if(s.free){s.hat <- approx.lattice.s.hat(f, M.hat, h=h)}
    else{s.hat <- approx.lattice.s.hat.constant(f, M.hat, h=h)}
    iter=iter+1
    if(iter>max.iters){
      break
      warning("Max iterations reached in iterate.approx.mle()")
    }
  }

  return(list(s.hat=s.hat, M.hat=M.hat, iter=iter))
}

## Helper function for estimating s and M - if M is fixed then just estimate s
## If s.free=FALSE, keep s constant. 
update.s.m <- function(f, fixed.M, h=0.5, M=NULL, s.free=TRUE){
  if(fixed.M){                          #Computing an intial guess fom M and s
    if(s.free){s <- approx.lattice.s.hat(f, M, h=h)}
    else{s <- approx.lattice.s.hat.constant(f, M, h=h)}
  }
  else{
    mle <- iterate.approx.mle(f, s.free=s.free, h=h)
    s <- mle$s.hat
    M <- mle$M.hat
  }
  return(list(s=s, M=M))
}

## Now estimate s and m, using the em steps. Don't iterate, just make
## one step
update.s.m.em <- function(call, fixed.M, h=0.5, M=NULL, s.free=TRUE){
  if(s.free){
    s <- approx.lattice.s.hat.em(call, M, h=h)
  } else {
    s <- approx.lattice.s.hat.constant.em(call, M, h=h)
  }

  if(!fixed.M){
    M <- approx.lattice.m.hat.em(call, h=h)
  }
  return(list(s=s, M=M))
}


## estimate s, M, and the frequency paths using an EM algorithm
## 1) Use the observed frequencies to get an inital estimate of
## s, M, and the path
## 2) conditional on s and M, and assuming migration from the
## estimated path, compute the Viterbi path
## 3) repeat until s converges.
##
## The M.update argument tells you how to update M at each stage. If it's a constant, then it just
## fixes M to be that value. if it's  "Simple", if linearly interpolates the data points and computes
## the MLE assuming the true observations. If it's "Simple EM" or "Stochastic EM", it will use either the
## witerbi path, or a sample of paths to estimate M, and s. Also, if you are using stochastic EM, it will
## use a sample of paths to estimate the frequency for migration, rather than the viterbi path. Stops when the
## log likelihood increases by less than "tol", or when it decreases, in which case return the best one found so far. 

s.m.estimate.em.lattice <- function( obs, Ne, M=NULL, h=0.5, update="Simple", tol=0.01, max.iters=100, verbose=FALSE, params=list(), grid=100, s.free=TRUE, initial.M=NULL ){
  if(is.null(params$n.paths)){params$n.paths <- 10}
  if(!(update %in% c("Simple", "Hard EM", "Stochastic EM", "Soft EM"))){stop("Unknown update rule")}
  
  fixed.M <- is.numeric(M)
  k1 <- dim(obs$N)[1]
  k2 <- dim(obs$N)[2]
  g <- dim(obs$N)[3]
  f <- interpolate.f(obs)               #Initial starting point: Linearly interpolate. 
  paths <- NULL

  if(!is.null(initial.M)){
    M <- initial.M
    s.m <- update.s.m(f, TRUE, h=h, M, s.free)
  }
  else{
    s.m <- update.s.m(f, fixed.M, h=h, M, s.free)
  }
  if(update=="Simple"){                  #No iteration here, just return the initial guess. 
    return( list(s=s.m$s, M=s.m$M, f=f, iterations=NA))
  }
  
  s.old <- s.m$s+2*tol
  M.old <- s.m$M+2*tol

  if(verbose){
    cat(paste( "Iteration 0; M =", s.m$M, "\n"))
    cat("s=\n")
    print(s.m$s)
    cat("\n")
  }

  iter=1
  result <- list(log.likelihood=-Inf)
  old.result <- list(log.likelihood=-Inf)
  while((is.infinite(result$log.likelihood) | result$log.likelihood-old.result$log.likelihood>tol | iter < 4) & iter<max.iters){
    old.result <- result
    s.old <- s.m$s
    M.old <- s.m$M

    if(update=="Hard EM"){
      hmm.solved <- wfhmm.lattice.call(obs, s.m$s, s.m$M, Ne, f,  h=h, viterbi=TRUE, grid=grid, forward.backward=FALSE)
      if(verbose){cat(paste("\nlog likelihood =",hmm.solved$log.likelihood , "\n"))}
      s.m <- update.s.m(f, fixed.M,  h=h, M=s.m$M, s.free=s.free)
      f <- hmm.solved$viterbi
    }
    else if(update=="Soft EM"){
      hmm.solved <- wfhmm.lattice.call(obs, s.m$s, s.m$M, Ne, f,  h=h, viterbi=TRUE, grid=grid, forward.backward=TRUE)
      if(verbose){cat(paste("\nlog likelihood =",hmm.solved$log.likelihood , "\n"))}
      s.m <- update.s.m.em(hmm.solved, fixed.M,  h=h, M=s.m$M, s.free=s.free)
      f <- hmm.solved$viterbi
    }
    else if(update=="Stochastic EM"){
      hmm.solved <- wfhmm.lattice.call(obs, s.m$s, s.m$M, Ne, f, h=h, viterbi=TRUE, paths=params$n.paths, grid=grid)
      if(verbose){cat(paste("\nlog likelihood =",hmm.solved$log.likelihood , "\n"))}
##       f <- hmm.solved$viterbi
      paths <- hmm.solved$paths
      f <- apply(paths,c(1,2,3),mean)
      ## get M hat and s.hat for each path
      M.hats <- rep(0,params$n.paths)
      s.hats <- array(0,dim=c(k1,k2,params$n.paths))
      for(i in 1:params$n.paths){
        this.s.m <- update.s.m(paths[,,,i], fixed.M, h=h, M=M, s.free=s.free)
        M.hats[i] <- this.s.m$M
        s.hats[,,i] <- this.s.m$s
      }
      ## Should be able to do this with apply? 
      this.s <- matrix(0,k1,k2)
      for(ii in 1:k1){
        for(jj in 1:k2){
          this.s[ii,jj] <- mean(s.hats[ii,jj,])
        }
      }
      s.m <- list(s=this.s, M=mean(M.hats))
    }
  
    iter <- iter+1
    result <- list(s=s.m$s, M=s.m$M, f=f, iterations=iter, paths=paths, log.likelihood=hmm.solved$log.likelihood)
    if(verbose){
      cat(paste( "Iteration", iter-1, "; M =", s.m$M, "\n"))
      cat("s=\n")
      print(s.m$s)
      cat(paste("gap =", max(max(abs(s.m$s-s.old)), abs(M.old-s.m$M)), "\n\n"))
    }
  }

  return(result)
}

## Wrap the function above and compute likelihood:

estimate.s.m <- function( obs, Ne, M=NULL, h=0.5, update="Simple", tol=0.001, max.iters=100, verbose=FALSE, params=list(), grid=100, s.free=TRUE, final.ll=TRUE, initial.M=NULL ){
  est <- s.m.estimate.em.lattice( obs, Ne, M=M, h=h, update=update, tol=tol, max.iters=max.iters, verbose=verbose, params=params, grid=grid, s.free=s.free, initial.M=initial.M )

  if(est$M<0){
    cat("M.hat < 0; setting to 0\n")
    est$M <- 0
  }
  
  if(final.ll|verbose){
    hmm.solved <- wfhmm.lattice.call(obs, est$s, est$M, Ne, est$f, h=h, viterbi=TRUE, grid=grid)
    est$log.likelihood <- hmm.solved$log.likelihood
  }
  return(est)
}
