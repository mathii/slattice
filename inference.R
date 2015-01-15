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

## Inference in the one population model, both exact and using the hmm
## includes general wrapper functions.

source("wfhmm.R")

## MLE for first order approximation...
## with haploid (or genic selection)
approx.MLE.haploid <- function(f){
  g=length(f)
  fi <- f[1:(g-1)]

  term1 <- f[g]-f[1]
  term2 <- sum(fi*(1-fi))

  if(term2){
    return(term1/term2)
  } else{
    return(0)                           #This assumes that if f(1-f)=0 then s=0, if f switched between 0 and 1 it would be wrong
  }                                     #but we are assuming that if f(1-f)=0 then f is identially 0 or 1. 
}

## EM update step, fb is the matrix of posterior probabilities,
## and states are the corresponding frequency values. fwd is the
## forward matrix, which is only used if h!=0.5 

approx.em.MLE <- function(call, h=0.5){
  ## if(h!=0.5){stop("not implemented")}
  fb <- call$fb$fb.mat
  states <- call$params$states

  g=dim(fb)[2]
  D=length(states)

  ht <- states+h*(1-2*states)
  EfT <- h*sum(fb[,g]*states)             #expected final f
  Ef0 <- h*sum(fb[,1]*states)             #expected initial f
  het <- sum(apply(fb[,1:(g-1)], 2, '*', ht*ht*states*(1-states))) #expected heterozygosity

  ## only compute this expensive term if we need to. 
  fwd.term <- 0
  if(0.5!=h){
    this.sum <- matrix(0,nrow=D, ncol=g-1)
    for(i in 1:D){
      for(j in 1:(g-1)){
        if(0==fb[i,j]){
          this.sum[i,j] <- 0
        }else{
          this.sum[i,j] <- fb[i,j]*states[i]*(sum(wfhmm.conditional.expectation(call, i, j)-states[i]))
        }
      }
      ## sometimes some of the conditional exp terms can underflow and become NaN. Ignore those there (they are 0)
      fwd.term <- sum(this.sum[!is.nan(this.sum)])*(1-2*h)
    }
  }
  
  return((EfT-Ef0+fwd.term)/het/2)
}
  

## approximate MLE for the general dominance case
## should agree with the haploid MLE when h=0.5

approx.MLE <- function(f, h=0.5){
  g=length(f)
  fi <- f[1:(g-1)]                      #f_i
  fip <- f[2:g]                         #f_{i+1}
  hi <- fi+h*(1-2*fi)

  term1 <- sum(hi*(fip-fi))
  term2 <- 2*sum(hi*hi*fi*(1-fi))

  if(term2){
    return(term1/term2)
  } else{
    return(0)                           #This assumes that if f(1-f)=0 then s=0, if f switched between 0 and 1 it would be wrong
  }                                     #but we are assuming that if f(1-f)=0 then f is identially 0 or 1. 
}

## MLE for second order approximation
## with haploid (or genic selection)

approx.second.order.MLE <- function(f){
  g=length(f)
  fi <- f[1:(g-1)]
  k1 <- sum(fi*(1-fi))
  k2 <- sum(fi*fi*(1-fi))
  df <- f[g]-f[1]

  s.hat <- (k1-sqrt(k1*k1-4*k2*df))/(2*k2)
  return(s.hat)
}

## derivative of the log-likelihood

dl.ds <- function(s, N, f){
  g=length(f)
  ni <- f[2:g]*N
  fi <- f[1:(g-1)]

  term1 <- sum(ni*fi/(fi+s*fi))
  term2 <- sum(N*fi/(1+s*fi))

  return(term1-term2)
}

## derivative of the likelihood in the form of equation 2

dl.ds.eqn2 <- function(s, f, multiply.through=FALSE){
  g=length(f)
  term1 <- sum(f[2:g])
  term2 <- (1+s)*sum(1/(1/f[1:(g-1)]+s))
  term3 <- sum((f[1:(g-1)]*(1+s))/(1+f[1:(g-1)]*s))
  if(multiply.through){
    return((term1-term3)*prod((1+f[1:(g-1)]*s)))
  }else{
    return(term2-term1)
  }
}

## calculate the likelihood of s, given the total population frequency and
## history

s.likelihood <- function(s, N, fr){
  g=length(fr)

  ni <- fr[2:g]*N
  fi <- fr[1:(g-1)]

  nonzero.terms <- ni!=0                #since logs go to infinity. 
  
  term1 <- ni*log((fi+s*fi)/(1+s*fi))
  term1 <- sum(term1[nonzero.terms])
  term2 <- sum((N-ni)*log(1+s*fi))

  return(term1-term2)
}

## just wrap the simple estimator with linear interpolation
## Most of the arguments are unused but just for the sake of it.

s.estimate.simple <- function(obs, Ne, h=0.5, tol=NA, max.iters=NA, verbose=FALSE, params=list()){
  g <- NROW(obs)
  f.obs <- obs$N.A/obs$N
  ## Linearly interpolate missing values. 
  f <- approx((1:g)[!is.na(f.obs)], f.obs[!is.na(f.obs)], 1:g, rule=2)$y
  s <- approx.MLE(f, h=h)
  return(list(s=s))
}

## EM estimator for s - alternate computing the MLE of s with
## computing the most likely path, given s. Admittedly this isn't
## usually what people mean when they think of an em estimator

s.estimate.hard.em <- function( obs, Ne, h=0.5, tol=0.001, max.iters=100, verbose=FALSE, params=list() ){
  if(is.null(params$transitions)){params$transitions <- "normal"}
  if(is.null(params$grid)){params$grid <- 100}

  ## Get an inital estimate of s from the observed frequencies.
  s <- s.estimate.simple(obs, Ne, h=h)$s
  s.old <- s+2*tol
 
  if(verbose){
    cat(paste("Initial s =", s, "\n"))
    plot.wright.fisher.observations(obs, main=paste("Iteration", 0,  "s =", round(s,4),"\n"))
  }
  
  
  iter=1
  while(abs(s-s.old)>tol & iter<max.iters){
    f <- wfhmm.call(obs, s, Ne, h=h, grid=params$grid, transitions=params$transitions)$viterbi
    s.old <- s
    s <- approx.MLE(f, h=h)
    
    if(verbose){
      plot.wright.fisher.observations(obs, main=paste("Iteration", iter,  "s =", round(s,4) ,"\n"))
      lines(f, col="#CC5500")
      cat(paste("Iteration", iter,  "s =", s, "\n"))
    }
    iter <- iter+1
  }

  return(list(s=s, path=f, iterations=iter-1))
}

## Estimate s using a soft em estimator - so the estimator is averaged
## over the posterior distribution of frequencies.

s.estimate.soft.em <-  function( obs, Ne, h=0.5, tol=0.001, max.iters=100, verbose=FALSE, params=list() ){
  if(is.null(params$transitions)){params$transitions <- "normal"}
  if(is.null(params$grid)){params$grid <- 100}

  ## Get an inital estimate of s from the observed frequencies.
  s <- s.estimate.simple(obs, Ne, h=h)$s
  s.old <- s+2*tol
 
  iter.values <- c(s)
  if(verbose){
    cat(paste("Initial s =", s, "\n"))
    plot.wright.fisher.observations(obs, main=paste("Iteration", 0,  "s =", round(s,4),"\n"))
  }
  
  
  iter=1
  while(abs(s-s.old)>tol & iter<max.iters){
    s.old <- s
    call <- wfhmm.call(obs, s, Ne, h=h, grid=params$grid, transitions=params$transitions, viterbi=FALSE, forward.backward=TRUE)
    states <- call$params$states
    s <- approx.em.MLE(call, h=h)
    iter.values <- c(iter.values, s)

    if(verbose){
      require("RColorBrewer")
      fb <- call$fb$fb.mat
      plot.wright.fisher.observations(obs, main=paste("Iteration", iter,  "s =", round(s,4),"\n"))
      contour(1:dim(fb)[2], states, t(fb), add=TRUE, col=brewer.pal(9, "Reds"), nlevels=9)
      cat(paste("Iteration", iter,  "s =", s,"\n"))
    }
    iter <- iter+1
  }

  return(list(s=s, posterior=call$fb$fb.mat, iterations=iter-1, iter.values=iter.values, call=call))
}

## stochastic EM estimator for s - start with an estimate of s,
## then simulate [n.paths] paths, compuse s hat for each of these
## paths and take the average as your next estimate of s. continue
## Until your estimate of s converges

s.estimate.stochastic.em <- function( obs, Ne, h=0.5, tol=0.001, max.iters=100, verbose=FALSE, params=list(n.paths=100)){
  ## Slightly hacky default value
  if(is.null(params$n.paths)){params$n.paths <- 100}
  if(is.null(params$transitions)){params$transitions <- "normal"}
  if(is.null(params$grid)){params$grid <- 100}

  
    ## Get an inital estimate of s from the observed frequencies.
  s <- s.estimate.simple(obs, Ne, h=h)$s
  s.old <- s+2*tol
 
  if(verbose){
    cat(paste("Initial s =", s, "\n"))
    plot.wright.fisher.observations(obs)
  }
  
  
  iter=1
  while(abs(s-s.old)>tol & iter<max.iters){
    paths <- wfhmm.call(obs, s, Ne, , h=h, viterbi=FALSE, paths=params$n.paths, grid=params$grid, transitions=params$transitions)$paths
    s.old <- s
    s.est <- apply(paths, 1, approx.MLE, h=h)
    s <- mean(s.est, na.rm=TRUE)        #In case one of the paths blows up... hope it's not too bad though!

    if(verbose){
      plot.wright.fisher.observations(obs, main=paste("Iteration", iter))
      for(i in 1:params$n.paths){lines(paths[i,], col="lightgrey")}
      cat(paste("Iteration", iter,  "s =", s, " (", sd(s.est),")\n"))
    }
    iter <- iter+1
  }

  return(list(s=s, paths=paths, iterations=iter-1))
}

## Estimate s using one of the specified methods, and then generate
## an estimate of the standard deviation of the estimate by simulating
## paths.
## Arguments:
## obs - a list with two vectors names N and N.A with the total number of observations and A alleles respectively
## Ne - Effective population size
## h - dominance coefficeint. h=0.5 is genic selection
## tol - tolerance: stop when subsequent changes are less than this value
## max.iters - stop after this many iterations.
## viterbi - return the most likely path

estimate.s <- function( obs, Ne, h=0.5, tol=0.001, max.iters=100, viterbi=TRUE, verbose=FALSE, method=c( "Soft EM", "Hard EM", "Stochastic EM", "Simple" ), params=list(), likelihood="Model"){

  if(!all(names(obs)==c("N", "N.A"))){stop("Observations must have 2 columns named N and N.A")}
  
  func.map <- list("Hard EM"=s.estimate.hard.em, "Stochastic EM"=s.estimate.stochastic.em, "Simple"=s.estimate.simple, "Soft EM"=s.estimate.soft.em)
  s.func <- func.map[method[1]][[1]]    #For need to get the first element of the extracted list to get the actual function
  est <- s.func( obs, Ne, h=h, tol=tol, max.iters=max.iters, verbose=verbose, params=params)

  sample.paths <- wfhmm.call(obs, est$s, Ne, h=h, viterbi=viterbi, likelihood=likelihood)
  est$params <- sample.paths$params
  est$fb <- sample.paths$fb
  if(viterbi){
    est$viterbi <- sample.paths$viterbi
  }
  est$log.likelihood <- sample.paths$log.likelihood
  est$path.log.likelihood <- sample.paths$path.log.likelihood
  est$obs.log.likelihood <- sample.paths$obs.log.likelihood

  return(est)
}

## Find confidence interval - compute a confidence interval by solving to find where
## the likelihood is less that pchisq(alpha) lower than the ML. Twice the difference in
## ll has a chisq 1 distribution

find.confidence.interval <- function(obs, Ne, s.hat, alpha=0.05 , h=0.5, lower.range=s.hat-0.1, upper.range=s.hat+0.1, tol=1e-4){
  diff <- 0.5*qchisq(alpha,1,lower.tail=FALSE)
  max.ll <- wfhmm.call(obs, s.hat, Ne, h=h, viterbi=TRUE, paths=0, likelihood="Model")$log.likelihood
  p <- list(obs=obs, Ne=Ne, diff.ll=diff, max.ll=max.ll)

  lower.ci <- uniroot(ci.objective, c(lower.range, s.hat), p, tol=tol, extendInt="downX")
  upper.ci <- uniroot(ci.objective, c(s.hat, upper.range), p, tol=tol, extendInt="upX")

  return(c(lower.ci$root, upper.ci$root))
}

## objective function for the likelihood, using the viterbi path

ci.objective <- function(s, p){
  solve <- wfhmm.call(p$obs, s, p$Ne, viterbi=FALSE, forward.backward=TRUE, paths=0, likelihood="Model")
  return(p$max.ll-solve$log.likelihood-p$diff.ll)
}

## p-value: approximate two tailed p-value using the log-likelihood limiting
## distribution, testing against H_0: s=0

p.value <- function(obs, Ne, s.hat, h=0.5){
  max.ll <- wfhmm.call(obs, s.hat, Ne, h=h, viterbi=FALSE, paths=0, likelihood="Model")$log.likelihood
  null.ll <- wfhmm.call(obs, 0, Ne, viterbi=FALSE, paths=0, likelihood="Model")$log.likelihood
  p.val <- 1-pchisq(2*(max.ll-null.ll),1)
  return(p.val)
}
