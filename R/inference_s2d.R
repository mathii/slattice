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

## Inference in the one population model with time-varying selection
## includes general wrapper functions.

approx.em.MLE.s2d <- function(call, h=0.5){
  if(h!=0.5){stop("Not implemented")}
  states <- call$params$states

  g=dim(call$fb$fb)[2]
  N.s.states <- length(call$params$s) 
  N.states <- length(states)
  
  updated.s <- rep(0, N.s.states)
  
  for(i in 1:N.s.states){
      start <- call$params$grid*(i-1)+1
      end <- call$params$grid*i
      these.fb <- call$fb$fb[start:end,]

      fwd.expectation <- rep(0, g-1)
      for(j in 1:(g-1)){
          this.sum <- 0
          for(k in 1:N.states){
              f1 <- these.fb[k,j]
              f2 <- wfhmm_s2d.conditional.expectation(call, start+k-1, j)
              ## sometimes some of the conditional exp terms can underflow and become NaN.
              ## Ignore those there (they are 0) 

              if(!is.nan(f1*f2)){
                  this.sum <- this.sum+f1*f2
              }
          }
          fwd.expectation[j] <- this.sum
      }

      mean.f <- colSums(apply(these.fb, 2, "*", states))[1:(g-1)]
      
      ## top <- sum(diff(mean.f)*call$fb$s.fb[i,1:(g-1)])
      top <- sum(fwd.expectation-mean.f)
      mean.het <- colSums(apply(these.fb[,1:(g-1)], 2, "*", states*(1-states)))
      ## bot <- sum(mean.het*call$fb$s.fb[i,1:(g-1)])
      bot <- sum(mean.het)
      updated.s[i] <- top/bot
  }
  
  return(updated.s)
}

#Estimator from the vitebi algorithm
approx.MLE.s2d <- function(call, h=0.5){
    if(h!=0.5){stop("Not implemented")}
    N.s.states <- length(call$params$s) 
    updated.s <- rep(0, N.s.states)

    f.path <- call$viterbi
    s.path <- call$s.viterbi
    g=length(f.path)
    
    for(i in 1:N.s.states){
        include <- s.path==i
        f.diff <- diff(f.path)
        top <- sum(f.diff[include[1:(g-1)]])
        bot <- sum((f.path*(1-f.path))[include[1:(g-1)]])

        if(bot>0){
            updated.s[i] <- top/bot
        }else{
            updated.s[i] <- call$params$s[i]
        }
    }
  
  return(updated.s)

    
}

## Estimate s using a soft em estimator - so the estimator is averaged
## over the posterior distribution of frequencies.

s.estimate.soft.em.s2d <-  function( obs, Ne, initial.s.states=0, h=0.5, s.trans.p=1/NROW(obs), tol=0.001, max.iters=100, verbose=FALSE, viterbi=FALSE, params=list() ){
  if(is.null(params$transitions)){params$transitions <- "normal"}
  if(is.null(params$grid)){params$grid <- 100}

  ## Use initial values of s.
  s <- initial.s.states
  s.old <- s+2*tol
 
  if(verbose){
    cat(paste("Initial s =", s, "\n"))
    plot.wright.fisher.observations(obs, main=paste("Iteration", 0,  "s =", round(s,4),"\n"))
  }
  
  likelihoods <- c()
  iter=1
  while(max(abs(s-s.old))>tol & iter<max.iters){
    s.old <- s
    call <- wfhmm_s2d.call(obs, s, Ne, h=h, s.trans.p=s.trans.p, grid=params$grid, transitions=params$transitions, viterbi=verbose|viterbi, forward.backward=TRUE)
    states <- call$params$states
    s <- approx.em.MLE.s2d(call, h=h)

    if(verbose){
      require("RColorBrewer")
      fb <- call$fb$f.fb
      plot.wright.fisher.observations(obs, main=paste("Iteration", iter,  "s =", round(s,4),"\n"), ylim=range(call$params$states))
      contour(1:dim(fb)[2], states, t(fb), add=TRUE, col=brewer.pal(9, "Reds"), nlevels=9)
      nseg <- dim(fb)[2]
      segments(1:(nseg-1), call$viterbi[1:(nseg-1)], 2:nseg, call$viterbi[2:nseg], col=brewer.pal(9, "Set1")[call$s.viterbi])
      cat(paste("Iteration", iter,  "s =", s,"\n"))
    }
    
    iter <- iter+1
    likelihoods <- c(likelihoods, call$log.likelihood)
  }

  call <- wfhmm_s2d.call(obs, s, Ne, h=h, s.trans.p=s.trans.p, grid=params$grid, transitions=params$transitions, viterbi=verbose|viterbi, forward.backward=TRUE)

  return(list(s=s, posterior=call$fb$fb.mat, iterations=iter-1,call=call, ll=likelihoods, call=call))
}

s.estimate.hard.em.s2d <- function( obs, Ne,  initial.s.states=0, h=0.5, tol=0.001, max.iters=100, verbose=FALSE, params=list() ){
  if(is.null(params$transitions)){params$transitions <- "normal"}
  if(is.null(params$grid)){params$grid <- 100}

  ## Get an inital estimate of s from the observed frequencies.
  s <-  initial.s.states
  s.old <- s+2*tol

  N.s.states <- length(call$params$s) 
  
  if(verbose){
    cat(paste("Initial s =", s, "\n"))
    plot.wright.fisher.observations(obs, main=paste("Iteration", 0,  "s =", round(s,4),"\n"))
  }
  
  likelihoods <- c()
  iter=1
  while(abs(min(s-s.old))>tol & iter<max.iters){
    call <- wfhmm_s2d.call(obs, s, Ne, h=h, grid=params$grid, transitions=params$transitions)
    s.old <- s
    s <- approx.MLE.s2d(call, h=h)
    print(s)
    if(verbose){
        plot.wright.fisher.observations(obs, main=paste("Iteration", iter,  "s =", round(s,4) ,"\n"))
        for(i in 1:N.s.states){
            g <- length(call$viterbi)
            segments(1:(g-1), call$viterbi[1:(g-1)], 2:g, call$viterbi[2:g], col=(1+1:N.s.states)[call$s.viterbi])

        }
        cat(paste("Iteration", iter,  "s =", s, "\n"))
    }
    iter <- iter+1
    likelihoods <- c(likelihoods, call$log.likelihood)
  }
  return(list(s=s, posterior=call$fb$fb.mat, iterations=iter-1, call=call, ll=likelihoods))
}


