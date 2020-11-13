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

## Functions for simulating selection in the Wright-Fisher model
## and generating observations from the model.

## Simulate a wf population of fixed size n for g
## generations, with initial allele frequency p0 and
## selection coefficient s

simulate.wright.fisher <- function(n, g, p0, s){
  frequency <- rep(0,g)
  frequency[1] <- p0
  f <- p0
  
  for(i in 2:g){
    p.sel <- f*(1+s)/(1+f*s)
    n.A <- rbinom(1, n, p.sel)
    f <- n.A/n
    frequency[i] <- f
  }

  return(frequency)
}

## Simulate a diploid wf population with dominance
simulate.wright.fisher.dominance <- function(n, g, p0, s, h){
  if(n%%2){stop("n must be divisible by 2")}
  frequency <- rep(0,g)
  frequency[1] <- p0
  f <- p0
  fit <- c(1, 1+2*h*s, 1+2*s)
  
  alleles <- matrix(rbinom(n,1,f),nrow=2, ncol=n/2)
  
  for(i in 2:g){
    gt <- colSums(alleles)
    gt.counts <- c(sum(gt==0), sum(gt==1), sum(gt==2))
    probs <- 2*gt.counts/n
    f <- sum(fit*probs*c(0,1,2))/sum(fit*probs*2)
    alleles <- matrix(rbinom(n,1,f),nrow=2, ncol=n/2)
    frequency[i] <- sum(alleles)/n
  }
  return(frequency)
}

## Simulate wright fisher until frequenct reaches either 0 or 1
## returns fixation (0/1), number of generations, and total heteozygosity. 
simulate.wright.fisher.fixation <- function(n, p0, s){
  f <- p0
  h <- 0
  i <- 0
  
  while(f>0 & f<1){
    p.sel <- f*(1+s)/(1+f*s)
    n.A <- rbinom(1, n, p.sel)
    f <- n.A/n
    i=i+1
    h=h+2*f*(1-f)
  }
  
  return(c(f,i,h))
}

## Simulate wright fisher conditional on fixation at
## frequency +1
simulate.wright.fisher.fixation.1 <- function(n,p0,s){
  result <- c(0,0,0)
  while(result[1]==0){
    result <- simulate.wright.fisher.fixation(n,p0,s)
  }
  return(result)
}

## generate observations from a path, with N observations made at each time-point

generate.observations.from.path <- function(f, N ){
  if(length(N)==1){
    N <- N+f*0
  }else if(length(N)!=length(f)){
    stop("length(N) must be either 1 (constant size), or of length(f)")
  }

  results <- data.frame(N,N.A=0)
  for(i in 1:length(N)){
    results[i,"N.A"] <- rbinom(1, N[i], f[i])
  }
  return(results)
}

## generate a path, and observations, simulating mising data,
## and random sample size. 

generate.observations <- function(Ne, g, p0, s, missing.p=0.1, size.method="binomial", size.params=list(N=100, p=0.5)){ 
  f <- simulate.wright.fisher(Ne, g, p0, s)
  if(size.method=="constant"){
    N <- rep(size.params$N, g)
  }else if(size.method=="binomial"){
    N <- rbinom(g, size.params$N, size.params$p)
  }else{
    stop(paste("Unknown size method", size.method))
  }

  ## Some are missing, i.e. N=0
  N[as.logical(rbinom(g,1,missing.p))] <- 0
  obs <- generate.observations.from.path(f, N)
  return(list(f=f, obs=obs))
}
  
## generate a random code to tag this run, and store the start time
## set the seed to a fixed value

initialize.simulation <- function(seed=12345){
  set.seed(as.integer(Sys.time()))
  code=paste(sample(letters, 6, replace=TRUE), sep="", collapse="")
  start.time <- date()
  nodename <- Sys.info()[4]
  set.seed(seed)
  results=list(seed=seed, code=code, nodename=nodename)
}

## write log of the simulation to some root directory. 

write.simulation.log <- function(root="./", simdata=list()){
  code=simdata$code
  if(is.null(code)){code <- "unknow"}
  logfile <- paste( root, code, "_log.txt", sep="")
  end.time <- date()
  cat(paste("Started:", simdata$start.time, "\n"), file=logfile, append=FALSE)
  cat(paste("Finished:", end.time, "\n"), file=logfile, append=TRUE)
  cat(paste("Machine:", simdata$nodename, "\n"), file=logfile, append=FALSE)
  cat(paste("Seed:", simdata$seed, "\n"), file=logfile, append=TRUE)
  cat(paste("R Version:", R.version.string,"\n"), file=logfile, append=TRUE)
  cat(paste(c( "git log -1:", system("git log -1", intern=TRUE)), collapse="\n"), file=logfile, append=TRUE)
}

## Simulate a wf population of fixed size n for g
## generations, with initial allele frequency p0 and
## selection coefficients s1 and s2 and changepoints 

simulate.wright.fisher.change <- function(n, g, p0, s1, s2, chg.pts){
        if(length(chg.pts)==0){
            return(simulate.wright.fisher(n, g, p0, s1))
        }
        else if(length(chg.pts>0) & max(chg.pts)<g){
            current.generation <- 1
            selcos <- rep(c(s1,s2), length.out=length(chg.pts)+1)
            epochs <- diff(c(0,chg.pts, g))
            frequency <- c()
            ft <- p0
            for(i in 1:length(epochs)){
                this <- simulate.wright.fisher(n, epochs[i], ft , selcos[i])
                frequency <- c(frequency, this)
                ft <- rev(frequency)[1]
            }
            return(frequency)
        }
        else{
            stop("length of chgpts must be >0 and max <g")
        }
}
