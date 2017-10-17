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

## Plotting in the  pop case

library(RColorBrewer)

library(plotrix)

########################################################################################################
##
## Smooth gaussian error bars - plots error bars with color intensity that looks like Gaussian
## density. You need to call plot(x,y) first - this just adds the error bars. sd.x and sd.y are
## standard deviations you want plotted, wd.x and wd.y control the width of the bars. 
##
## x<-1:100
## sd<-abs(rnorm(100))
## y<-rnorm(100, mean=0, sd=sd)
## plot(x,y,pch=20, col="red", bty="n")
## inm.gaussian.errorbars(x,y,sd.x=NULL,sd.y=sd, wd.x=0.25)
##
## y.lim doesn't quite do the right thing - squashes the ci in, rather than cutting it off. 
#########################################################################################################

inm.gaussian.errorbars <- function(x, y, sd.x=NULL, sd.y=NULL, col="#377EBA", N=100, sd=3, wd.x=0.025, wd.y=0.025, y.lim=c(-Inf, Inf), white=c(255,255,255), ...)
{
  col.rgb <- col2rgb(col)
  h.1 <- dnorm(0)
  h.0 <- 0
  
  step=white-col.rgb
  sd.pts <- sd*(((-N):N))/N
  rgb.fade.cols <-white - step %*% matrix(dnorm(sd.pts)/h.1,nrow=1)
  str.fade.cols <- mapply(rgb, rgb.fade.cols[1,], rgb.fade.cols[2,], rgb.fade.cols[3,], MoreArgs=list(maxColorValue=255))
  
  if(!all(is.null(sd.x))){
    for(i in 1:length(sd.x)){
      gradient.rect(x[i]-sd*sd.x[i], y[i]-wd.y, x[i]+sd*sd.x[i], y[i]+wd.y,col=str.fade.cols, border=NA, gradient="x")
    }
  }

  if(!all(is.null(sd.y))){
    for(i in 1:length(sd.y)){
      sd.down <- sd
      sd.up <- sd
      if( y[i]-sd.down*sd.y[i] < y.lim[1]){ sd.down <- (y[i]-y.lim[1])/sd.y[i]}
      if( y[i]+sd.up*sd.y[i] > y.lim[2]){ sd.up <- (y.lim[2]-y[i])/sd.y[i]}

      gradient.rect(x[i]-wd.x, y[i]-sd.down*sd.y[i], x[i]+wd.x, y[i]+sd.up*sd.y[i],col=str.fade.cols, border=NA, gradient="y")
    }
  }
}

########################################################################################################
##
## Smooth binomial error bars - plots error bars with color intensity that looks like Binomial
## density. You need to call plot(x,y) first - this just adds the error bars. Vertical error bars
## only, unlike the gaussian error bars. y.a is the sample succeses, and y.n the sample sizes, so
## the colour of the error bar at y is proportional to P(Y=y.a | p=y) where Y.a~B(y.n,p).
## 
#########################################################################################################

inm.binomial.errorbars <- function(x, y.a, y.n, sd.x=NULL, sd.y=NULL, col="#377EBA", N=1000, wd.x=0.025, y.lim=c(0, 1), white=c(255,255,255), N.intensities=100, ...){
    col.rgb <- col2rgb(col)

    range=seq(y.lim[1],y.lim[2],length.out=N)

    step=white-col.rgb

    for(i in 1:length(x)){
      if(0!=y.n[i]){
        p <- y.a/y.n
        intensities <- dbinom(y.a[i], y.n[i], range)
        intensities <- intensities/max(intensities) #Think this is TRTD?
        this.range <- range[intensities>1/256]
        intensities <- intensities[intensities>1/256]
        lo.y <- min(this.range)
        hi.y <- max(this.range)
        intensities <- approx(this.range, intensities, seq(lo.y, hi.y, length.out=N.intensities))$y

        rgb.fade.cols <-white - step %*% matrix(intensities,nrow=1)
        str.fade.cols <- mapply(rgb, rgb.fade.cols[1,], rgb.fade.cols[2,], rgb.fade.cols[3,], MoreArgs=list(maxColorValue=255))
        gradient.rect(x[i]-wd.x, lo.y, x[i]+wd.x, hi.y, col=str.fade.cols, border=NA, gradient="y")      
      }
    }
}

########################################################################################################
##
## Plot contours of a posterior HMM density in such a way that all the points of equal probability are
## connected, i.e. at each line, x% of the probability lies between the two lines.  
## 
#########################################################################################################

inm.plot.isoprobs <- function(probs, states, levels=seq(0.1,0.9,0.1), cols=brewer.pal(9, "Reds"), xoffset=0){
  n <- dim(probs)[2]                    #number of timepoints

  l=1
  for( level in levels){
    a <- level/2
    upper <- rep(0,n)
    lower <- rep(0,n)

    for(i in 1:n){
      these.probs <- cumsum(probs[,i])
      cdf <- approxfun(states, these.probs)
      upperfun <- function(x){return(cdf(x)-a)}
      lowerfun <- function(x){return(cdf(x)-1+a)}
    
      lower[i] <- uniroot(lowerfun, range(states))$root
      upper[i] <- uniroot(upperfun, range(states))$root
    }
    lines(xoffset+1:n, lower, col=cols[l])
    lines(xoffset+1:n, upper, col=cols[l])

    l=l+1
  }
}  

########################################################################################################
## 
## Lighten colours by shifting them towards white
## x=0 gives pure white and x=1 givs 
## 
##########################################################################################################

inm.lighten.cols <- function(cols, x=0.5){
  numbers<-col2rgb(cols)
  numbers <- x*numbers+(1-x)*255
  cols <- mapply(rgb, numbers[1,], numbers[2,], numbers[3,], MoreArgs=list(maxColorValue=255))
  return(cols)
}


## plot the likelihood of s, given the total population size and
## history
plot.likelihood <- function(f, N, slim=c(0,0.1), grid=1000){
  x <- seq(slim[1], slim[2], length.out=grid+1)
  ll <- sapply(x, s.likelihood, N=N, f=f)

  plot(x, ll, bty="n", lwd=2, col="#377EBA", xlab="s", ylab="log-likelihood", type="l")
}

## plot observations and true values of a wright-fisher path, in a cute way.

plot.wright.fisher.observations <- function(obs, f=NULL, est=NULL, true.s=NULL, est.s=NULL, paths=NULL, wd.x=0.1, xoffset=0, N=1000, ...){
  g <- length(obs$N)
  plot(xoffset+1:g, obs$N.A/obs$N, ylab="Frequency", col="#377EBA", pch=16, bty="n", xlab="Generation", ...)

  if(!all(is.null(paths))){
    for(i in 1:NROW(paths)){
      lines(xoffset+1:g, paths[i,], col="grey")
    }
  }
  
  inm.binomial.errorbars((xoffset+1:g), obs$N.A, obs$N, wd.x=wd.x, N=N )
  axis(1)
  
  if(!is.null(f)){
    lines(xoffset+1:g, f, lwd=2, col="#4DAF4A")
  }
  
  if(!is.null(est)){
    lines(xoffset+1:g, est, lwd=2, col="#CC5500")
  }

  ## todo text for s. 
}

## Plot the distribution of various estimators, compared to the real thing
## estimates should be a data frame, with column names being the estimators
plot.estimate.distribution <- function(true.s, estimators, xlab="s", ...){
  n.estimators <- NCOL(estimators)

  cols <- brewer.pal( n.estimators, "Set1")
  
  dens <- density(estimators[,1])
  plot(dens$x, dens$y, col=cols[1], lwd=2, type="l", bty="n", xlab=xlab, ylab="Density", ...)
  abline(v=mean(estimators[,1]), col=cols[1], lwd=2, lty=2)
  if(n.estimators>1){
    for(i in 2:n.estimators){
      dens <- density(estimators[,i])
      lines(dens$x, dens$y, col=cols[i], lwd=2, lty=1)
      abline(v=mean(estimators[,i]), col=cols[i], lwd=2, lty=2)
    }
  }
  abline(v=true.s, col="black", lty=2, lwd=2, ... )

  legend( "topright", c("True s", colnames(estimators)), bty="n", col=c("black", cols), lwd=rep(2,n.estimators+1), lty=c(2,rep(1,n.estimators+1)))
  
}

