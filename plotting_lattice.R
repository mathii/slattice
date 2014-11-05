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

## plotting in the lattice case

source("plotting.R")
source("lattice.R")
source("inference_lattice.R")

## Make a cute plot of the allele frequencies for the WF lattice model

plot.wright.fisher.lattice.afs <- function(n, f, s=0, scale=c(-0.1,0.1), ...){
  k1=dim(f)[1]
  k2=dim(f)[2]
  g=dim(f)[3]
  if(!is.matrix(s)==1){s=matrix(s, nrow=k1, ncol=k2)}

  pal<-brewer.pal(11, "RdYlGn")
  pal.seq <- seq(-0.1,0.1,length.out=11+1)
  
  ## draw grid
  plot(0,0, col="white", xlim=c(0,1), ylim=c(-0.16,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)

  for(i in 1:k1){
    for(j in 1:k2){
      rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1), col=pal[findInterval(s[i,j], pal.seq, all.inside=TRUE)])
    }
  }
  
  for(i in 1:k1){
    for(j in 1:k2){
      slice=f[i,j,]
      xpos=(j-1)/k2+(0:(g-1))/((g-1)*k2)
      ypos=(1-i/k1)+slice/k1
      lines(xpos, ypos, col= "black", lwd=2, type="s")
    }
  }

  for(i in 1:11){
    rect((i-1)/11, -0.2, i/11, -0.1, col=pal[i])
  }
  mtext(scale, 1, at=c(0,1))
  mtext("Selection coefficient", 1, at=c(0.5))
}

## Plot the profile likelihoods around the MLE
## if m is specified, use that, else find the MLE of
## m as well. 
plot.lattice.likelihoods <- function(fr, N, m=NULL, true.s=NULL, true.m=NULL, extend=0.1, grid=100, ... ){
  k1 <- dim(fr)[1]
  k2 <- dim(fr)[2]
  if(!is.matrix(true.s)){true.s <- matrix(true.s,k1,k2)}
  
  pal<-brewer.pal(11, "RdYlGn")
  pal.seq <- seq(-0.1,0.1,length.out=11+1)

  cat("Finding MLE\n\n")
  if(is.null(m)){
    full.result <- optim(c(rep(0,k1*k2),0.01*N), s.m.lattice.likelihood.vector, N=N, fr=fr, method="BFGS")
    result <- full.result$par
    maximised.l <- -full.result$value
    s.hat <- matrix(result[1:(k1*k2)], nrow=k1, byrow=TRUE)
    m.hat <- result[k1*k2+1]
  }
  else{
    full.result <- optim(rep(0,k1*k2), s.lattice.likelihood.vector, m=m, N=N, fr=fr, method="BFGS")
    result <- full.result$par
    maximised.l <- -full.result$value
    s.hat <- matrix(result, nrow=k1, byrow=TRUE)
    m.hat <- m
  }
  approx.s.hat <- approx.lattice.s.hat(fr, m.hat/N)
  approx.m.hat <- approx.lattice.m.hat(fr, s.hat)*N

  
  cat(paste("m.hat=", m.hat, "\n", sep=""))
  cat("s.hat=\n")
  print(s.hat)
  cat("approx.s.hat=\n")
  print(approx.s.hat)
  
  s.min <- min(c(min(true.s),min(s.hat),min(approx.s.hat)))-extend*(max(s.hat)-min(s.hat))
  s.max <- max(c(max(true.s),max(s.hat),min(approx.s.hat)))+extend*(max(s.hat)-min(s.hat))
  s.x <- seq(s.min, s.max, length.out=grid)
  
  ll <- array(0, c(k1,k2,grid) )
  for(i in 1:k1){
    for(j in 1:k2){
      for(l in 1:grid){
        s.tmp <- s.hat
        s.tmp[i,j] <- s.x[l]
        ll[i,j,l] <- s.m.lattice.likelihood(s.tmp, m.hat, N, fr)
      }
    }
  }
  
  l.max <- max(ll)+extend*(max(ll)-min(ll))
  l.min <- min(ll)-extend*(max(ll)-min(ll))

  xpos <- (s.x-s.min)/(s.max-s.min)/k2
  ypos <- (ll-l.min)/(l.max-l.min)/k1
  
  ## draw grid
  plot(0,0, col="white", xlim=c(0,1), ylim=c(-0.16,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  
  for(i in 1:k1){
    for(j in 1:k2){
      rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1))

      ## This will go wrong if the ci is not a single interval.
      ci <- (maximised.l-ll[i,j,]<1.92)
      lines(xpos+(j-1)/k2, ypos[i,j,]+(1-i/k1), col="#377EBA", lwd=2)
      lines(xpos[ci]+(j-1)/k2, ypos[i,j,ci]+(1-i/k1), col="#E41A1C", lwd=2)
      
      if(!all(is.null(true.s))){
##         real.col <- pal[findInterval(true.s[i,j], pal.seq, all.inside=TRUE)]
        real.col <- "#4D8F4A"
        lines(rep((j-1)/k2+(true.s[i,j]-s.min)/(s.max-s.min)/k2,2), c(1-i/k1, 1-(i-1)/k1), col=real.col, lty=2)
      }
      lines(rep((j-1)/k2+(s.hat[i,j]-s.min)/(s.max-s.min)/k2,2), c(1-i/k1, 1-(i-1)/k1), col="blue", lty=2)
      lines(rep((j-1)/k2+(approx.s.hat[i,j]-s.min)/(s.max-s.min)/k2,2), c(1-i/k1, 1-(i-1)/k1), col="#E41A1C", lty=2)
    }
  }

  ## Scale for s on the bottom left panel
  text(0,0,round(s.min,3),pos=1)
  text(1/k2,0,round(s.max,3),pos=1)
  mtext(format(l.max,digits=3,scientific=TRUE), 2, at=0)
  mtext(format(l.max,digits=3,scientific=TRUE), 2, at=1/k1)
  
  ## Draw a little panel for m
  m.range <- seq(round(m.hat*0.8,1), round(m.hat*1.2,1), length.out=grid)
  l.m <- sapply(m.range, m.lattice.likelihood, s=s.hat, N=N, fr=fr)
  rect(0.5, -0.2, 1,-0.05)
  m.x <- 0.5+0.5*(m.range-min(m.range))/(max(m.range)-min(m.range))
  l.min <- min(l.m)
  l.max <- max(l.m)+extend*(l.max-l.min)
  m.y <- -0.2+0.15*(l.m-l.min)/(l.max-l.min)

  lines(m.x, m.y, col="#377EBA", lwd=2)
  ci <- (maximised.l-l.m<1.92)
  lines(m.x[ci], m.y[ci], col="#E41A1C", lwd=2)

  lines(rep(0.5+0.5*(m.hat-min(m.range))/(max(m.range)-min(m.range)),2),c(-0.2,-0.05), col="blue", lty=2)
  lines(rep(0.5+0.5*(approx.m.hat-min(m.range))/(max(m.range)-min(m.range)),2),c(-0.2,-0.05), col="#E41A1C", lty=2)

  if(!is.null(true.m)){
    lines(rep(0.5+0.5*(true.m-min(m.range))/(max(m.range)-min(m.range)),2),c(-0.2,-0.05), col=real.col, lty=2)
  }
  
  mtext(c(m.range[1], "m", m.range[grid]), 1, at=c(0.5,0.75,1))

  real.col <- "#4D8F4A"
  ## legend
  legend(0.05,-0.03, c("Likelihood", "95% CI", "MLE", "Approx MLE", "True value" ), col=c( "#377EBA","#E41A1C","blue","#E41A1C",  real.col), lty=c(1,1,2,2,2), lwd=c(2,2,1,1,1), bty="n")
}

## Now plot, on a grid
## f and est should be k1*k2*t matrices. paths is a k1*k2*t*n.paths matrix

plot.wright.fisher.lattice.observations <- function(obs, f=NULL, est=NULL, paths=NULL, true.s=NULL, est.s=NULL, wd.x=0.002, error.bars=FALSE, scale.max=NULL, sd.colour.boxes=10, sd.est.s=NULL, brewerpal="RdYlGn", reversepal=FALSE, palsize=11, lighten=0.5, f.col="#4DAF4A", draw.scale=TRUE, ...){
  k1=dim(obs$N)[1]
  k2=dim(obs$N)[2]
  g=dim(obs$N)[3]
  
  ## draw grid
  ylim <- c(0,1)
  if(!all(is.null(est.s))){ylim <- c(-0.16,1) }
  plot(0,0, col="white", xlim=c(0,1), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)

  ## Palette for plotting estimated selection coefficient. 
  pal<-brewer.pal(palsize, brewerpal)
  n.pals <- length(pal)
  if(reversepal){pal <- rev(pal)}

  if(is.null(scale.max)){
    scale.max <- 0.1
    if(!all(is.null(est.s))){
      scale.max <- max(abs(est.s))
      scale <- c(-scale.max, scale.max)
    }
  }else if(length(scale.max)==2){
    scale <- scale.max
  }else if(length(scale.max)==1){
    scale <- c(-scale.max, scale.max)
  }
  
  pal.seq <- seq(scale[1], scale[2],length.out=n.pals+1)
  
  for(i in 1:k1){
    for(j in 1:k2){

      interval.col <- "#FFFFFF"
      if(!all(is.null(est.s))){         
        if(all(is.null(sd.est.s))){ #If we gave an estimate but no sd. 
          interval.col <- pal[findInterval(est.s[i,j], pal.seq, all.inside=TRUE)]
          interval.col <- inm.lighten.cols(interval.col, lighten)
          rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1), col=interval.col)
        }else{                          #If we gave estimates with sd
          ns <- sd.colour.boxes
          ss <- matrix(rnorm(ns*ns,mean=est.s[i,j], sd=sd.est.s[i,j]),nrow=ns,ncol=ns)
          for(ii in 1:ns){
            for(jj in 1:ns){
              interval.col <- pal[findInterval(ss[ii,jj], pal.seq, all.inside=TRUE)]
              interval.col <- inm.lighten.cols(interval.col, lighten)
              rect((j-1)/k2+(jj-1)/(k2*ns), (1-i/k1)+(ii-1)/(k1*ns), (j-1)/k2+(jj)/(k2*ns), (1-i/k1)+(ii)/(k1*ns), col=interval.col, lwd=0, border=interval.col)
            }
          }
          rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1)) #Finally draw the border
        }
      }else{                         #if we gave no estimate just draw a white box
        rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1), col=interval.col)
      }

      if(!is.null(paths)){
        for(l in 1:dim(paths)[4]){
          slice=paths[i,j,,l]
          xpos=(j-1)/k2+(1:(g))/((g)*k2)
          ypos=(1-i/k1)+slice/k1
          lines(xpos, ypos, col= "grey", lwd=1)
        }
      }
      
      f.o <- obs$N.A[i,j,]/obs$N[i,j,]
      xpos=(j-1)/k2+(1:(g))/((g)*k2)
      ypos=(1-i/k1)+f.o/k1

      
      points(xpos, ypos, col= "#377EBA", pch=16, cex=3/max(k1,k2))
      if(error.bars){
        sdo <- sqrt(f.o*(1-f.o)/obs$N[i,j,])
        real <- !is.na(sdo)
        inm.gaussian.errorbars(xpos[real], ypos[real], sd.y=sdo[real]/k1, y.lim=c(1-i/k1, 1-(i-1)/k1), wd.x=wd.x, white=c(col2rgb(interval.col)) )
        ## Overwrite a bit below with a white box, in case we run over
        if(i==k1){
          rect((j-1)/k2, -0.1, j/k2, 0, col="white", border="white")
        }
        if(j==k2){
         rect(1, (1-i/k1), 1.1, j/k1, col="white", border="white")
        } 
        ## Redraw the borders in case we overwrote them with the error bars
        rect((j-1)/k2, (1-i/k1), j/k2, (1-(i-1)/k1))
      }
      
      if(!all(is.null(f))){
        slice=f[i,j,]
        xpos=(j-1)/k2+(1:(g))/((g)*k2)
        ypos=(1-i/k1)+slice/k1
        lines(xpos, ypos, col= f.col, lwd=2)
      }

      if(!all(is.null(est))){
        slice=est[i,j,]
        xpos=(j-1)/k2+(1:(g))/((g)*k2)
        ypos=(1-i/k1)+slice/k1
        lines(xpos, ypos, col= "#CC5500", lwd=2)
      }


    }
  }

  if(!all(is.null(est.s)) && draw.scale){
    for(i in 1:n.pals){
      rect((i-1)/n.pals, -0.2, i/n.pals, -0.1, col=inm.lighten.cols(pal[i],lighten))
    }
    scale <- round(scale,2)
    mtext(scale, 1, at=c(0,1))
    mtext("Estimated selection coefficient", 1, at=c(0.5))
  }
}
