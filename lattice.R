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

## Helper functions for the lattice model
.mrm.store <- list()

## Migration rates for grid size k
mig.rate.mat <- function(k1, k2=k1){
  tag <- paste0(as.character(k1),",",as.character(k2))
  if(!is.null(.mrm.store[[tag]])){
      return(.mrm.store[[tag]])
  }


  mrm <- matrix(4, nrow=k1, ncol=k2)
  mrm[1,] <- mrm[1,]-1
  mrm[k1,] <- mrm[k1,]-1
  mrm[,1] <- mrm[,1]-1
  mrm[,k2] <- mrm[,k2]-1

  tmp <- .mrm.store
  tmp[[tag]] <- mrm
  assign(".mrm.store", tmp, .GlobalEnv)
  
  return(mrm)
}

## Linearly interpolate f in each deme. Assumes there is at least one observation

interpolate.f <- function(obs){
  f <- obs$N.A/obs$N
  k1 <- dim(f)[1]
  k2 <- dim(f)[2]
  g <- dim(f)[3]
  ## Linearly interpolate missing values. - Add: Interpolate in space as well as time.  
  for(i in 1:k1){
    for(j in 1:k2){
      if(sum(!is.na(f[i,j,]))==0){
        f[i,j,] <- 0.5
      }
      if(sum(!is.na(f[i,j,]))==1){
        f[i,j,] <- f[i,j,which(!is.na(f[i,j,]))]
      }
      else{
        f[i,j,] <- approx((1:g)[!is.na(f[i,j,])], f[i,j,!is.na(f[i,j,])], 1:g, rule=2)$y
      }
    }
  }
  return(f)
}

