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

## Miscellaneous functions - nothing should link this.

## Hope this is correct!
expected.total.heterozygosity <- function(N,s,p){
  term1 <- (1-p)*(1-exp(-4*N*s*(p+1)))
  term2 <- (1+p)*(exp(-4*N*s)-exp(-4*N*s*p))
  term3 <- (1-exp(-4*N*s))*(1-exp(-4*N*s*p))
  return(2*(term1+term2)/(s*term3))
}

## expected total heterozygosity, conditional on starting at
## 0 and hitting fixation

expected.total.heterozygosity.p0 <- function(N,s){
  term1 <- (1+exp(-4*N*s))/(1-exp(-4*N*s))
 term2 <- 1/(2*N*s*s)
  return((2/s)*term1-term2)
}
  

