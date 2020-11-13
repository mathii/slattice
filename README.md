## Estimate selection coefficients from time series data

This code implements the selection coefficient estimators described in "Estimating selection coefficients in spatially structured populations from time series data in spatially structured populations" - Mathieson & McVean; Genetics (2013), and "Estimating time-varying selection coefficients from time series of allele frequency data" (TBD). 

All the code runs in R, and was developed and tested on R version 2.15.1 and 4.0.2.

To install the package, either (1) install directly from github, or (2) clone the repository and install directly. If you don't need to look at the examples, or it causes problems, you can build without the vignette by setting `build_vignettes=FALSE` in option 1 or `R CMD build --no-build-vignettes` in option 2 (building the vignette takes a few minutes and requires some extra packages). You can find the example code that's used in the vignette in vignettes/slattice.Rmd. 

1. 

```
require(devtools)
install_github("mathii/slattice", build_vignettes=TRUE)
```

2.
```
git clone https://github.com/mathii/slattice.git slattice
R CMD build slattice
R CMD INSTALL slattice_1.0.tar.gz
```

For some brief examples of how to run the estimator, look at the vignette (in R): 

```
vignette("slattice")
vignette("timevar")
```

The files in this package are described below. Most of the libraries are in three parts; xxx.R containing code for the single popuation, xxx_lattice.R amd xxx_s2d.R which depend on xxx.R and contain code for the lattice population and time-varying selection, respectively: 

- inference.R, inference_lattice.R and inference_s2d.R 
Wrapper functions to manage the actual inference, and MLEs for the case of complete observations

- simulation.R and simulation_lattice.R
Code to run simple simulations under the WF (lattice) model. 

- plotting.R and plotting_lattice.R
Make various plots of allele frequencies etc. 

- wfhmm.R, wfhmm_lattice.R and wfhmm_s2d.R
Implement the HMMs described in the papers, with generic emission and transition densities. 

- INM 18Dec12 and 13Nov20

