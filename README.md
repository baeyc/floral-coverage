# Bayesian estimation of floral coverage for different landuse covers

R code used for the following paper:

A model to account for data dependency when estimating floral cover in different land use types over a season. Baey, Sahlin, Clough, Smith. Environmental and Ecological Statistics, 24(4):505-527, 2017. [doi:10.1007/s10651-017-0387-x](https://doi.org/10.1007/s10651-017-0387-x)


## Quick start
There are two ways to run the code: either using a (bash-like) script or directly in R. 

To run it via a script, uncomment the following lines:

```r
options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
```

and comment the following line if it appears:
```r
#options(echo=FALSE)
args <- c("flowStrip","gumbel","1","0","0.7","0","0.7","150000","1000")

```

To run it in R, do the reverse (comment the first two lines above and explicitly assign values to `args`)

### Arguments
`args` contains all the arguments that will be used, in the following order (it has to be respected):
 - landuseCat: the landuse category, which can be one of `"flowStrip"`, 
 - copFamily: the Archimedean copula family, which can be one of `"clayton"`, `"joe"`, `"gumbel"` or `"frank"`
 - nbrep: repetition number (to be used to name the output files)
 - m1l: the lower bound of the uniform prior on $\mu_1$
 - m1u: the upper bound of the uniform prior on $\mu_1$
 - m2l: the lower bound of the uniform prior on $\mu_2$
 - m2u: the upper bound of the uniform prior on $\mu_2$
 - nb_mcmc: the size of the chain (number of iterations after the burn-in period)
 - nb_burn: the number of burn-in iterations
