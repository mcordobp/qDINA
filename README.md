# qDINA

Code and supporting materials for the qDINA model.

The qDINA model extends the DINA model by introducing item-specific dispersion parameters, \(\phi_j\), which modify the contribution of each item through power weights \(w_j = 1/\phi_j\).

Repository: https://github.com/mcordobp/qDINA

## Repository structure

- `R/phi_calibration.R`  
  Calibration of the item-specific dispersion parameters.

- `scr/qdina_gibbs.cpp`  
  Gibbs sampler for the qDINA model.

- `simulation/generate_simulation_data.R`  
  Generates the simulated datasets used in the simulation study.

- `simulation/fit_qdina_simulation.R`  
  Fits the baseline DINA model, calibrates the dispersion parameters, and fits the qDINA model for a simulation condition.

- `empirical/fit_qdina_empirical.R`  
  Fits the DINA and qDINA models to the empirical dataset.

- `input/Simulation5/`  
  Simulated datasets used in the simulation study.

## R dependencies

The code requires the following R packages:

```r
install.packages(c(
  "Rcpp",
  "RcppArmadillo",
  "rgen",
  "dina",
  "simcdm",
  "CDM",
  "dplyr",
  "tidyr"
))
```

## Simulation study

The simulated datasets can be regenerated with:

```r
source("simulation/generate_simulation_data.R")
```

To fit one simulation condition, edit the input and output paths in:

```text
simulation/fit_qdina_simulation.R
```

and run the script from the root directory of the repository.

The fitting procedure follows three steps:

1. Fit the baseline DINA model.
2. Estimate the item-specific dispersion parameters.
3. Fit the qDINA model using the calibrated \(\phi_j\) values.

## Empirical application

The empirical application uses the `data.cdm05` dataset available in the `CDM` R package.

Run:

```r
source("empirical/fit_qdina_empirical.R")
```

from the root directory of the repository.

## References

Culpepper, S. A. (2015). Bayesian estimation of the DINA model with Gibbs sampling. *Journal of Educational and Behavioral Statistics, 40*(5), 454–476. https://doi.org/10.3102/1076998615595403

Culpepper, S. A., & Balamuta, J. J. (2019). *dina: Bayesian estimation of DINA model*. R package. https://cran.r-project.org/package=dina
