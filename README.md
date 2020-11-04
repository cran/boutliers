
# boutliers package


## Outlier detection and influence diagnostics for meta-analysis

A R package for implementing outlier detection and influence diagnostics for meta-analysis:

- Studentized residuals by leave-one-out analysis
- Likelihood ratio test using a mean-shifted model
- The VRATIO statistic (relative change of the variance of the overall estimator) by leave-one-out analysis
- The TAU2RATIO statistic (relative change of the heterogeneity variance) by leave-one-out analysis

Bootstrap distributions of the influence statistics are calculated by simple codes, and the thresholds to determine outliers are clearly provided.



## Installation

``` r
# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("nomahi/boutliers")
```

Downloads: [please see the release page](https://github.com/nomahi/boutliers/releases)
