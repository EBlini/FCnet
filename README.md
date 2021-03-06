
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FCnet

<!-- badges: start -->

<!-- badges: end -->

An R package for the analysis of Functional Connectivity matrices,
lesional maps, or disconnection maps through elastic NETs.

## Installation

The package is currently available through
[GitHub](https://github.com/). The installation requires the R package
`devtools`.

``` r
# install.packages("devtools")
devtools::install_github("EBlini/FCnet")
```

## Scope

The analysis of (Functional Connectivity) neuroimaging data can be
daunting due to the very high dimensionality of the features involved.
In time, several approaches to the problem have been devised. `FCnet`
allows one to easily implement a three steps procedure consisting of:

1)  **Feature reduction**: the functional connectivity matrices (or
    volumes with lesion/disconnection mapping) are first summarized
    through data reduction techniques such as Principal Component
    Analysis or Independent Components Analysis.

2)  **Robust regression**: the reduced matrix of Weights is then entered
    into a robust regression model (with either ridge or LASSO penalty).
    The model is crossvalidated internally by means of Leave-One-Out
    (possibly nested) crossvalidation.

3)  **Back-projection**: models’ coefficients can be back-projected onto
    the original space, in order to rank the most predictive edges of a
    matrix or voxels.

For useful references, see: [Siegel et
al., 2016](https://www.pnas.org/content/113/30/E4367); [Salvalaggio et
al., 2020](https://academic.oup.com/brain/article/143/7/2173/5861020);
[Calesella et
al., 2020](https://link.springer.com/chapter/10.1007%2F978-3-030-59277-6_3).

## Overview

An overview of the package is available
[here](https://eblini.github.io/FCnet/articles/FCnet_overview_of_package.html).
